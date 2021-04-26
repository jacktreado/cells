/* 

	MAIN FILE FOR JAMMING VIA ISOTROPIC COMPRESSION
		OF NCELLS USING LINKED-LIST SPEED UP

	** CELLS ARE BIDISPERSE, 
		WITH PURELY REPULSIVE INTERACTIONS

	** WILL ALSO PRINT VDOS INFORMATION TO FILE
		AT JAMMING ONSET

	Jack Treado
	08/28/2020, in the time of covid

*/



// preprocessor directives
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

# define NDIM 2
# define NNN 4

// namespaces
using namespace Eigen;
using namespace std;

// GLOBAL CONSTANTS
const double PI 			= 4.0*atan(1);
const int w 				= 10;
const int wnum 				= 25;
const int pnum 				= 14;

// simulation constants
const double phiInit 		= 0.4;
const double phiJMin 		= 0.6;
const double timeStepMag 	= 0.005;
const double sizeRatio 		= 1.4;
const double sizeFraction 	= 0.5;


// FIRE constants for initial minimizations (SP + DP)
const double alpha0      	= 0.2;
const double finc        	= 1.1;
const double fdec        	= 0.5;
const double falpha      	= 0.99;

const int NSKIP 			= 2e4;
const int NMIN        		= 10;
const int NNEGMAX     		= 1000;
const int NDELAY      		= 20;
const int itmax       		= 5e7;


// DP force constants
const double eint 			= 1.0;			// interaction energy
const double del 			= 1.0;			// radius of vertices in units of l0


// FUNCTION PROTOTYPES

// indexing
int gindex(int ci, int vi, vector<int>& szList);
void cindices(int& ci, int& vi, int gi, int NCELLS, vector<int>& szList);

// particle geometry
double area(vector<double>& dpos, int ci, vector<double>& L, vector<int>& nv, vector<int>& szList);
double perimeter(vector<double>& dpos, int ci, vector<double>& L, vector<int>& nv, vector<int>& szList);

// remove rattlers from contact network, return rattler number
int removeRattlers(vector<int>& cij);

// print to file
void printPos(ofstream& posout, vector<double>& vpos, vector<double>& a0, vector<double>& l0, vector<double>& L, vector<int>& cij, vector<int>& nv, vector<int>& szList, double phi, int NCELLS); 


// MAIN
int main(int argc, char const *argv[]){
	// variables for indexing loops
	int i, ci, cj, vi, vj, gi, gj, d;

	// parameters to be read in 
	int NCELLS, smallNV, largeNV, smallN, largeN, NVTOT, NVSMALL, cellDOF, vertDOF, seed;
	double Ptol, Ftol, phi0, dphi, T0, kl, kb, kbb, calA0, smallCalA0, largeCalA0;

	// read in parameters from command line input
	// g++ -O3 -I src sequential/jamming/bidRepulsiveCellJamming.cpp -o jam.o
	// ./jam.o 12 24 1.04 1e-3 1.0 0.01 0 1e-7 1e-12 1 pos.test vdos.test
	string NCELLS_str 		= argv[1];
	string smallNV_str 		= argv[2];
	string calA0_str 		= argv[3];
	string dphi_str 		= argv[4];
	string kl_str 			= argv[5];
	string kb_str 			= argv[6];
	string kbb_str 			= argv[7];
	string Ptol_str 		= argv[8];
	string Ftol_str 		= argv[9];
	string seed_str 		= argv[10];
	string positionFile 	= argv[11];
	string vdosFile  		= argv[12];

	stringstream NCELLSss(NCELLS_str);
	stringstream smallNVss(smallNV_str);
	stringstream calA0ss(calA0_str);
	stringstream dphiss(dphi_str);
	stringstream klss(kl_str);
	stringstream kbss(kb_str);
	stringstream kbbss(kbb_str);
	stringstream Ptolss(Ptol_str);
	stringstream Ftolss(Ftol_str);
	stringstream seedss(seed_str);

	NCELLSss >> NCELLS;
	smallNVss >> smallNV;
	calA0ss >> calA0;
	dphiss >> dphi;
	klss >> kl;
	kbss >> kb;
	kbbss >> kbb;
	Ptolss >> Ptol;
	Ftolss >> Ftol;
	seedss >> seed;

	// open position file (jammed config)
	cout << "opening files" << endl;
	ofstream posout;
	posout.open(positionFile.c_str());
	if (!posout.is_open()){
		cout << "	** ERROR: position file " << positionFile << " could not be opened, ending." << endl;
		return 1;
	}

	// open vdos info file
	ofstream vdosout;
	vdosout.open(vdosFile.c_str());
	if (!vdosout.is_open()){
		cout << "	** ERROR: vdos file " << vdosFile << " could not be opened, ending." << endl;
		return 1;
	}

	// number of vertices on large particles
	largeNV = round(sizeRatio*smallNV);

	// total number of vertices
	smallN 	= round(sizeFraction*NCELLS);
	largeN 	= NCELLS - smallN;
	NVSMALL = smallNV*smallN;
	NVTOT 	= NVSMALL + largeNV*largeN;

	// szList and nv (keep track of global vertex indices)
	cout << "makin vectors" << endl;
	vector<int> szList(NCELLS,0);
	vector<int> nv(NCELLS,0);
	vector<int> hw(NCELLS,0);
	nv.at(0) = smallNV;
	for (ci=1; ci<NCELLS; ci++){
		if (ci < smallN){
			nv.at(ci) = smallNV;
			szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);
		}
		else{
			nv.at(ci) = largeNV;
			szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);
		}

		// determine half-way belt pair depending on number of vertices (odd or even matters)
		if (nv.at(ci) % 2 == 0)
			hw.at(ci) = nv.at(ci)/2;
		else		
			hw.at(ci) = ceil((nv.at(ci) - 1)/2);
	}

	// degree of freedom counts
	cellDOF = NDIM*NCELLS;
	vertDOF = NDIM*NVTOT;

	// save list of adjacent vertices
	vector<int> im1(NVTOT,0);
	vector<int> ip1(NVTOT,0);
	int vim1, vip1;
	for (ci=0; ci<NCELLS; ci++){
		for (vi=0; vi<nv.at(ci); vi++){
			// wrap local indices
			vim1 = (vi - 1 + nv.at(ci)) % nv.at(ci);
			vip1 = (vi + 1) % nv.at(ci);

			// get global wrapped indices
			gi 			= gindex(ci,vi,szList);
			im1.at(gi) 	= gindex(ci,vim1,szList);
			ip1.at(gi) 	= gindex(ci,vip1,szList);
		}
	}

	// shape parameters
	smallCalA0 = calA0*smallNV*tan(PI/smallNV)/PI;
	largeCalA0 = calA0*largeNV*tan(PI/largeNV)/PI;

	cout << "getting times" << endl;

	// fundamental MD time units
	double dtMD, dt0, dt;

	dtMD 	= 1.0;
	dt0 	= timeStepMag*dtMD;
	dt 		= dt0;

	// output opening statement to console
	cout << "=======================================================" << endl << endl;

	cout << "		cellNVE.cpp 									" << endl;
	cout << "		Jack Treado, 2020   							" << endl;
	cout << "		NVE ensemble of deformable particles 			" << endl << endl;

	cout << "		NCELLS 		= " << NCELLS << "					" << endl;
	cout << "		# small 	= " << smallN << "					" << endl;
	cout << "		# large 	= " << largeN << "					" << endl << endl;

	cout << "       small NV 	= " << smallNV << "					" << endl;
	cout << "		large NV 	= " << largeNV << "       			" << endl;
	cout << "		NVTOT 		= " << NVTOT << "					" << endl << endl;

	cout << "		calA0 		= " << calA0 << "					" << endl;
	cout << "		small calA0 = " << smallCalA0 << "				" << endl;
	cout << "		large calA0 = " << largeCalA0 << "				" << endl << endl;

	cout << "		phi0 		= " << phiInit << " 					" << endl;
	cout << "		kl 			= " << kl << "						" << endl;
	cout << "		kb 			= " << kb << "						" << endl;
	cout << "		seed 		= " << seed << "					" << endl << endl;

	cout << "		pos file 	= " << positionFile << "			" << endl;
	cout << "		vdos file 	= " << vdosFile << " 				" << endl << endl;
	
	cout << "=======================================================" << endl << endl;

	// seed random number generator
	srand48(seed);


	/* * * * * * * * * * * * * * * * * * 

				INITIALIZE

					PARTICLES

	 * * * * * * * * * * * * * * * * * */

	// initialization variables
	int nvtmp, hwtmp;
	double a0tmp, lenscale, calA0tmp, areaSum = 0.0;

	// initialize vectors for storing coordinates, shape information
	vector<double> vrad(NVTOT,1.0);
	vector<double> drad(NCELLS,1.0);

	vector<double> vpos(vertDOF,0.0);
	vector<double> dpos(cellDOF,0.0);

	vector<double> a0(NCELLS,1.0);
	vector<double> l0(NCELLS,1.0);
	vector<double> dc(NCELLS,1.0);

	// initialize effective disk radius (for minimization), and l0 parameter
	for (ci=0; ci<NCELLS; ci++){
		// set initial area
		if (ci < smallN){
			lenscale = 1.0;
			a0tmp = 1.0;
			calA0tmp = smallCalA0;
			nvtmp = smallNV;
		}
		else{
			lenscale = sizeRatio;
			a0tmp = sizeRatio*sizeRatio;
			calA0tmp = largeCalA0;
			nvtmp = largeNV;
		}

		// store preferred area
		a0.at(ci) 		= a0tmp;

		// set disk radius
		drad.at(ci) 	= 1.1*sqrt((2.0*a0tmp)/(nvtmp*sin(2.0*PI/nvtmp)));

		// set l0, vector radius
		l0.at(ci) 	= 2.0*lenscale*sqrt(PI*calA0tmp)/nvtmp;
		dc.at(ci) 	= l0.at(ci)/(sin(PI/nvtmp));
		gi 			= szList.at(ci);
		for (vi=0; vi<nvtmp; vi++)
			vrad.at(gi+vi)	= 0.5*l0.at(ci)*del;

		// add to sum of particle areas (including contribution from vertices)
		areaSum 		+= a0tmp + 0.25*PI*pow(l0.at(ci)*del,2.0)*(0.5*nvtmp - 1);
		cout << "drad = " << drad.at(ci) << ", disk area = " << PI*pow(drad.at(ci),2.0) << ", l0 = " << l0.at(ci) << ", a0 = " << a0tmp << ", vrad = " << vrad.at(gi) << endl;
	}

	// determine box lengths from particle sizes and input packing fraction
	vector<double> L(NDIM,1.0);
	for (d=0; d<NDIM; d++)
		L.at(d) = sqrt(areaSum/phiInit);
	phi0 = phiInit;

	// initialize cell centers randomly
	for (ci=0; ci<cellDOF; ci += 2)
		dpos.at(ci) = L[ci % 2]*drand48();
	for (ci=cellDOF-1; ci>0; ci -= 2)
		dpos.at(ci) = L[ci % 2]*drand48();

	// initialize contact network
	int NCTCS = 0.5*NCELLS*(NCELLS-1);
	vector<int> cij(NCTCS,0);



	/* * * * * * * * * * * * * * * * * * 

			INITIALIZE

				BOX-LINKED-LIST

	 * * * * * * * * * * * * * * * * * */


	// Cell-linked-list variables
	double boxLengthScale = 2.0;
	double llscale = l0.at(NCELLS-1);

	// box lengths in each direction
	vector<int> sb(NDIM,0);
	vector<double> lb(NDIM,0.0);
	int NBX = 1;
	for (d=0; d<NDIM; d++){
		// determine number of cells along given dimension by rmax
		sb[d] = round(L[d]/(boxLengthScale*llscale));

		// just in case, if < 3, change to 3 so box neighbor checking will work
		if (sb[d] < 3)
			sb[d] = 3;

		// determine box length by number of cells
		lb[d] = L[d]/sb[d];

		// count total number of cells
		NBX *= sb[d];
	}

	// neighboring boxes for each box (4 neighbors / box)
	int nntmp;
	int scx = sb[0];
	vector< vector<int> > nn;
	nn.resize(NBX);

	// loop over cells, save forward neighbors for each box
	for (i=0; i<NBX; i++){
		// reshape entry
		nn[i].resize(NNN);
		
		// neighbors
		nn[i][0] 			= (i + 1) % NBX; 			// right neighbor (i+1)
		nn[i][1] 			= (i + scx) % NBX;			// top neighbor (j+1)
		nntmp 				= (i + NBX - scx) % NBX;	// bottom neighbor (j-1)
		nn[i][2] 			= (nn[i][1] + 1) % NBX;		// top-right neighbor (i+1, j+1)
		nn[i][3] 			= nntmp + 1;				// bottom-right neighbor (i+1, j-1)

		// right-hand bc (periodic)
		if ((i+1) % scx == 0){
			nn[i][0] = i - scx + 1;
			nn[i][2] = nn[i][1]  - scx + 1;
			nn[i][3] = nntmp - scx + 1;
		}
	} 

	// linked-list variables
	vector<int> head(NBX,0);
	vector<int> last(NBX,0);
	vector<int> list(NVTOT+1,0);












	/* * * * * * * * * * * * * * * * * * 

			INITIAL FIRE

					MINIMZATION

	 * * * * * * * * * * * * * * * * * */


	// ----------------------------
	//
	// S P  M I N I M I Z A T I O N
	// 
	// ----------------------------

	// FIRST SET ALL RADII TO INITIAL VALUES
	int lenReset = 0;
	lenscale = drad.at(NCELLS-1);
	for (ci=1; ci<NCELLS; ci++)
		drad.at(ci) = drad.at(0);

	// initialize disk velocity and force vectors
	vector<double> dv(cellDOF,0.0);
	vector<double> dF(cellDOF,0.0);
	vector<double> dFold(cellDOF,0.0);

	// FIRE VARIABLES
	double P  		= 0;	
	double fnorm 	= 0;
	double vnorm 	= 0;
	double alpha   	= alpha0;

	double dtmax   	= 10*dt0;
	double dtmin   	= 1e-8*dt0;

	int npPos      	= 0;
	int npNeg      	= 0;
	int npPMin      = 0;

	int fireit    	= 0;
	double fcheck  	= 10*Ftol;

	// interaction variables
	double xi, yi, xj, yj, dx, dy, fx, fy, rij, sij, ftmp;

	// loop until force relaxes
	while ((fcheck > Ftol || npPMin < NMIN) && fireit < itmax){
		// VV POSITION UPDATE
		for (i=0; i<cellDOF; i++){
			dpos[i] += dt*dv[i] + 0.5*dt*dt*dF[i];
			dF[i] = 0;
		}


		// FORCE UPDATE
		for (ci=0; ci<NCELLS; ci++){
			xi = dpos[NDIM*ci];
			yi = dpos[NDIM*ci + 1];
			for (cj=ci+1; cj<NCELLS; cj++){
				xj = dpos[NDIM*cj];
				yj = dpos[NDIM*cj + 1];

				// contact distance
				sij = drad[ci] + drad[cj];

				// true distance
				dx = xj - xi;
				dx = dx - L[0]*round(dx/L[0]);
				if (dx < sij){
					dy = yj - yi;
					dy = dy - L[1]*round(dy/L[1]);
					if (dy < sij){
						rij = sqrt(dx*dx + dy*dy);
						if (rij < sij){
							ftmp 				= eint*(1.0 - (rij/sij))/sij;
							fx 					= ftmp*(dx/rij);
							fy 					= ftmp*(dy/rij);

							dF[NDIM*ci] 		-= fx;
							dF[NDIM*ci + 1] 	-= fy;

							dF[NDIM*cj]			+= fx;
							dF[NDIM*cj + 1] 	+= fy;
						}
					}
				}
			}
		}


		// VV VELOCITY UPDATE
		for (i=0; i<cellDOF; i++){
			dv[i] += 0.5*(dF[i] + dFold[i])*dt;
			dFold[i] = dF[i];
		}



		// FIRE UPDATE
		// compute fnorm, vnorm and P
		fnorm = 0.0;
		vnorm = 0.0;
		P = 0.0;
		for (i=0; i<cellDOF; i++){
			fnorm 	+= dF[i]*dF[i];
			vnorm 	+= dv[i]*dv[i];
			P 		+= dv[i]*dF[i];
		}
		fnorm = sqrt(fnorm);
		vnorm = sqrt(vnorm);

		// update fcheck based on fnorm (= force per degree of freedom)
		fcheck = fnorm/sqrt(NDIM*NVTOT);

		// update npPMin
		if (fcheck < Ftol && fireit > NDELAY){
			npPMin++;

			if (lenReset == 0){
				cout << "\t ** Resetting disk radii to avoid segregation..." << endl;
				lenReset = 1;

				// set radii to be larger, avoids segregation
				for (ci=smallN; ci<NCELLS; ci++)
					drad.at(ci) = lenscale;
			}
		}
		else
			npPMin = 0;

		// print to console
		if (fireit % NSKIP == 0){
			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << "		I N I T I A L  S P			" << endl;
			cout << " 	F I R E 						" << endl;
			cout << "		M I N I M I Z A T I O N 	" << endl;
			cout << "===========================================" << endl;
			cout << endl;
			cout << "	** fireit = " << fireit << endl;
			cout << "	** fcheck = " << fcheck << endl;
			cout << "	** vnorm = " << vnorm << endl;
			cout << "	** dt = " << dt << endl;
			cout << "	** P = " << P << endl;
			cout << "	** Pdir = " << P/(fnorm*vnorm) << endl;
			cout << "	** alpha = " << alpha << endl;
		}

		// Step 1. adjust simulation based on net motion of degrees of freedom
		if (P > 0){
			// increase positive counter
			npPos++;

			// reset negative counter
			npNeg = 0;

			// alter simulation if enough positive steps have been taken
			if (npPos > NMIN){
				// change time step
				if (dt*finc < dtmax)
					dt *= finc;

				// decrease alpha
				alpha *= falpha;
			}
		}
		else{
			// reset positive counter
			npPos = 0;

			// increase negative counter
			npNeg++;

			// check if simulation is stuck
			if (npNeg > NNEGMAX){
				cout << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
				return 1;
			}

			// take half step backwards, reset velocities
			for (i=0; i<cellDOF; i++){
				// take half step backwards
				dpos[i] -= 0.5*dt*dv[i] + 0.25*dt*dt*dF[i];

				// reset velocities
				dv[i] = 0.0;
			}

			// decrease time step if past initial delay
			if (fireit > NDELAY){
				// decrease time step 
				if (dt*fdec > dtmin)
					dt *= fdec;

				// reset alpha
				alpha = alpha0;
			}
		}


		// update velocities (s.d. vs inertial dynamics) only if forces are acting
		if (fnorm > 0){
			for (i=0; i<cellDOF; i++)
				dv[i] = (1 - alpha)*dv[i] + alpha*(vnorm/fnorm)*dF[i];
		}

		// update iterator
		fireit++;
	}
	// check if FIRE converged
	if (fireit == itmax){
		cout << "	** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
		return 1;
	}
	else{
		cout << endl << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===========================================" << endl;
		cout << " 	F I R E 						" << endl;
		cout << "		M I N I M I Z A T I O N 	" << endl;
		cout << "	C O N V E R G E D! 				" << endl << endl;

		cout << "	(for initial disk minimization) " << endl;
		cout << "===========================================" << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << endl;
		cout << "	** fireit = " << fireit << endl;
		cout << "	** fcheck = " << fcheck << endl;
		cout << "	** vnorm = " << vnorm << endl;
		cout << "	** dt = " << dt << endl;
		cout << "	** P = " << P << endl;
		cout << "	** Pdir = " << P/(fnorm*vnorm) << endl;
		cout << "	** alpha = " << alpha << endl;
	}

	// initialize vertex positions based on cell centers
	for (ci=0; ci<NCELLS; ci++){
		for (vi=0; vi<nv.at(ci); vi++){
			// get global vertex index
			gi = gindex(ci,vi,szList);

			// get distance from cell center to vertex if reg poly
			lenscale = sqrt((2.0*a0.at(ci))/(nv.at(ci)*sin((2.0*PI)/nv.at(ci))));

			// set vertex positions
			vpos.at(NDIM*gi) 		= lenscale*cos((2.0*PI*vi)/nv.at(ci)) + dpos.at(NDIM*ci) + 1e-2*l0[ci]*drand48();
			vpos.at(NDIM*gi + 1)	= lenscale*sin((2.0*PI*vi)/nv.at(ci)) + dpos.at(NDIM*ci + 1) + 1e-2*l0[ci]*drand48();
		}
	}























	/* * * * * * * * * * * * * * * * * * 

			COMPRESS TO

				JAMMING ONSET

	 * * * * * * * * * * * * * * * * * */

	// initialize disk velocity and force vectors
	vector<double> vvel(vertDOF,0.0);
	vector<double> vF(vertDOF,0.0);

	// jamming check variables
	bool undercompressed, overcompressed, jammed;
	int k, kmax, xind, yind;
	double rH, r0, rL, drgrow, drshrink, scaleFactor;
	double pcheck;

	// max number of jamming iteractions
	k = 0;
	kmax = 1e4;

	// initial boolean variables
	undercompressed = 0;
	overcompressed = 0;
	jammed = 0;

	// jamming bounds
	rH = -1;
	rL = -1;

	// compute scale factors
	drgrow = sqrt((phi0 + dphi)/phi0);
	drshrink = sqrt((phi0 - 0.5*dphi)/phi0);

	// saved state
	vector<double> vposSave(vertDOF,0.0);
	vector<double> vradSave(NVTOT,0.0);
	vector<double> a0Save(NCELLS,0.0);
	vector<double> l0Save(NCELLS,0.0);
	vector<double> dcSave(NCELLS,0.0);


	// save initial state
	vposSave = vpos;
	vradSave = vrad;
	a0Save = a0;
	l0Save = l0;

	// linked list variables
	int boxid, bi, bj, pi, pj, sbtmp;
	int d0, dend;

	// temporary tolerance (to speed initial compression)
	double Ftoltmp;
	if (Ftol < 0.01*Ptol)
		Ftoltmp = 0.01*Ptol;

	// total potential energy
	double U = 0.0;

	// length unit variable
	double rho0 = 0.0;

	// shape force variables
	double fa, fl, fb, fbb, l0tmp, dctmp, atmp, li, lim1, kappai, cx, cy;
	double da, dli, dlim1;
	double lim2x, lim2y, lim1x, lim1y, lix, liy, lip1x, lip1y;
	double rim2x, rim2y, rim1x, rim1y, rix, riy, rip1x, rip1y, rip2x, rip2y;
	double lbbx, lbby, lbb, meanL;
	double ua, ul, ub;

	// contact variables
	int Nvv, Ncc, nr;

	// compress to jamming, relax U and F using FIRE
	while (!jammed && k < kmax){
		// update iterator
		k++;

		// check to see if cell linked-list needs to be redrawn
		if ((l0[NCELLS-1]/llscale) > 0.9*boxLengthScale){
			// print to console
			cout << "\t ** Resetting linked-list: old llscale = " << llscale << ", new llscale = " << l0[NCELLS-1] << endl;
			cout << "\t ** Resetting linked-list: old NBX = " << NBX; 


			// reset linked list for larger particles
			llscale = l0[NCELLS-1];

			// reset box length info
			fill(sb.begin(), sb.end(), 0);
			fill(lb.begin(), lb.end(), 0.0);
			for (i=0; i < NBX; i++)
				nn[i].clear();


			NBX = 1;
			for (d=0; d<NDIM; d++){
				// determine number of cells along given dimension by rmax
				sb[d] = round(L[d]/(boxLengthScale*llscale));

				// just in case, if < 3, change to 3 so box neighbor checking will work
				if (sb[d] < 3)
					sb[d] = 3;

				// determine box length by number of cells
				lb[d] = L[d]/sb[d];

				// count total number of cells
				NBX *= sb[d];
			}
			cout << ", new NBX = " << NBX << endl;

			// neighboring boxes for each box (4 neighbors / box)
			scx = sb[0];
			nn.clear();
			nn.resize(NBX);

			// loop over cells, save forward neighbors for each box
			for (i=0; i<NBX; i++){
				// reshape entry
				nn[i].resize(NNN);
				
				// neighbors
				nn[i][0] 			= (i + 1) % NBX; 			// right neighbor (i+1)
				nn[i][1] 			= (i + scx) % NBX;			// top neighbor (j+1)
				nntmp 				= (i + NBX - scx) % NBX;	// bottom neighbor (j-1)
				nn[i][2] 			= (nn[i][1] + 1) % NBX;		// top-right neighbor (i+1, j+1)
				nn[i][3] 			= nntmp + 1;				// bottom-right neighbor (i+1, j-1)

				// right-hand bc (periodic)
				if ((i+1) % scx == 0){
					nn[i][0] = i - scx + 1;
					nn[i][2] = nn[i][1]  - scx + 1;
					nn[i][3] = nntmp - scx + 1;
				}
			}

			// reset head/last/list
			head.clear();
			head.resize(NBX);

			last.clear();
			last.resize(NBX);

			fill(list.begin(), list.end(), 0);
		}

		// // update tolerance
		if (phi0 > 0.6)
			Ftoltmp = Ftol;

		// RESET FIRE VARIABLES
		P  			= 0;	
		fnorm 		= 0;
		vnorm 		= 0;
		alpha   	= alpha0;

		dtmax   	= 10.0*dt0;
		dtmin   	= 1e-2*dt0;
		dt 			= dt0;

		npPos      	= 0;
		npNeg      	= 0;
		npPMin      = 0;

		fireit    	= 0;
		fcheck  	= 10*Ftoltmp;

		// reset forces
		fill(vF.begin(), vF.end(), 0.0);
		fill(vvel.begin(), vvel.end(), 0.0);

		// set length constant (variable due to particle growth)
		rho0 = sqrt(a0.at(0));

		// RELAX FORCES USING FIRE
		while (fcheck > Ftoltmp && fireit < itmax){
			// compute P
			P = 0.0;
			for (i=0; i<vertDOF; i++)
				P += vvel[i]*vF[i];

			// print to console
			if (fireit % NSKIP == 0){
				cout << endl << endl;
				cout << "===========================================" << endl;
				cout << " 	F I R E 						" << endl;
				cout << "		M I N I M I Z A T I O N 	" << endl;
				cout << "===========================================" << endl;
				cout << endl;
				cout << "	** fireit 	= " << fireit << endl;
				cout << "	** fcheck 	= " << fcheck << endl;
				cout << "	** pcheck 	= " << pcheck << endl;
				cout << "	** U 		= " << U << endl;

				cout << "	** dt 		= " << dt << endl;
				cout << "	** P 		= " << P << endl;
				cout << "	** alpha 	= " << alpha << endl;
				cout << "	** npPos 	= " << npPos << endl;
				cout << "	** npNeg 	= " << npNeg << endl;
			}

			// Step 1. adjust simulation based on net motion of degrees of freedom
			if (P > 0){
				// increase positive counter
				npPos++;

				// reset negative counter
				npNeg = 0;

				// alter simulation if enough positive steps have been taken
				if (npPos > NDELAY){
					// change time step
					if (dt*finc < dtmax)
						dt *= finc;

					// decrease alpha
					alpha *= falpha;
				}
			}
			else{
				// reset positive counter
				npPos = 0;

				// increase negative counter
				npNeg++;

				// check if simulation is stuck
				if (npNeg > NNEGMAX){
					cout << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
					return 1;
				}

				// take half step backwards, reset velocities
				for (i=0; i<vertDOF; i++){
					// take half step backwards
					vpos[i] -= 0.5*dt*vvel[i];

					// reset vertex velocities
					vvel[i] = 0.0;
				}

				// decrease time step if past initial delay
				if (fireit > NDELAY){
					// decrease time step 
					if (dt*fdec > dtmin)
						dt *= fdec;

					// reset alpha
					alpha = alpha0;
				}
			}

			// VV VELOCITY UPDATE #1
			for (i=0; i<vertDOF; i++)
				vvel[i] += 0.5*dt*vF[i];

			// compute fnorm, vnorm and P
			fnorm = 0.0;
			vnorm = 0.0;
			for (i=0; i<vertDOF; i++){
				fnorm 	+= vF[i]*vF[i];
				vnorm 	+= vvel[i]*vvel[i];
			}
			fnorm = sqrt(fnorm);
			vnorm = sqrt(vnorm);

			// update velocities (s.d. vs inertial dynamics) only if forces are acting
			if (fnorm > 0){
				for (i=0; i<vertDOF; i++)
					vvel[i] = (1 - alpha)*vvel[i] + alpha*(vF[i]/fnorm)*vnorm;
			}


			// VV POSITION UPDATE
			for (i=0; i<vertDOF; i++){
				// update position
				vpos[i] += dt*vvel[i];

				// recenter in box
				if (vpos[i] > L[i % NDIM])
					vpos[i] -= L[i % NDIM];
				else if (vpos[i] < 0)
					vpos[i] += L[i % NDIM];

				// reset forces
				vF[i] = 0;
			}

			// reset linked list info
			fill(list.begin(), list.end(), 0);
			fill(head.begin(), head.end(), 0);
			fill(last.begin(), last.end(), 0);

			// sort vertices into linked list
			for (gi=0; gi<NVTOT; gi++){
				// 1. get cell id of current particle position
				boxid = 0;
				sbtmp = 1;
				for (d=0; d<NDIM; d++){
					// add d index to 1d list
					boxid += floor(vpos[NDIM*gi + d]/lb[d])*sbtmp;

					// increment dimensional factor
					sbtmp *= sb[d];
				}

				// 2. add to head list or link within list
				// NOTE: particle ids are labelled starting from 1, setting to 0 means end of linked list
				if (head[boxid] == 0){
					head[boxid] = gi + 1;
					last[boxid] = gi + 1;
				}
				else{
					list[last[boxid]] = gi + 1;
					last[boxid] = gi + 1;
				}
			}

			// reset contact network
			fill(cij.begin(), cij.end(), 0);

			// FORCE UPDATE

			// interaction forces (USE BOX LINKED LIST)
			U = 0.0;
			pcheck = 0.0;
			for (bi=0; bi<NBX; bi++){

				// get start of list of particles
				pi = head[bi];

				// loop over linked list
				while (pi > 0){
					// real particle index
					gi = pi - 1;

					// next particle in list
					pj = list[pi];

					// loop down neighbors of pi in same cell
					while (pj > 0){
						// real index of pj
						gj = pj - 1;

						if (gj == ip1[gi] || gj == im1[gi]){
							pj = list[pj];
							continue;
						}

						// contact distance
						sij = vrad[gi] + vrad[gj];

						// particle distance
						dx = vpos[NDIM*gj] - vpos[NDIM*gi];
						dx -= L[0]*round(dx/L[0]);
						if (dx < sij){
							dy = vpos[NDIM*gj + 1] - vpos[NDIM*gi + 1];
							dy -= L[1]*round(dy/L[1]);
							if (dy < sij){
								rij = sqrt(dx*dx + dy*dy);
								if (rij < sij){
									// force scale
									ftmp 				= eint*(1 - (rij/sij))*(rho0/sij);
									fx 					= ftmp*(dx/rij);
									fy 					= ftmp*(dy/rij);

									// add to forces
									vF[NDIM*gi] 		-= fx;
									vF[NDIM*gi + 1] 	-= fy;

									vF[NDIM*gj] 		+= fx;
									vF[NDIM*gj + 1] 	+= fy;

									// increae potential energy
									U += 0.5*eint*pow((1 - (rij/sij)),2.0);

									// add to virial expression for pressure
									pcheck += dx*fx + dy*fy;

									// add to contacts
									cindices(ci, vi, gi, NCELLS, szList);
									cindices(cj, vj, gj, NCELLS, szList);

									if (ci > cj)
										cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
									else if (ci < cj)
										cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++; 
								}
							}
						}

						// update pj
						pj = list[pj];
					}

					// test overlaps with forward neighboring cells
					for (bj=0; bj<NNN; bj++){
						// get first particle in neighboring cell
						pj = head[nn[bi][bj]];

						// loop down neighbors of pi in same cell
						while (pj > 0){
							// real index of pj
							gj = pj - 1;

							if (gj == ip1[gi] || gj == im1[gi]){
								pj = list[pj];
								continue;
							}
							// contact distance
							sij = vrad[gi] + vrad[gj];

							// particle distance
							dx = vpos[NDIM*gj] - vpos[NDIM*gi];
							dx -= L[0]*round(dx/L[0]);
							if (dx < sij){
								dy = vpos[NDIM*gj + 1] - vpos[NDIM*gi + 1];
								dy -= L[1]*round(dy/L[1]);
								if (dy < sij){
									rij = sqrt(dx*dx + dy*dy);
									if (rij < sij){
										// force scale
										ftmp 				= eint*(1 - (rij/sij))*(rho0/sij);
										fx 					= ftmp*(dx/rij);
										fy 					= ftmp*(dy/rij);

										// add to forces
										vF[NDIM*gi] 		-= fx;
										vF[NDIM*gi + 1] 	-= fy;

										vF[NDIM*gj] 		+= fx;
										vF[NDIM*gj + 1] 	+= fy;

										// increae potential energy
										U += 0.5*eint*pow((1 - (rij/sij)),2.0);

										// add to virial expression for pressure
										pcheck += dx*fx + dy*fy;

										// add to contacts
										cindices(ci, vi, gi, NCELLS, szList);
										cindices(cj, vj, gj, NCELLS, szList);

										if (ci > cj)
											cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
										else if (ci < cj)
											cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++; 
									}
								}
							}

							// update pj
							pj = list[pj];
						}
					}

					// update pi index to be next
					pi = list[pi];
				}
			}

			// normalize pressure by box area (make dimensionless with extra factor of rho)
			pcheck *= (rho0/(2.0*L[0]*L[1]));


			// shape forces (loop over global vertex labels)
			ci = 0;
			for (gi=0; gi<NVTOT; gi++){

				// -- Area force (and get cell index ci)
				if (ci < NCELLS){
					if (gi == szList[ci]){
						// compute shape parameter
						nvtmp = nv[ci];
						hwtmp = hw[ci];
						a0tmp = a0[ci];
						l0tmp = l0[ci];
						dctmp = dc[ci];

						// compute area deviation
						atmp = area(vpos,ci,L,nv,szList);
						da = (atmp/a0tmp) - 1.0;

						// shape force parameters (kl and kl are unitless energy ratios)
						fa = da*(rho0/a0tmp);		// derivation from the fact that rho0^2 does not necessarily cancel a0tmp
						fl = kl*(rho0/l0tmp);
						fb = kb*(rho0/(l0tmp*l0tmp));
						fbb = kbb*(rho0/(hwtmp*dctmp));
						
						// compute cell center of mass
						xi = vpos[NDIM*gi];
						yi = vpos[NDIM*gi + 1];
						cx = xi; 
						cy = yi;
						for (vi=1; vi<nvtmp; vi++){
							dx = vpos.at(NDIM*(gi+vi)) - xi;
							dx -= L[0]*round(dx/L[0]);

							dy = vpos.at(NDIM*(gi+vi) + 1) - yi;
							dy -= L[1]*round(dy/L[1]);

							xi += dx;
							yi += dy;

							cx += xi;
							cy += yi;
						}
						cx /= nvtmp;
						cy /= nvtmp;

						// compute mean belt length
						if (kbb > 0){
							meanL = 0.0;
							for (vi=0; vi<hwtmp; vi++){
								// get current belt vector
								lbbx = vpos[NDIM*(gi + vi + hwtmp)] - vpos[NDIM*(gi + vi)];
								lbbx -= L[0]*round(lbbx/L[0]);

								lbby = vpos[NDIM*(gi + vi + hwtmp) + 1] - vpos[NDIM*(gi + vi) + 1];
								lbby -= L[1]*round(lbby/L[1]);

								meanL += sqrt(lbbx*lbbx + lbby*lbby);
							}
							meanL /= hwtmp;
						}

						// get coordinates relative to center of mass
						rix = vpos[NDIM*gi] - cx;
						riy = vpos[NDIM*gi + 1] - cy;

						// get (prior) adjacent vertices
						rim1x = vpos[NDIM*im1[gi]] - cx;
						rim1x -= L[0]*round(rim1x/L[0]);

						rim1y = vpos[NDIM*im1[gi] + 1] - cy;
						rim1y -= L[1]*round(rim1y/L[1]);

						rim2x = vpos[NDIM*im1[im1[gi]]] - cx;
						rim2x -= L[0]*round(rim2x/L[0]);

						rim2y = vpos[NDIM*im1[im1[gi]] + 1] - cy;
						rim2y -= L[1]*round(rim2y/L[1]);

						// increment cell index
						ci++;
					}
				}


				// get next adjacent vertices
				rip1x = vpos.at(NDIM*ip1[gi]) - cx;
				rip1x -= L[0]*round(rip1x/L[0]);

				rip1y = vpos.at(NDIM*ip1[gi] + 1) - cy;
				rip1y -= L[1]*round(rip1y/L[1]);


				// -- Area force
				vF[NDIM*gi] 		+= 0.5*fa*(rim1y - rip1y);
				vF[NDIM*gi + 1] 	+= 0.5*fa*(rip1x - rim1x);


				// -- Perimeter force

				// segment vector elements
				lim1x 	= rix - rim1x;
				lim1y 	= riy - rim1y;

				lix 	= rip1x - rix;
				liy 	= rip1y - riy;

				// segment lengths
				lim1 	= sqrt(lim1x*lim1x + lim1y*lim1y);
				li 		= sqrt(lix*lix + liy*liy);

				// segment deviations
				dlim1  	= (lim1/l0tmp) - 1.0;
				dli 	= (li/l0tmp) - 1.0;

				// add to forces
				vF[NDIM*gi] 		+= fl*(dli*(lix/li) - dlim1*(lim1x/lim1));
				vF[NDIM*gi + 1] 	+= fl*(dli*(liy/li) - dlim1*(lim1y/lim1));


				// -- Bending force
				if (kb > 0){
					// segment vectors for ip2
					rip2x = vpos[NDIM*ip1[ip1[gi]]] - cx;
					rip2x -= L[0]*round(rip2x/L[0]);

					rip2y = vpos[NDIM*ip1[ip1[gi]] + 1] - cy;
					rip2y -= L[1]*round(rip2y/L[1]);

					lip1x = rip2x - rip1x;
					lip1y = rip2y - rip1y;

					lim2x = rim1x - rim2x;
					lim2y = rim1y - rim2y;

					// add to force
					vF[NDIM*gi] 		+= fb*(3.0*(lix - lim1x) + lim2x - lip1x);
					vF[NDIM*gi + 1] 	+= fb*(3.0*(liy - lim1y) + lim2y - lip1y);
				}

				// -- Belt force

				// segment vectors across particle (bc of int algebra, (nv-1)/2 either nv/2 (even) or (nv-1)/2 (odd))
				if (kbb > 0 && gi - szList[ci-1] < hwtmp){
					// get current belt vector
					lbbx = vpos[NDIM*(gi + hwtmp)] - vpos[NDIM*gi];
					lbbx -= L[0]*round(lbbx/L[0]);

					lbby = vpos[NDIM*(gi + hwtmp) + 1] - vpos[NDIM*gi + 1];
					lbby -= L[1]*round(lbby/L[1]);

					lbb = sqrt(lbbx*lbbx + lbby*lbby);

					// add to force
					vF[NDIM*gi] 					+= fbb*((meanL/dctmp) - 1.0)*(lbbx/lbb);
					vF[NDIM*gi + 1] 				+= fbb*((meanL/dctmp) - 1.0)*(lbby/lbb);

					vF[NDIM*(gi + hwtmp)] 	-= fbb*((meanL/dctmp) - 1.0)*(lbbx/lbb);
					vF[NDIM*(gi + hwtmp) + 1] -= fbb*((meanL/dctmp) - 1.0)*(lbby/lbb);
				}

				// update old coordinates
				rim2x = rim1x;
				rim1x = rix;
				rix = rip1x;

				rim2y = rim1y;
				rim1y = riy;
				riy = rip1y;
			}


			// VV VELOCITY UPDATE
			for (i=0; i<vertDOF; i++)
				vvel[i] += 0.5*vF[i]*dt;

			// update fcheck based on fnorm (= force per degree of freedom)
			fcheck = 0.0;
			for (i=0; i<vertDOF; i++)
				fcheck += vF[i]*vF[i];
			fcheck = sqrt(fcheck/vertDOF);

			// update iterator
			fireit++;
		}
		// check if FIRE converged
		if (fireit == itmax){
			cout << "	** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
			return 1;
		}
		else{
			cout << endl;
			cout << "===========================================" << endl;
			cout << " 	F I R E 						" << endl;
			cout << "		M I N I M I Z A T I O N 	" << endl;
			cout << "	C O N V E R G E D! 				" << endl;
			cout << "	** at k = " << k << " 			" << endl;
			cout << "===========================================" << endl;
			cout << endl;
			cout << "	** fireit 	= " << fireit << endl;
			cout << "	** fcheck 	= " << fcheck << endl;
			cout << "	** pcheck 	= " << pcheck << endl;
			cout << "	** U 		= " << U << endl;

			cout << "	** vnorm 	= " << vnorm << endl;
			cout << "	** dt 		= " << dt << endl;
			cout << "	** P 		= " << P << endl;
			cout << "	** Pdir 	= " << P/(fnorm*vnorm) << endl;
			cout << "	** alpha 	= " << alpha << endl;
			cout << endl << endl;
		}

		// remove rattlers
		nr = removeRattlers(cij);

		// update number of contacts
		Nvv = 0;
		Ncc = 0;
		for (i=0; i<NCTCS; i++){
			Nvv += cij[i];
			if (cij[i] > 0)
				Ncc++;
		}

		// boolean check for jamming
		undercompressed = ((pcheck < 2.0*Ptol && rH < 0) || (pcheck < Ptol && rH > 0));
		overcompressed = (pcheck > 2.0*Ptol && phi0 > phiJMin);
		jammed = (pcheck < 2.0*Ptol && pcheck > Ptol && rH > 0 && rL > 0);

		// output to console
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===============================================" << endl << endl;
		cout << " 	Q U A S I S T A T I C  				" << endl;
		cout << " 	  	I S O T R O P I C 				" << endl;
		cout << "			C O M P R E S S I O N 		" << endl << endl;
		cout << "===============================================" << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* dphi 			= " << dphi << endl;
		cout << "	* phi 			= " << phi0 << endl;
		cout << "	* r0 			= " << r0 << endl;
		cout << "	* rH 			= " << rH << endl;
		cout << "	* rL 			= " << rL << endl;
		cout << "	* fcheck 		= " << fcheck << endl;
		cout << "	* pcheck 		= " << pcheck << endl;
		cout << "	* U 		 	= " << U << endl;
		cout << "	* Nvv 			= " << Nvv << endl;
		cout << "	* Ncc 			= " << Ncc << endl;
		cout << "	* # of rattlers = " << nr << endl;
		cout << "	* undercompressed = " << undercompressed << endl;
		cout << "	* overcompressed = " << overcompressed << endl;
		cout << "	* jammed = " << jammed << endl << endl;
		cout << endl;

		// cout << " 	* Printing positions to file" << endl;
		// printPos(posout, vpos, a0, l0, L, cij, nv, szList, phi0, NCELLS);

		// update particle sizes based on target check
		if (rH < 0){
			// if still undercompressed, then grow until overcompressed found
			if (undercompressed){
				r0 = rho0;
				scaleFactor = sqrt((phi0 + dphi)/phi0);
			}
			// if first overcompressed, decompress by dphi/2 until unjamming
			else if (overcompressed){
				// current = upper bound length scale r
	            rH = rho0;

	            // save first overcompressed state
				r0 = rH;
				vposSave = vpos;
				vradSave = vrad;
				a0Save = a0;
				l0Save = l0;
				dcSave = dc;

	            // compute new scale factor
	            scaleFactor = drshrink;

	            // print to console
				cout << "	-- -- overcompressed for the first time, scaleFactor = " << scaleFactor << endl;
			}
		}
		else{
			if (rL < 0){
				// if first undercompressed, save last overcompressed state, beging root search
				if (undercompressed){
					// current = new lower bound length scale r
					rL = rho0;

					// load state
					vpos = vposSave;
					vrad = vradSave;
					a0 = a0Save;
					l0 = l0Save;
					dc = dcSave;

					// compute new scale factor by root search
		            scaleFactor = 0.5*(rH + rL)/r0;

					// print to console
					cout << "	-- -- undercompressed for the first time, scaleFactor = " << scaleFactor << endl;
					cout << "	-- -- BEGINNING ROOT SEARCH IN ENTHALPY MIN PROTOCOL..." << endl;
				}
				// if still overcompressed, decrement again
				else if (overcompressed){
					// current = upper bound length scale r
		            rH = rho0;

		            // save overcompressed state
					r0 = rH;
					vposSave = vpos;
					vradSave = vrad;
					a0Save = a0;
					l0Save = l0;
					dcSave = dc;

		            // keep shrinking at same rate until unjamming
		            scaleFactor = drshrink;

		            // print to console
					cout << "	-- -- overcompressed, still no unjamming, scaleFactor = " << scaleFactor << endl;
				}
			}
			else{
				// if found undercompressed state, go to state between undercompressed and last overcompressed states (from saved state)
				if (undercompressed){
					// current = new lower bound length scale r
					rL = rho0;

					// load state
					vpos = vposSave;
					vrad = vradSave;
					a0 = a0Save;
					l0 = l0Save;
					dc = dcSave;

					// compute new scale factor
		            scaleFactor = 0.5*(rH + rL)/r0;

					// print to console
					cout << "	-- -- undercompressed, scaleFactor = " << scaleFactor << endl;

				}
				else if (overcompressed){
					// current = new upper bound length scale r
		            rH = rho0;

					// load state
					vpos = vposSave;
					vrad = vradSave;
					a0 = a0Save;
					l0 = l0Save;
					dc = dcSave;

					// compute new scale factor
		            scaleFactor = 0.5*(rH + rL)/r0;

					// print to console
					cout << "	-- -- overcompressed, scaleFactor = " << scaleFactor << endl;
				}
				else if (jammed){
					cout << "	** At k = " << k << ", target pressure found!" << endl;
					cout << "	** fcheck = " << fcheck << endl;
					cout << "	** pcheck = " << pcheck << endl;
					cout << "	** U = " << U << endl;
					cout << "	** Nvv = " << Nvv << endl;
					cout << "	** Ncc = " << Ncc << endl;
					cout << " WRITING ENTHALPY-MINIMIZED CONFIG TO .jam FILE" << endl;
					cout << " ENDING COMPRESSION SIMULATION" << endl;
					printPos(posout, vpos, a0, l0, L, cij, nv, szList, phi0, NCELLS);
					break;
				}
			}
		}

		// grow or shrink particles by scale factor
		phi0 = 0.0;
		for (ci=0; ci<NCELLS; ci++){
			// scale preferred lengths
			l0[ci] *= scaleFactor;
			dc[ci] *= scaleFactor;
			a0[ci] *= scaleFactor*scaleFactor;

			// first global index for ci
			gi = szList.at(ci);

			// compute cell center of mass
			xi = vpos[NDIM*gi];
			yi = vpos[NDIM*gi + 1];
			cx = xi; 
			cy = yi;
			for (vi=1; vi<nv.at(ci); vi++){
				dx = vpos.at(NDIM*(gi+vi)) - xi;
				dx -= L[0]*round(dx/L[0]);

				dy = vpos.at(NDIM*(gi+vi) + 1) - yi;
				dy -= L[1]*round(dy/L[1]);

				xi += dx;
				yi += dy;

				cx += xi;
				cy += yi;
			}
			cx /= nv.at(ci);
			cy /= nv.at(ci);

			for (vi=0; vi<nv.at(ci); vi++){
				// x and y inds
				xind = NDIM*(gi+vi);
				yind = xind + 1;

				// closest relative position
				dx = vpos[xind] - cx;
				dx -= L[0]*round(dx/L[0]);

				dy = vpos[yind] - cy;
				dy -= L[1]*round(dy/L[1]);

				// update vertex positions
				vpos[xind] 		+= (scaleFactor - 1.0)*dx;
				vpos[yind] 		+= (scaleFactor - 1.0)*dy;

				// scale vertex radii
				vrad[gi+vi] *= scaleFactor;
			}

			// update packing fraction
			phi0 += a0[ci] + 0.25*PI*pow(l0[ci]*del,2.0)*(0.5*nv[ci] - 1);
		}
		phi0 /= L[0]*L[1];
	}
	if (k == kmax){
		cout << "	** ERROR: IN 2d cell jamming finding, k reached kmax without finding jamming. Ending." << endl;
		return 1;
	}


	/* * * * * * * * * * * * * * * * * * 

			   COMPUTE VDOS

	 * * * * * * * * * * * * * * * * * */


	// LOCAL VARIABLES

	// integers
	int l, lxm1, lym1, lxp1, lyp1, lxbb, lybb, kxbb, kybb;
	int vim2, vip2, vjm1, vjp1;
	int kxm2, kym2, kxm1, kym1, kx, kxp1, kxp2, ky, kyp1, kyp2, lx, ly;
	int mxi, myi, mxj, myj;
	int inContact;

	// doubles
	double delim1, deli, delA;
	double ljm1x, ljm1y, ljx, ljy;
	double ulim1x, ulim1y, ulix, uliy;
	double da_dxi, da_dyi, da_dxj, da_dyj;
	double kapim1, kapi, kapip1;
	double dkapi_dxi, dkapi_dyi, dkapip1_dxi, dkapip1_dyi, dkapim1_dxi, dkapim1_dyi;
	double dkapi_dxip1, dkapi_dyip1, dkapip1_dxip1, dkapip1_dyip1, dkapip1_dxip2, dkapip1_dyip2;
	double dL_dxi, dL_dyi, dL_dxj, dL_dyj; // should really be called D instead of L to match PRM
	double delBB, Lbb_kx, Lbb_ky, Lbb_k, uLbb_kx, uLbb_ky, Lbb_lx, Lbb_ly, Lbb_l, uLbb_lx, uLbb_ly;
	double Kl1, Kl2, Kb1, Kb2;
	double kij, h;
	double uxij, uyij;


	// initialize matrices
	// NOTE: surface tension VDOS not yet supported
	Eigen::MatrixXd Ha(vertDOF,vertDOF);		// stiffness matrix for area term
	Eigen::MatrixXd Sa(vertDOF,vertDOF);		// stress matrix for area term
	Eigen::MatrixXd Hl(vertDOF,vertDOF);		// stiffness matrix for perimeter term
	Eigen::MatrixXd Sl(vertDOF,vertDOF);		// stress matrix for perimeter term
	Eigen::MatrixXd Hb(vertDOF,vertDOF);		// stiffness matrix for bending energy
	Eigen::MatrixXd Sb(vertDOF,vertDOF);		// stress matrix for bending term
	Eigen::MatrixXd Hbb(vertDOF,vertDOF);		// stiffness matrix for belt term
	Eigen::MatrixXd Sbb(vertDOF,vertDOF);		// stiffness matrix for belt term
	Eigen::MatrixXd Hvv(vertDOF,vertDOF);		// stiffness matrix for interaction terms
	Eigen::MatrixXd Svv(vertDOF,vertDOF);		// stress matrix for interaction terms
	Eigen::MatrixXd H(vertDOF,vertDOF);			// stiffness matrix
	Eigen::MatrixXd S(vertDOF,vertDOF);			// stress matrix
	Eigen::MatrixXd M(vertDOF,vertDOF);			// full dynamical matrix	

	// initialize all matrices to be 0 initially
	for (k=0; k<vertDOF; k++){
		for (l=0; l<vertDOF; l++){
			Ha(k,l) = 0.0;
			Sa(k,l) = 0.0;
			Hl(k,l) = 0.0;
			Sl(k,l) = 0.0;
			Hb(k,l) = 0.0;
			Sb(k,l) = 0.0;
			Hbb(k,l) = 0.0;
			Sbb(k,l) = 0.0;
			Hvv(k,l) = 0.0;
			Svv(k,l) = 0.0;
			S(k,l) = 0.0;
			H(k,l) = 0.0;
			M(k,l) = 0.0;
		}
	}

	// Loop over cells, compute shape forces for each individual cell and contributions from
	// vertex-vertex interactions
	cout << "	** COMPUTING VDOS ... " << endl;
	rho0 = sqrt(a0.at(0));
	for (ci=0; ci<NCELLS; ci++){

		// print statement
		cout << "		-- Computing dynamical matrix elements for cell ci = " << ci << endl;
		

		// ------------------------------------------
		// 
		// 				SHAPE
		// 					CONTRIBUTIONS
		//
		// ------------------------------------------

		// number of vertices of cell ci
		nvtmp = nv[ci];
		hwtmp = hw[ci];

		// geometric factors
		l0tmp = l0[ci];
		a0tmp = a0[ci];
		dctmp = dc[ci];

		// area deviations
		delA = (area(vpos,ci,L,nv,szList)/a0tmp) - 1.0;

		// dimensionless stiffness constants
		Kl1 = kl*(rho0*rho0)/l0tmp;				// units = L
		Kl2 = Kl1/l0tmp;						// units = 1

		Kb1 = kb*(rho0*rho0);					// units = L^2
		Kb2 = Kb1/(l0tmp*l0tmp);				// units = 1

		// compute mean belt length
		delBB = 0.0;
		meanL = 0.0;
		if (kbb > 0){
			gi = szList[ci];
			for (vi=0; vi<hwtmp; vi++){
				lbbx = vpos[NDIM*(gi + vi + hwtmp)] - vpos[NDIM*(gi + vi)];
				lbbx -= L[0]*round(lbbx/L[0]);

				lbby = vpos[NDIM*(gi + vi + hwtmp) + 1] - vpos[NDIM*(gi + vi) + 1];
				lbby -= L[1]*round(lbby/L[1]);

				meanL += sqrt(lbbx*lbbx + lbby*lbby);
			}
			meanL /= hwtmp;
			delBB = (meanL/dctmp) - 1.0;
		}

		// loop over vertices, compute each DM element
		for (vi=0; vi<nvtmp; vi++){

			// wrap vertices
			vim2 		= (vi - 2 + nvtmp) % nvtmp;
			vim1 		= (vi - 1 + nvtmp) % nvtmp;
			vip1 		= (vi + 1) % nvtmp;
			vip2 		= (vi + 2) % nvtmp;			

			// vertex elements
			kxm2 		= NDIM*(szList.at(ci) + vim2);
			kym2 		= NDIM*(szList.at(ci) + vim2) + 1;

			kxm1 		= NDIM*(szList.at(ci) + vim1);
			kym1 		= NDIM*(szList.at(ci) + vim1) + 1;

			kx 			= NDIM*(szList.at(ci) + vi);
			ky 			= NDIM*(szList.at(ci) + vi) + 1;

			kxp1 		= NDIM*(szList.at(ci) + vip1);
			kyp1 		= NDIM*(szList.at(ci) + vip1) + 1;

			kxp2 		= NDIM*(szList.at(ci) + vip2);
			kyp2 		= NDIM*(szList.at(ci) + vip2) + 1;

			// segment length vector components
			lim2x 		= vpos[kxm1] - vpos[kxm2];
			lim2x		-= L[0]*round(lim2x/L[0]);

			lim2y 		= vpos[kym1] - vpos[kym2];
			lim2y 		-= L[1]*round(lim2y/L[1]);


			lim1x 		= vpos[kx] - vpos[kxm1];
			lim1x		-= L[0]*round(lim1x/L[0]);

			lim1y 		= vpos[ky] - vpos[kym1];
			lim1y 		-= L[1]*round(lim1y/L[1]);


			lix 		= vpos[kxp1] - vpos[kx];
			lix			-= L[0]*round(lix/L[0]);

			liy 		= vpos[kyp1] - vpos[ky];
			liy 		-= L[1]*round(liy/L[1]);


			lip1x 		= vpos[kxp2] - vpos[kxp1];
			lip1x		-= L[0]*round(lip1x/L[0]);

			lip1y 		= vpos[kyp2] - vpos[kyp1];
			lip1y 		-= L[1]*round(lip1y/L[1]);


			// segment lengths
			lim1 		= sqrt(lim1x*lim1x + lim1y*lim1y);
			li 			= sqrt(lix*lix + liy*liy);

			// segment strains
			deli 		= (li/l0tmp) - 1.0;
			delim1 		= (lim1/l0tmp) - 1.0;


			// -- PERIMETER SPRINGS

			// unit vectors
    		ulim1x   = lim1x/lim1;
    		ulim1y   = lim1y/lim1;
    
    		ulix   = lix/li;
    		uliy   = liy/li;

			// 	STIFFNESS MATRIX

			// main diagonal
		    Hl(kx,kx)       = Kl2*(ulix*ulix + ulim1x*ulim1x);
    		Hl(ky,ky)       = Kl2*(uliy*uliy + ulim1y*ulim1y);
    
    		Hl(kx,ky)       = Kl2*(ulix*uliy + ulim1x*ulim1y);
    		Hl(ky,kx)       = Hl(kx,ky);
		    
		    // 1off diagonal
		    Hl(kx,kxp1)     = -Kl2*ulix*ulix;
    		Hl(ky,kyp1)     = -Kl2*uliy*uliy;
    
    		Hl(kx,kyp1)     = -Kl2*ulix*uliy;
    		Hl(ky,kxp1)     = Hl(kx,kyp1);
		    
		    // enforce symmetry in lower triangle
		    Hl(kxp1,kx)     = Hl(kx,kxp1);
		    Hl(kyp1,ky)     = Hl(ky,kyp1);
		    
		    Hl(kyp1,kx)     = Hl(kx,kyp1);
		    Hl(kxp1,ky)     = Hl(ky,kxp1);


		    // 	STRESS MATRIX

		    // main diagonal
		    Sl(kx,kx)       = Kl1*(delim1*((ulim1y*ulim1y)/lim1) + deli*((uliy*uliy)/li));
		    Sl(ky,ky)       = Kl1*(delim1*((ulim1x*ulim1x)/lim1) + deli*((ulix*ulix)/li));
		    
		    Sl(kx,ky)       = -Kl1*(delim1*((ulim1x*ulim1y)/lim1) + deli*((ulix*uliy)/li));
		    Sl(ky,kx)       = Sl(kx,ky);
		    
		    // 1off diagonal
		    Sl(kx,kxp1)     = -Kl1*deli*((uliy*uliy)/li);
		    Sl(ky,kyp1)     = -Kl1*deli*((ulix*ulix)/li);
		    
		    Sl(kx,kyp1)     = Kl1*deli*((ulix*uliy)/li);
		    Sl(ky,kxp1)     = Sl(kx,kyp1);

		    // enforce symmetry in lower triangle
		    Sl(kxp1,kx)     = Sl(kx,kxp1);
    		Sl(kyp1,ky)     = Sl(ky,kyp1);
    
    		Sl(kyp1,kx)     = Sl(kx,kyp1);
    		Sl(kxp1,ky)     = Sl(ky,kxp1);




    		// -- CURVATURE SPRINGS


    		// curvatures
    		kapim1 			= sqrt(pow(lim1x - lim2x,2.0) + pow(lim1y - lim2y,2.0))/l0tmp;
    		kapi 			= sqrt(pow(lix - lim1x,2.0) + pow(liy - lim1y,2.0))/l0tmp;
    		kapip1 			= sqrt(pow(lip1x - lix,2.0) + pow(lip1y - liy,2.0))/l0tmp;

    		// curvature derivatives

    		// derivatives of kapim1
    		dkapim1_dxi 	= (lim1x - lim2x)/(kapim1*l0tmp*l0tmp);
    		dkapim1_dyi 	= (lim1y - lim2y)/(kapim1*l0tmp*l0tmp);

    		// derivatives of kapi
    		dkapi_dxip1 	= (lix - lim1x)/(kapi*l0tmp*l0tmp);
    		dkapi_dyip1 	= (liy - lim1y)/(kapi*l0tmp*l0tmp);
    		dkapi_dxi 		= -2.0*dkapi_dxip1;
    		dkapi_dyi 		= -2.0*dkapi_dyip1;	

    		// derivatives of kapip1
    		dkapip1_dxi 	= (lip1x - lix)/(kapip1*l0tmp*l0tmp);
    		dkapip1_dyi 	= (lip1y - liy)/(kapip1*l0tmp*l0tmp);
    		dkapip1_dxip1 	= -2.0*dkapip1_dxi;
    		dkapip1_dyip1 	= -2.0*dkapip1_dyi;
    		dkapip1_dxip2	= dkapip1_dxi;
    		dkapip1_dyip2 	= dkapip1_dyi;


    		// 	STIFFNESS MATRIX

    		// block-diagonal terms
		    Hb(kx,kx)       = Kb1*(dkapim1_dxi*dkapim1_dxi + dkapi_dxi*dkapi_dxi + dkapip1_dxi*dkapip1_dxi);
		    Hb(ky,ky)       = Kb1*(dkapim1_dyi*dkapim1_dyi + dkapi_dyi*dkapi_dyi + dkapip1_dyi*dkapip1_dyi);
		    
		    Hb(kx,ky)       = Kb1*(dkapim1_dxi*dkapim1_dyi + dkapi_dxi*dkapi_dyi + dkapip1_dxi*dkapip1_dyi);
		    Hb(ky,kx)       = Hb(kx,ky);
		    
		    // 1off block-diagonal terms
		    Hb(kx,kxp1)     = Kb1*(dkapi_dxi*dkapi_dxip1 + dkapip1_dxi*dkapip1_dxip1);
		    Hb(ky,kyp1)     = Kb1*(dkapi_dyi*dkapi_dyip1 + dkapip1_dyi*dkapip1_dyip1);
		    
		    Hb(kx,kyp1)     = Kb1*(dkapi_dxi*dkapi_dyip1 + dkapip1_dxi*dkapip1_dyip1);
		    Hb(ky,kxp1)     = Kb1*(dkapi_dyi*dkapi_dxip1 + dkapip1_dyi*dkapip1_dxip1);
		    
		    // 2off block-diagonal terms
		    Hb(kx,kxp2)     = Kb1*dkapip1_dxi*dkapip1_dxip2;
		    Hb(ky,kyp2)     = Kb1*dkapip1_dyi*dkapip1_dyip2;
		    
		    Hb(kx,kyp2)     = Kb1*dkapip1_dxi*dkapip1_dyip2;
		    Hb(ky,kxp2)     = Kb1*dkapip1_dyi*dkapip1_dxip2;
		    
		    // enforce symmetry in lower triangle
		    Hb(kxp1,kx)     = Hb(kx,kxp1);
		    Hb(kyp1,ky)     = Hb(ky,kyp1);
		    
		    Hb(kxp1,ky)     = Hb(ky,kxp1);
		    Hb(kyp1,kx)     = Hb(kx,kyp1);
		    
		    Hb(kxp2,kx)     = Hb(kx,kxp2);
		    Hb(kyp2,ky)     = Hb(ky,kyp2);
		    
		    Hb(kyp2,kx)     = Hb(kx,kyp2);
		    Hb(kxp2,ky)     = Hb(ky,kxp2);
		    
		    
		    // 	STRESS MATRIX
		    
		    // block diagonal
		    Sb(kx,kx)       = Kb2*(6.0 - (l0tmp*dkapim1_dxi)*(l0tmp*dkapim1_dxi) - (l0tmp*dkapi_dxi)*(l0tmp*dkapi_dxi) - (l0tmp*dkapip1_dxi)*(l0tmp*dkapip1_dxi));
		    Sb(ky,ky)       = Kb2*(6.0 - (l0tmp*dkapim1_dyi)*(l0tmp*dkapim1_dyi) - (l0tmp*dkapi_dyi)*(l0tmp*dkapi_dyi) - (l0tmp*dkapip1_dyi)*(l0tmp*dkapip1_dyi));
		    
		    Sb(kx,ky)       = -Kb1*(dkapim1_dxi*dkapim1_dyi + dkapi_dxi*dkapi_dyi + dkapip1_dxi*dkapip1_dyi);
		    Sb(ky,kx)       = Sb(kx,ky);
		    
		    // 1off block diagonal
		    Sb(kx,kxp1)     = -2.0*Kb2*(2.0 - (l0tmp*dkapi_dxip1)*(l0tmp*dkapi_dxip1) - (l0tmp*dkapip1_dxi)*(l0tmp*dkapip1_dxi));
		    Sb(ky,kyp1)     = -2.0*Kb2*(2.0 - (l0tmp*dkapi_dyip1)*(l0tmp*dkapi_dyip1) - (l0tmp*dkapip1_dyi)*(l0tmp*dkapip1_dyi));
		    
		    Sb(kx,kyp1)     = -Kb1*(dkapi_dxi*dkapi_dyip1 + dkapip1_dxi*dkapip1_dyip1);
		    Sb(ky,kxp1)     = -Kb1*(dkapi_dyi*dkapi_dxip1 + dkapip1_dyi*dkapip1_dxip1);
		    
		    // 2off block diagonal
		    Sb(kx,kxp2)     = Kb2*(1.0 - (l0tmp*dkapip1_dxi)*(l0tmp*dkapip1_dxi));
		    Sb(ky,kyp2)     = Kb2*(1.0 - (l0tmp*dkapip1_dyi)*(l0tmp*dkapip1_dyi));
		    
		    Sb(kx,kyp2)     = -Kb1*dkapip1_dxi*dkapip1_dyip2;
		    Sb(ky,kxp2)     = -Kb1*dkapip1_dyi*dkapip1_dxip2;
		    
		    // enforce symmetry in lower triangle
		    Sb(kxp1,kx)     = Sb(kx,kxp1);
		    Sb(kyp1,ky)     = Sb(ky,kyp1);
		    
		    Sb(kxp1,ky)     = Sb(ky,kxp1);
		    Sb(kyp1,kx)     = Sb(kx,kyp1);
		    
		    Sb(kxp2,kx)     = Sb(kx,kxp2);
		    Sb(kyp2,ky)     = Sb(ky,kyp2);
		    
		    Sb(kxp2,ky)     = Sb(ky,kxp2);
		    Sb(kyp2,kx)     = Sb(kx,kyp2);



		    // -- BELT SPRING

		    // unit vectors in belt direction
		    if (kbb > 0 && vi < hwtmp){
		    	// belt indices
		    	kxbb = kx + NDIM*hwtmp;
		    	kybb = kxbb + 1;

		    	Lbb_kx = vpos[kxbb] - vpos[kx];
		    	Lbb_kx -= L[0]*round(Lbb_kx/L[0]);

		    	Lbb_ky = vpos[kybb] - vpos[ky];
		    	Lbb_ky -= L[1]*round(Lbb_ky/L[1]);

		    	// unit vectors
		    	Lbb_k = sqrt(Lbb_kx*Lbb_kx + Lbb_ky*Lbb_ky);
		    	uLbb_kx = Lbb_kx/Lbb_k;
		    	uLbb_ky = Lbb_ky/Lbb_k;

		    	// derivatives for stiffness elements
		    	dL_dxi = -uLbb_kx/(dctmp*hwtmp);
		    	dL_dyi = -uLbb_ky/(dctmp*hwtmp);

		    	// -- stiffness elements
	        
		        // main diagonal (vi)
		        Hbb(kx,kx)          = kbb*(rho0*rho0)*dL_dxi*dL_dxi;
		        Hbb(ky,ky)          = kbb*(rho0*rho0)*dL_dyi*dL_dyi;
		        Hbb(kx,ky)          = kbb*(rho0*rho0)*dL_dxi*dL_dyi;
		        Hbb(ky,kx)          = Hbb(kx,ky);
		        
		        // off diagonals
		        Hbb(kx,kxbb)       = -Hbb(kx,kx);
		        Hbb(ky,kybb)       = -Hbb(ky,ky);
		        Hbb(kx,kybb)       = -Hbb(kx,ky);
		        Hbb(ky,kxbb)       = -Hbb(ky,kx);
		        
		        Hbb(kxbb,kx)       = Hbb(kx,kxbb);
		        Hbb(kybb,ky)       = Hbb(ky,kybb);
		        Hbb(kxbb,ky)       = Hbb(ky,kxbb);
		        Hbb(kybb,kx)       = Hbb(kx,kybb);
		        
		        // main diagonal (belt pair)
		        Hbb(kxbb,kxbb)    = Hbb(kx,kx);
		        Hbb(kybb,kybb)    = Hbb(ky,ky);
		        Hbb(kxbb,kybb)    = Hbb(kx,ky);
		        Hbb(kybb,kxbb)    = Hbb(ky,kx);
		        
		        
		        // -- stress elements
		        
		        // main diagonal (vi)
		        Sbb(kx,kx)          = ((kbb*delBB*rho0*rho0)/(hwtmp*dctmp*Lbb_k))*uLbb_ky*uLbb_ky;
		        Sbb(ky,ky)          = ((kbb*delBB*rho0*rho0)/(hwtmp*dctmp*Lbb_k))*uLbb_kx*uLbb_kx;
		        Sbb(kx,ky)          = -((kbb*delBB*rho0*rho0)/(hwtmp*dctmp*Lbb_k))*uLbb_kx*uLbb_ky;
		        Sbb(ky,kx)          = Sbb(kx,ky);
		        
		        // off diagonals
		        Sbb(kx,kxbb)       = -Sbb(kx,kx);
		        Sbb(ky,kybb)       = -Sbb(ky,ky);
		        Sbb(kx,kybb)       = -Sbb(kx,ky);
		        Sbb(ky,kxbb)       = -Sbb(ky,kx);
		        
		        Sbb(kxbb,kx)       = Sbb(kx,kxbb);
		        Sbb(kybb,ky)       = Sbb(ky,kybb);
		        Sbb(kxbb,ky)       = Sbb(ky,kxbb);
		        Sbb(kybb,kx)       = Sbb(kx,kybb);
		        
		        // main diagonal (belt pair)
		        Sbb(kxbb,kxbb)    = Sbb(kx,kx);
		        Sbb(kybb,kybb)    = Sbb(ky,ky);
		        Sbb(kxbb,kybb)    = Sbb(kx,ky);
		        Sbb(kybb,kxbb)    = Sbb(ky,kx);
		    }



    		

		    // -- AREA SPRING (stress matrix)
		    Sa(kx,kyp1) = 0.5*delA*((rho0*rho0)/a0tmp);
    		Sa(ky,kxp1) = -0.5*delA*((rho0*rho0)/a0tmp);

    		Sa(kyp1,kx) = Sa(kx,kyp1);
    		Sa(kxp1,ky) = Sa(ky,kxp1);

    		// area derivatives (for stiffness matrix)
    		da_dxi      = 0.5*(liy + lim1y);
    		da_dyi      = -0.5*(lim1x + lix);

    		// loop over other vertices, for area elasticity stiffness matrix
    		for (vj=vi; vj<nvtmp; vj++){

    			// wrap jp1 and jm1
    			vjp1 		= (vj + 1) % nvtmp;
    			vjm1 		= (vj - 1 + nvtmp) % nvtmp;

    			// dof elements
    			lxm1 		= NDIM*(szList.at(ci) + vjm1);
    			lym1 		= lxm1 + 1;

    			lx 			= NDIM*(szList.at(ci) + vj);
    			ly 			= lx + 1;

    			lxp1		= NDIM*(szList.at(ci) + vjp1);
    			lyp1		= lxp1 + 1;

    			// j segments
    			ljm1x 		= vpos[lx] - vpos[lxm1];
				ljm1x		-= L[0]*round(ljm1x/L[0]);

				ljm1y 		= vpos[ly] - vpos[lym1];
				ljm1y 		-= L[1]*round(ljm1y/L[1]);


				ljx 		= vpos[lxp1] - vpos[lx];
				ljx			-= L[0]*round(ljx/L[0]);

				ljy 		= vpos[lyp1] - vpos[ly];
				ljy 		-= L[1]*round(ljy/L[1]);

    			// area derivatives
    			da_dxj      = 0.5*(ljy + ljm1y);
    			da_dyj      = -0.5*(ljm1x + ljx);

    			// 	STIFFNESS MATRIX
    			Ha(kx,lx) = da_dxi*da_dxj*((rho0*rho0)/pow(a0tmp,2.0));
		        Ha(kx,ly) = da_dxi*da_dyj*((rho0*rho0)/pow(a0tmp,2.0));
		        
		        Ha(ky,lx) = da_dyi*da_dxj*((rho0*rho0)/pow(a0tmp,2.0));
		        Ha(ky,ly) = da_dyi*da_dyj*((rho0*rho0)/pow(a0tmp,2.0));
		        
		        Ha(lx,kx) = Ha(kx,lx);
		        Ha(ly,kx) = Ha(kx,ly);
		        
		        Ha(lx,ky) = Ha(ky,lx);
		        Ha(ly,ky) = Ha(ky,ly);

		        // also compute stiffness matrix for Belt energy if applicable
		        if (kbb > 0 && vj < hwtmp && vj > vi){
		        	// vj bb indexing
    				lxbb 		= lx + NDIM*hwtmp;
    				lybb 		= lxbb + 1;

    				// vj bb vectors
		        	Lbb_lx = vpos[lxbb] - vpos[lx];
			    	Lbb_lx -= L[0]*round(Lbb_lx/L[0]);

			    	Lbb_ly = vpos[lybb] - vpos[ly];
			    	Lbb_ly -= L[1]*round(Lbb_ly/L[1]);

			    	// unit vectors
			    	Lbb_l = sqrt(Lbb_lx*Lbb_lx + Lbb_ly*Lbb_ly);
			    	uLbb_lx = Lbb_lx/Lbb_l;
			    	uLbb_ly = Lbb_ly/Lbb_l;

			    	// derivatives for vj stiffness elements
			    	dL_dxj = -uLbb_lx/(dctmp*hwtmp);
			    	dL_dyj = -uLbb_ly/(dctmp*hwtmp);


			    	// Stiffness matrix elements (off-diagonals)
			    	// main off-diagonal block
		            Hbb(kx,lx) = kbb*(rho0*rho0)*dL_dxi*dL_dxj;
		            Hbb(ky,ly) = kbb*(rho0*rho0)*dL_dyi*dL_dyj;
		            Hbb(kx,ly) = kbb*(rho0*rho0)*dL_dxi*dL_dyj;
		            Hbb(ky,lx) = kbb*(rho0*rho0)*dL_dyi*dL_dxj;
		            
		            // enforce symmetry
		            Hbb(lx,kx) = Hbb(kx,lx);
		            Hbb(ly,ky) = Hbb(ky,ly);
		            Hbb(lx,ky) = Hbb(ky,lx);
		            Hbb(ly,kx) = Hbb(kx,ly);
		            
		            // use r + hw parity
		            Hbb(kx,lxbb) = -Hbb(kx,lx);
		            Hbb(ky,lybb) = -Hbb(ky,ly);
		            Hbb(kx,lybb) = -Hbb(kx,ly);
		            Hbb(ky,lxbb) = -Hbb(ky,lx);
		            
		            Hbb(kxbb,lx) = -Hbb(kx,lx);
		            Hbb(kybb,ly) = -Hbb(ky,ly);
		            Hbb(kxbb,ly) = -Hbb(kx,ly);
		            Hbb(kybb,lx) = -Hbb(ky,lx);
		            
		            Hbb(kxbb,lxbb) = Hbb(kx,lx);
		            Hbb(kybb,lybb) = Hbb(ky,ly);
		            Hbb(kxbb,lybb) = Hbb(kx,ly);
		            Hbb(kybb,lxbb) = Hbb(ky,lx);
		            
		            // enforce symmetry
		            Hbb(lxbb,kx) = Hbb(kx,lxbb);
		            Hbb(lybb,ky) = Hbb(ky,lybb);
		            Hbb(lxbb,ky) = Hbb(ky,lxbb);
		            Hbb(lybb,kx) = Hbb(kx,lybb);
		            
		            Hbb(lx,kxbb) = Hbb(kxbb,lx);
		            Hbb(ly,kybb) = Hbb(kybb,ly);
		            Hbb(lx,kybb) = Hbb(kybb,lx);
		            Hbb(ly,kxbb) = Hbb(kxbb,ly);
		            
		            Hbb(lxbb,kxbb) = Hbb(kxbb,lxbb);
		            Hbb(lybb,kybb) = Hbb(kybb,lybb);
		            Hbb(lxbb,kybb) = Hbb(kybb,lxbb);
		            Hbb(lybb,kxbb) = Hbb(kxbb,lybb);
		      	}
    		}
		}






		// ------------------------------------------
		// 
		// 			INTERACTION
		// 					CONTRIBUTIONS
		//
		// ------------------------------------------

		// off-diagonal components
		for (cj=ci+1; cj<NCELLS; cj++){

			// only check overlaps if contact is force-bearing, 
			// 	i.e. if both ci and cj are non-rattlers
			inContact = cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]; 

			if (inContact > 0){

				// loop over pairs of vertices on both cells, check for overlap, compute matrix elements
				for (vi=0; vi<nvtmp; vi++){

					// matrix element indices (cell ci, vertex vi)
					gi = szList.at(ci) + vi;
					mxi = NDIM*gi;
					myi = mxi + 1;

					for (vj=0; vj<nv[cj]; vj++){

						// matrix element indices (cell cj, vertex vj)
						gj = szList.at(cj) + vj;
						mxj = NDIM*gj;
						myj = mxj + 1;

						// contact distance
						sij = vrad[gi] + vrad[gj];

						// get distance between vertices
						// particle distance
						dx = vpos[mxj] - vpos[mxi];
						dx -= L[0]*round(dx/L[0]);
						if (dx < sij){
							dy = vpos[myj] - vpos[myi];
							dy -= L[1]*round(dy/L[1]);
							if (dy < sij){
								rij = sqrt(dx*dx + dy*dy);
								if (rij < sij){
									// spring constant
									kij = (eint*rho0*rho0)/(sij*rij);

									// dimensionless overlap
									h = rij/sij;

									// derivatives of distance w.r.t. coordinates
									uxij = dx/rij;
									uyij = dy/rij;

									// compute stiffness and stress matrices (off diagonal, enforce symmetry in lower triangles)

									// -- stiffness matrix
									Hvv(mxi,mxj) = -((eint*rho0*rho0)/(sij*sij))*(uxij*uxij);
									Hvv(myi,myj) = -((eint*rho0*rho0)/(sij*sij))*(uyij*uyij);
									Hvv(mxi,myj) = -((eint*rho0*rho0)/(sij*sij))*(uxij*uyij);
									Hvv(myi,mxj) = -((eint*rho0*rho0)/(sij*sij))*(uyij*uxij);

									Hvv(mxj,mxi) = Hvv(mxi,mxj);
									Hvv(myj,myi) = Hvv(myi,myj);
									Hvv(mxj,myi) = Hvv(myi,mxj);
									Hvv(myj,mxi) = Hvv(mxi,myj);



									// -- stress matrix
									Svv(mxi,mxj) = kij*(1.0 - h)*(uyij*uyij);
									Svv(myi,myj) = kij*(1.0 - h)*(uxij*uxij);
									Svv(mxi,myj) = -kij*(1.0 - h)*(uxij*uyij);
									Svv(myi,mxj) = -kij*(1.0 - h)*(uxij*uyij);

									Svv(mxj,mxi) = Svv(mxi,mxj);
					                Svv(myj,myi) = Svv(myi,myj);
					                Svv(mxj,myi) = Svv(myi,mxj);
					                Svv(myj,mxi) = Svv(mxi,myj);


					                
					                // add to diagonal, using off diagonals and reciprocity

					                // -- stiffness matrix
					                Hvv(mxi,mxi) -= Hvv(mxi,mxj);
					                Hvv(myi,myi) -= Hvv(myi,myj);
					                Hvv(mxi,myi) -= Hvv(mxi,myj);
					                Hvv(myi,mxi) -= Hvv(myi,mxj);
					                
					                Hvv(mxj,mxj) -= Hvv(mxi,mxj);
					                Hvv(myj,myj) -= Hvv(myi,myj);
					                Hvv(mxj,myj) -= Hvv(mxi,myj);
					                Hvv(myj,mxj) -= Hvv(myi,mxj);


					                // -- stress matrix
					                Svv(mxi,mxi) -= Svv(mxi,mxj);
					                Svv(myi,myi) -= Svv(myi,myj);
					                Svv(mxi,myi) -= Svv(mxi,myj);
					                Svv(myi,mxi) -= Svv(myi,mxj);
					                
					                Svv(mxj,mxj) -= Svv(mxi,mxj);
					                Svv(myj,myj) -= Svv(myi,myj);
					                Svv(mxj,myj) -= Svv(mxi,myj);
					                Svv(myj,mxj) -= Svv(myi,mxj);
								}
							}
						}
					}
				}
			}

		}	
	}

	// compute D from sum of other dynamical matrices
	// initialize all matrices to be 0 initially
	for (k=0; k<vertDOF; k++){
		for (l=0; l<vertDOF; l++){
			H(k,l) = Ha(k,l) + Hl(k,l) + Hb(k,l) + Hbb(k,l) + Hvv(k,l);
			S(k,l) = -Sa(k,l) - Sl(k,l) - Sb(k,l) - Sbb(k,l) - Svv(k,l);
			M(k,l) = H(k,l) - S(k,l);
		}
	}

	// compute eigenvalues
	cout << "\t** Computing eigenvalues and eigenvectors of M, H and S matrices" << endl;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> allModes(M);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> hModes(H);

	// // define eigenvector matrix
	Eigen::MatrixXd evecs = allModes.eigenvectors();

	// print eigenvalues to vdos file
	cout << "\t** Printing evals to file" << endl;
	vdosout << vertDOF << endl;
	vdosout << allModes.eigenvalues() << endl;
	// vdosout << sModes.eigenvalues() << endl;

	// computing projections onto H
	// cout << "\t** Computing eigenvector projections on stiffness matrix H" << endl;
	// double evecProj = 0.0;
	// for (i=0; i<vertDOF; i++){
	// 	evecProj = 0.0;
	// 	for (k=0; k<vertDOF; k++){
	// 		evecProj += M(k,k)*evecs(k,i)*evecs(k,i);
	// 		for (l=(k+1); l<vertDOF; l++){
	// 			evecProj += 2.0*M(k,l)*evecs(k,i)*evecs(l,i);
	// 		}
	// 	}
	// 	vdosout << evecProj << endl;
	// }
	vdosout << hModes.eigenvalues() << endl;

	// compute vertex participation ratio of each mode
	double pnum, pdenom;
	for (k=0; k<vertDOF; k++){
		// p.r. denominator ONLY, numerator always = 1
		pdenom = 0.0;
		for (l=0; l<(vertDOF-1); l += 2){
			pdenom += pow(evecs(l,k)*evecs(l,k) + evecs(l+1,k)*evecs(l+1,k),2.0);
		}
		pdenom *= NVTOT;

		// print
		vdosout << 1.0/pdenom << endl;
	}


	// close open objects
	posout.close();
	vdosout.close();


	// print to console, return
	cout << "\n\n\nFINISHED MAIN FOR bidRepulsiveCellJamming.cpp, ENDING." << endl << endl << endl;
	return 0;
}








/* 

	&&&&&&&&&&&&&&&&&&&&&&& FUNCTION DEFINITIONS &&&&&&&&&&&&&&&&&&&&&

	FUNCTIONS DEFINED

	gindex 			: returns global vertex index (gi) given cell (ci) and local vertex index (vi)
	cindex 			: returns cell index (ci) given global vertex index (gi)

	area 			: returns area of cell ci
	perimeter 		: returns perimeter of cell ci

	removeRattlers	: remove all rattlers from a contact network

	printPos 		: output vertex positions to .pos file for processing and visualization

	&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 

*/



// -- INDEXING


// get global vertex index gi given input cell index ci and vertex index vi
int gindex(int ci, int vi, vector<int>& szList){
	return szList[ci] + vi;
} 


// get cell index ci and vertex index 
void cindices(int& ci, int& vi, int gi, int NCELLS, vector<int>& szList){
	if (gi >= szList[NCELLS-1]){
		ci = NCELLS - 1;
		vi = gi - szList[NCELLS-1];
	}
	else{
		for (int i=1; i<NCELLS; i++){
			if (szList[i] > gi){
				ci = i-1;
				vi = gi - szList[i-1];
				break;
			}
		}
	}
}








// -- CELL SHAPE


// get cell area
double area(vector<double>& vpos, int ci, vector<double>& L, vector<int>& nv, vector<int>& szList){
	// local variables
	int vi, vip1, gi, gip1;
	double dx, dy, xi, yi, xip1, yip1, areaVal = 0.0;

	// initial position: vi = 0
	gi = gindex(ci,0,szList);
	xi = vpos[NDIM*gi];
	yi = vpos[NDIM*gi + 1];

	// loop over vertices of cell ci, get area by shoe-string method
	for (vi=0; vi<nv.at(ci); vi++){
		// next vertex
		vip1 = (vi + 1) % nv.at(ci);
		gip1 = gindex(ci,vip1,szList);

		// get positions (check minimum images)
		dx 		= vpos[NDIM*gip1] - xi;
		dx 		-= L[0]*round(dx/L[0]);
		xip1 	= xi + dx;

		dy 		= vpos[NDIM*gip1 + 1] - yi;
		dy 		-= L[1]*round(dy/L[1]);
		yip1 	= yi + dy;

		// increment area
		areaVal += xi*yip1 - xip1*yi;

		// set next coordinates
		xi = xip1;
		yi = yip1;
	}
	areaVal *= 0.5;

	return abs(areaVal);
}


// get cell perimeter
double perimeter(vector<double>& vpos, int ci, vector<double>& L, vector<int>& nv, vector<int>& szList){
		// local variables
	int vi, vip1, gi, gip1;
	double dx, dy, xi, yi, xip1, yip1, l, perimVal = 0.0;

	// initial position: vi = 0
	gi = gindex(ci,0,szList);
	xi = vpos[NDIM*gi];
	yi = vpos[NDIM*gi + 1];

	// loop over vertices of cell ci, get perimeter
	for (vi=0; vi<nv.at(ci); vi++){
		// next vertex
		vip1 = (vi + 1) % nv.at(ci);
		gip1 = gindex(ci,vip1,szList);

		// get positions (check minimum images)
		dx 		= vpos[NDIM*gip1] - xi;
		dx 		-= L[0]*round(dx/L[0]);
		xip1 	= xi + dx;

		dy 		= vpos[NDIM*gip1 + 1] - yi;
		dy 		-= L[1]*round(dy/L[1]);
		yip1 	= yi + dy;

		// compute segment length
		l = sqrt(dx*dx + dy*dy);

		// add to perimeter
		perimVal += l;

		// update coordinates
		xi = xip1;
		yi = yip1;
	}

	// return perimeter
	return perimVal;
}






// RECURSIVELY REMOVE RATTLERS FROM CONTACT NETWORK

int removeRattlers(vector<int>& cij){
	// local variables
	int NCTCS, NCELLS, ci, cj, ctmp, rvv, rcc, nr, nm;

	// total size of contact space
	NCTCS = cij.size();

	// get NCELLS from contact network
	NCELLS = (1 + sqrt(1 + 8*NCTCS))/2;

	// number of rattlers
	nr = 0;

	// number of "marginal" rattlers to be removed
	nm = 0;

	// loop over rows, eliminate contacts to rattlers
	for (ci=0; ci<NCELLS; ci++) {
		// get number of contacts on cell ci
		rvv = 0;
		rcc = 0;
		for (cj=0; cj<NCELLS; cj++){
			if (ci != cj){
				if (ci > cj)
					ctmp = cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2];
				else
					ctmp = cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2];
			}
			else
				ctmp = 0;
			

			rvv += ctmp;
			if (ctmp > 0)
				rcc++;
		}

		// check to see if particle should be removed from network
		if (rcc <= NDIM && rvv <= 3) {
			// increment # of rattlers
			nr++;

			// if in contact, remove contacts
			if (rvv > 0) {
				nm++;

				for (cj=0; cj<NCELLS; cj++) {
					// delete contact between ci and cj
					if (ci != cj){
						if (ci > cj)
							cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2] = 0;
						else
							cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2] = 0; 
					}
				}
			}
		}
	}

	// recursively check once rattlers are removed
	if (nm == 0)
		return nr;
	else
		return removeRattlers(cij);
}




// -- PRINT TO FILE


// print cell positions
void printPos(ofstream& posout, vector<double>& vpos, vector<double>& a0, vector<double>& l0, vector<double>& L, vector<int>& cij, vector<int>& nv, vector<int>& szList, double phi, int NCELLS){
	// local variables
	int ci, cj, vi, gi, ctmp, zc, zv;
	double xi, yi, dx, dy, Lx, Ly;

	// save box sizes
	Lx = L.at(0);
	Ly = L.at(1);

	// print information starting information
	posout << setw(w) << left << "NEWFR" << " " << endl;
	posout << setw(w) << left << "NUMCL" << setw(w) << right << NCELLS << endl;
	posout << setw(w) << left << "PACKF" << setw(wnum) << setprecision(pnum) << right << phi << endl;

	// print box sizes
	posout << setw(w) << left << "BOXSZ";
	posout << setw(wnum) << setprecision(pnum) << right << Lx;
	posout << setw(wnum) << setprecision(pnum) << right << Ly;
	posout << endl;

	// print coordinate for rest of the cells
	for (ci=0; ci<NCELLS; ci++){

		// get cell contact data
		zc = 0;
		zv = 0;
		for (cj=0; cj<NCELLS; cj++){
			if (ci != cj){
				// grab contact info from entry ci, cj
				ctmp = 0;
				if (ci > cj)
					ctmp = cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2];
				else
					ctmp = cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]; 

				// add to contact information
				zv += ctmp;
				if (ctmp > 0)
					zc++;
			}
		}

		// cell information
		posout << setw(w) << left << "CINFO";
		posout << setw(w) << right << nv.at(ci);
		posout << setw(w) << right << zc;
		posout << setw(w) << right << zv;
		posout << setw(wnum) << right << a0.at(ci);
		posout << setw(wnum) << right << l0.at(ci);
		posout << endl;

		// get initial vertex positions
		gi = gindex(ci,0,szList);
		xi = vpos.at(NDIM*gi);
		yi = vpos.at(NDIM*gi + 1);

		// place back in box center
		xi = fmod(xi,Lx);
		yi = fmod(yi,Ly);

		posout << setw(w) << left << "VINFO";
		posout << setw(w) << left << ci;
		posout << setw(w) << left << 0;

		// output initial vertex information
		posout << setw(wnum) << setprecision(pnum) << right << xi;
		posout << setw(wnum) << setprecision(pnum) << right << yi;
		posout << endl;

		// vertex information for next vertices
		for (vi=1; vi<nv.at(ci); vi++){
			// get global vertex index for next vertex
			gi++;

			// get next vertex positions (use MIC)
			dx = vpos.at(NDIM*gi) - xi;
			dx -= Lx*round(dx/Lx);
			xi += dx;

			dy = vpos.at(NDIM*gi + 1) - yi;
			dy -= Ly*round(dy/Ly);
			yi += dy;

			// Print indexing information
			posout << setw(w) << left << "VINFO";
			posout << setw(w) << left << ci;
			posout << setw(w) << left << vi;

			// output vertex information
			posout << setw(wnum) << setprecision(pnum) << right << xi;
			posout << setw(wnum) << setprecision(pnum) << right << yi;
			posout << endl;
		}
	}

	// print end frame
	posout << setw(w) << left << "ENDFR" << " " << endl;
}





