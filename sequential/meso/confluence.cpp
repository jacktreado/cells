/* 

	MAIN FILE FOR 2D CONFLUENCE CODE

		-- Compresses packing to arbitrary phi = sum(a0)/L^2
		-- Prints position / shape info every dphiPrint
		-- Supports bending energy (kb) and vertex-vertex attraction

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

# define NDIM 2
# define NNN 4

// namespace
using namespace std;

// GLOBAL CONSTANTS
const double PI 			= 4*atan(1);
const int w 				= 10;
const int wnum 				= 25;
const int pnum 				= 14;

// simulation constants
const int nvmin 			= 12;
const double timeStepMag 	= 0.005;
const double dphiGrow 		= 0.001;
const double dphiPrint 		= 0.01;
const double sizeRatio 		= 1.4;
const double sizeFraction 	= 0.5;

// FIRE constants for initial minimizations (SP + DP)
const double alpha0      	= 0.2;
const double finc        	= 1.1;
const double fdec        	= 0.5;
const double falpha      	= 0.99;
const double Ftol 			= 1e-8;

const int NSKIP 			= 2e4;
const int NMIN        		= 10;
const int NNEGMAX     		= 1000;
const int NDELAY      		= 10;
const int itmax       		= 5e7;

// DP force constants
const double ka 			= 1.0;			// area spring (should be = 1)
const double eint 			= 1.0;			// baseline interaction energy 
const double del 			= 1.0;			// radius of vertices in units of l0

// FUNCTION PROTOTYPES

// indexing
int gindex(int ci, int vi, vector<int>& szList);
void cindices(int& ci, int& vi, int gi, int NCELLS, vector<int>& szList);

// particle geometry
double area(vector<double>& vpos, int ci, vector<double>& L, vector<int>& nv, vector<int>& szList);
double perimeter(vector<double>& vpos, int ci, vector<double>& L, vector<int>& nv, vector<int>& szList);

// print to file
void printPos(ofstream& posout, vector<double>& vpos, vector<double>& vrad, vector<double>& a0, vector<double>& calA0, vector<double>& L, vector<int>& cij, vector<int>& nv, vector<int>& szList, double phi, int NCELLS); 



// MAIN
int main(int argc, char const *argv[]){
	// variables for indexing loops
	int i, ci, cj, vi, vj, gi, gj, d;

	// parameters to be read in 
	int NCELLS, largeNV, smallNV, smallN, largeN, NVTOT, cellDOF, vertDOF, seed;
	double calA0Input, phi, phiMax, phiMin, kl, kb, att;

	// read in parameters from command line input
	// test: g++ -O3 sequential/meso/confluence.cpp -o conf.o
	// test: ./conf.o 12 20 1.08 0.2 1.0 1.0 0 0 1 pos.test shape.test
	// 
	// PARAMETERS:
	// 1. NCELLS 		= # of dpm particles
	// 2. NV 			= # of vertices on smaller cells, larger cells = round(NV*1.4)
	// 3. calA0 		= shape parameter of all cells, initial a0 = 1
	// 4. phiMin 		= min packing fraction (start)
	// 5. phiMax 		= max packing fraction (end)
	// 6. kl 			= perimeter spring constant
	// 7. kb 			= bending spring constant
	// 8. att 			= short-range attraction parameter
	// 9. seed 			= random number generation seed
	// 10. positionFile	= position data file string
	// 11. shapeFile	= shape data file string 

	string NCELLS_str 			= argv[1];
	string NV_str 				= argv[2];
	string calA0_str 			= argv[3];
	string phiMin_str 			= argv[4];
	string phiMax_str  			= argv[5];
	string kl_str 				= argv[6];
	string kb_str 				= argv[7];
	string att_str 				= argv[8];
	string seed_str 			= argv[9];
	string positionFile 		= argv[10];
	string shapeFile 			= argv[11];

	stringstream NCELLSss(NCELLS_str);
	stringstream NVss(NV_str);
	stringstream calA0ss(calA0_str);
	stringstream phiMaxss(phiMax_str);
	stringstream phiMinss(phiMin_str);
	stringstream klss(kl_str);
	stringstream kbss(kb_str);
	stringstream attss(att_str);
	stringstream seedss(seed_str);

	NCELLSss >> NCELLS;
	NVss >> smallNV;
	calA0ss >> calA0Input;
	phiMaxss >> phiMax;
	phiMinss >> phiMin;
	klss >> kl;
	kbss >> kb;
	attss >> att;
	seedss >> seed;

	// seed random number generator
	srand48(seed);

	// open xyz file
	ofstream posout;
	posout.open(positionFile.c_str());
	if (!posout.is_open()){
		cout << "	** ERROR: position file " << positionFile << " could not be opened, ending." << endl;
		return 1;
	}

	ofstream shapeout;
	shapeout.open(shapeFile.c_str());
	if (!shapeout.is_open()){
		cout << "	** ERROR: shape file " << shapeFile << " could not be opened, ending." << endl;
		return 1;
	}

	// total number of vertices
	largeNV = round(sizeRatio*smallNV);

	// total number of vertices
	smallN 	= round(sizeFraction*NCELLS);
	largeN 	= NCELLS - smallN;
	NVTOT 	= smallNV*smallN + largeNV*largeN;

	// szList and nv (keep track of global vertex indices)
	cout << "makin vectors" << endl;
	vector<int> szList(NCELLS,0);
	vector<int> nv(NCELLS,0);
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

	// fundamental MD time units
	double dtMD, dt0, dt;

	dtMD 	= 1.0;
	dt0 	= timeStepMag*dtMD;
	dt 		= dt0;

	// output opening statement to console
	cout << "=======================================================" << endl << endl;

	cout << "		confluence.cpp 								" << endl;
	cout << "		Jack Treado, 2021   							" << endl;

	cout << "		NCELLS 		= " << NCELLS << "					" << endl;
	cout << "		NV (small)	= " << smallNV << "						" << endl;
	cout << "		NV (large)	= " << largeNV << "						" << endl;
	cout << "		NVTOT 		= " << NVTOT << "					" << endl;
	cout << endl;

	cout << "		calA0 		= " << calA0Input << "				" << endl;

	cout << "		phiMin 		= " << phiMin << "					" << endl << endl;
	cout << "		phiMax 		= " << phiMax << " 					" << endl;

	cout << "		kl 			= " << kl << "						" << endl;
	cout << "		kb 			= " << kb << "						" << endl;
	cout << "		att 		= " << att << " 					" << endl;
	cout << "		seed 		= " << seed << "					" << endl << endl;

	cout << "		pos file 	= " << positionFile << "			" << endl;
	cout << "		shape file 	= " << shapeFile << "				" << endl << endl;
	
	cout << "=======================================================" << endl << endl;


















	/* * * * * * * * * * * * * * * * * * 

				INITIALIZE

					PARTICLES

	 * * * * * * * * * * * * * * * * * */

	// initialization variables
	int nvtmp;
	double a0tmp, lenscale, calA0tmp, areaSum = 0.0;

	// initialize vectors for storing coordinates, shape information
	vector<double> vrad(NVTOT,1.0);
	vector<double> drad(NCELLS,1.0);

	vector<double> vpos(vertDOF,0.0);
	vector<double> dpos(cellDOF,0.0);

	vector<double> a0(NCELLS,1.0);
	vector<double> l0(NCELLS,1.0);
	vector<double> calA0(NCELLS,1.0);

	// shape parameters
	double smallCalA0 = calA0Input*smallNV*tan(PI/smallNV)/PI;
	double largeCalA0 = calA0Input*largeNV*tan(PI/largeNV)/PI;

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
		calA0.at(ci) = calA0tmp;

		// store preferred area
		a0.at(ci) 		= a0tmp;

		// set disk radius
		drad.at(ci) 	= 1.1*sqrt((2.0*a0tmp)/(nvtmp*sin(2.0*PI/nvtmp)));

		// set l0, vector radius
		l0.at(ci) 	= 2.0*lenscale*sqrt(PI*calA0tmp)/nvtmp;
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
		L.at(d) = sqrt(areaSum/phiMin);
	phi = phiMin;

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
	double boxLengthScale = 3.0;
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
		fcheck = fnorm/(NDIM*NCELLS);

		// update npPMin
		if (fcheck < Ftol && fireit > NDELAY)
			npPMin++;
		else
			npPMin = 0;

		// print to console
		if (fireit % NSKIP == 0){
			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << "		I N I T I A L 				" << endl;
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

			INITIAL DENSE STATE

	 * * * * * * * * * * * * * * * * * */

	// initialize disk velocity and force vectors
	vector<double> vvel(vertDOF,0.0);
	vector<double> vF(vertDOF,0.0);
	vector<double> vFold(vertDOF,0.0);

	// jamming check variables
	int k, kmax, xind, yind;
	double pcheck;

	// max number of jamming iteractions
	k = 0;
	kmax = 1e4;

	// compute scale factors
	double scaleFactor = sqrt((phiMin + dphiGrow)/phiMin);

	// linked list variables
	int boxid, bi, bj, pi, pj, sbtmp;

	// total potential energy
	double U = 0.0;

	// length unit variable
	double rho0 = 0.0;

	// shape force variables
	double fa, fl, fb, l0tmp, atmp, li, lim1, cx, cy;
	double da, dli, dlim1;
	double lim2x, lim2y, lim1x, lim1y, lix, liy, lip1x, lip1y;
	double rim2x, rim2y, rim1x, rim1y, rix, riy, rip1x, rip1y, rip2x, rip2y;
	double cutij, shellij;
	double lastPrintPhi = phi;

	// print frame info
	int printF = 0;

	// compress to jamming, relax U and F using FIRE
	while (phi < phiMax && k < kmax){
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

		// RESET FIRE VARIABLES
		P  			= 0;	
		fnorm 		= 0;
		vnorm 		= 0;
		alpha   	= alpha0;

		dtmax   	= 5.0*dt0;
		dtmin   	= 1e-8*dt0;
		dt 			= dt0;

		npPos      	= 0;
		npNeg      	= 0;
		npPMin      = 0;

		fireit    	= 0;
		fcheck  	= 10*Ftol;

		// reset forces
		for (i=0; i<vertDOF; i++){
			vF[i] = 0.0;
			vFold[i] = 0.0;
			vvel[i] = 0.0;
		}

		// set length constant (variable due to particle growth)
		rho0 = sqrt(a0.at(0));

		// RELAX FORCES USING FIRE
		while ((fcheck > Ftol || npPMin < NMIN) && fireit < itmax){
			// VV POSITION UPDATE
			for (i=0; i<vertDOF; i++){
				// update position
				vpos[i] += dt*vvel[i] + 0.5*dt*dt*vF[i];

				// recenter in box
				if (vpos[i] > L[i % NDIM])
					vpos[i] -= L[i % NDIM];
				else if (vpos[i] < 0)
					vpos[i] += L[i % NDIM];

				// reset forces
				vF[i] = 0;
			}

			// reset linked list variables
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

			// reset contact networks
			fill(cij.begin(), cij.end(), 0);

			// FORCE UPDATE

			// interaction forces (USE BOX LINKED LIST)
			pcheck = 0.0;
			U = 0.0;
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

						// attraction distances
						shellij = (1.0 + att)*sij;
						cutij = (1.0 + 0.5*att)*sij;

						// particle distance
						dx = vpos[NDIM*gj] - vpos[NDIM*gi];
						dx -= L[0]*round(dx/L[0]);
						if (dx < shellij){
							dy = vpos[NDIM*gj + 1] - vpos[NDIM*gi + 1];
							dy -= L[1]*round(dy/L[1]);
							if (dy < shellij){
								rij = sqrt(dx*dx + dy*dy);
								// if two tumor cells, check for attraction
								if (rij < shellij && rij > cutij && att > 0.0){
									// force scale
									ftmp 				= eint*((rij/sij) - 1.0 - att)/sij;
									fx 					= ftmp*(dx/rij);
									fy 					= ftmp*(dy/rij);

									// add to forces
									vF[NDIM*gi] 		-= fx;
									vF[NDIM*gi + 1] 	-= fy;

									vF[NDIM*gj] 		+= fx;
									vF[NDIM*gj + 1] 	+= fy;

									// add to virial expression for pressure
									pcheck += dx*fx + dy*fy;

									// add to potential energy
									U -= 0.5*eint*pow((rij/sij) - 1.0 - att,2.0);

									// add to contacts
									cindices(ci, vi, gi, NCELLS, szList);
									cindices(cj, vj, gj, NCELLS, szList);

									if (ci > cj)
										cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
									else if (ci < cj)
										cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++;  
								}
								else if (rij < cutij && rij > sij && att > 0.0){
									// force scale
									ftmp 				= eint*(1 - (rij/sij))/sij;
									fx 					= ftmp*(dx/rij);
									fy 					= ftmp*(dy/rij);

									// add to forces
									vF[NDIM*gi] 		-= fx;
									vF[NDIM*gi + 1] 	-= fy;

									vF[NDIM*gj] 		+= fx;
									vF[NDIM*gj + 1] 	+= fy;

									// add to virial expression for pressure
									pcheck += dx*fx + dy*fy;

									// add to potential energy
									U += 0.5*eint*(pow(1.0 - (rij/sij),2.0) - 0.5*att*att);

									// add to contacts
									cindices(ci, vi, gi, NCELLS, szList);
									cindices(cj, vj, gj, NCELLS, szList);

									if (ci > cj)
										cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
									else if (ci < cj)
										cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++;  
								}
								// otherwise, check for two overlapping purely-repulsive vertices
								else if (rij < sij){
									// force scale
									ftmp 				= eint*(1 - (rij/sij))/sij;
									fx 					= ftmp*(dx/rij);
									fy 					= ftmp*(dy/rij);

									// add to forces
									vF[NDIM*gi] 		-= fx;
									vF[NDIM*gi + 1] 	-= fy;

									vF[NDIM*gj] 		+= fx;
									vF[NDIM*gj + 1] 	+= fy;

									// add to virial expression for pressure
									pcheck += dx*fx + dy*fy;

									// add to potential energy
									U += 0.5*eint*(pow(1.0 - (rij/sij),2.0) - 0.5*att*att);

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

							// attraction distances
							shellij = (1.0 + att)*sij;
							cutij = (1.0 + 0.5*att)*sij;

							// particle distance
							dx = vpos[NDIM*gj] - vpos[NDIM*gi];
							dx -= L[0]*round(dx/L[0]);
							if (dx < shellij){
								dy = vpos[NDIM*gj + 1] - vpos[NDIM*gi + 1];
								dy -= L[1]*round(dy/L[1]);
								if (dy < shellij){
									rij = sqrt(dx*dx + dy*dy);
									// if two tumor cells, check for attraction
									if (rij < shellij && rij > cutij && att > 0.0){
										// force scale
										ftmp 				= eint*((rij/sij) - 1.0 - att)/sij;
										fx 					= ftmp*(dx/rij);
										fy 					= ftmp*(dy/rij);

										// add to forces
										vF[NDIM*gi] 		-= fx;
										vF[NDIM*gi + 1] 	-= fy;

										vF[NDIM*gj] 		+= fx;
										vF[NDIM*gj + 1] 	+= fy;

										// add to virial expression for pressure
										pcheck += dx*fx + dy*fy;

										// add to potential energy
										U -= 0.5*eint*pow((rij/sij) - 1.0 - att,2.0);

										// add to contacts
										cindices(ci, vi, gi, NCELLS, szList);
										cindices(cj, vj, gj, NCELLS, szList);

										if (ci > cj)
											cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
										else if (ci < cj)
											cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++;  
									}
									else if (rij < cutij && rij > sij && att > 0.0){
										// force scale
										ftmp 				= eint*(1 - (rij/sij))/sij;
										fx 					= ftmp*(dx/rij);
										fy 					= ftmp*(dy/rij);

										// add to forces
										vF[NDIM*gi] 		-= fx;
										vF[NDIM*gi + 1] 	-= fy;

										vF[NDIM*gj] 		+= fx;
										vF[NDIM*gj + 1] 	+= fy;

										// add to virial expression for pressure
										pcheck += dx*fx + dy*fy;

										// add to potential energy
										U += 0.5*eint*(pow(1.0 - (rij/sij),2.0) - 0.5*att*att);
										
										// add to contacts
										cindices(ci, vi, gi, NCELLS, szList);
										cindices(cj, vj, gj, NCELLS, szList);

										if (ci > cj)
											cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
										else if (ci < cj)
											cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++;  
									}
									// otherwise, check for two overlapping purely-repulsive vertices
									else if (rij < sij){
										// force scale
										ftmp 				= eint*(1 - (rij/sij))/sij;
										fx 					= ftmp*(dx/rij);
										fy 					= ftmp*(dy/rij);

										// add to forces
										vF[NDIM*gi] 		-= fx;
										vF[NDIM*gi + 1] 	-= fy;

										vF[NDIM*gj] 		+= fx;
										vF[NDIM*gj + 1] 	+= fy;

										// add to virial expression for pressure
										pcheck += dx*fx + dy*fy;

										// add to potential energy
										U += 0.5*eint*(pow(1.0 - (rij/sij),2.0) - 0.5*att*att);

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
						a0tmp = a0[ci];
						l0tmp = l0[ci];

						// compute area deviation
						atmp = area(vpos,ci,L,nv,szList);
						da = (atmp/a0tmp) - 1.0;

						// add to potential energy
						U += 0.5*pow(da,2.0);

						// shape force parameters (kl and kl are unitless energy ratios)
						fa = da*(rho0/a0tmp);		// derivation from the fact that rho0^2 does not necessarily cancel a0tmp
						fl = kl*(rho0/l0tmp);
						fb = kb*(rho0/(l0tmp*l0tmp));
						
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

				// add to potential energy
				U += 0.5*kl*pow(dli,2.0);

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

					// add to potential energy
					U += 0.5*kb*(pow(lix - lim1x,2.0) + pow(liy - lim1y,2.0))/(l0tmp*l0tmp);

					// add to force
					vF[NDIM*gi] 		+= fb*(3.0*(lix - lim1x) + lim2x - lip1x);
					vF[NDIM*gi + 1] 	+= fb*(3.0*(liy - lim1y) + lim2y - lip1y);
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
			for (i=0; i<vertDOF; i++){
				vvel[i] += 0.5*(vF[i] + vFold[i])*dt;
				vFold[i] = vF[i];
			}



			// FIRE UPDATE
			// compute fnorm, vnorm and P
			fnorm = 0.0;
			vnorm = 0.0;
			P = 0.0;
			for (i=0; i<vertDOF; i++){
				fnorm 	+= vF[i]*vF[i];
				vnorm 	+= vvel[i]*vvel[i];
				P 		+= vvel[i]*vF[i];
			}
			fnorm = sqrt(fnorm);
			vnorm = sqrt(vnorm);

			// update fcheck based on fnorm (= force per degree of freedom)
			fcheck = fnorm/sqrt(NCELLS);

			// update npPMin
			if (fcheck < Ftol)
				npPMin++;
			else
				npPMin = 0;

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

				cout << "	** vnorm 	= " << vnorm << endl;
				cout << "	** dt 		= " << dt << endl;
				cout << "	** P 		= " << P << endl;
				cout << "	** Pdir 	= " << P/(fnorm*vnorm) << endl;
				cout << "	** alpha 	= " << alpha << endl;
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
					vpos[i] -= 0.5*dt*vvel[i] + 0.25*dt*dt*vF[i];

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


			// update velocities (s.d. vs inertial dynamics) only if forces are acting
			if (fnorm > 0){
				for (i=0; i<vertDOF; i++)
					vvel[i] = (1 - alpha)*vvel[i] + alpha*(vnorm/fnorm)*vF[i];
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
		cout << "	* dr 			= " << scaleFactor << endl;
		cout << "	* phi 			= " << phi << endl;
		cout << "	* fcheck 		= " << fcheck << endl;
		cout << "	* pcheck 		= " << pcheck << endl;
		cout << "	* U 		 	= " << U << endl;
		cout << endl;

		// print position, energy + shape information
		if (abs(lastPrintPhi - phi) > dphiPrint){
			// print shape info
			cout << "\t** PRINTING SHAPE INFO TO FILE... " << endl;
			shapeout << setw(6) << printF;
			shapeout << setw(15) << fireit;
			shapeout << setw(15) << phi;
			shapeout << setw(15) << pcheck;
			shapeout << setw(15) << U;
			for (ci=0; ci<NCELLS; ci++){
				shapeout << setw(15) << perimeter(vpos,ci,L,nv,szList);
				shapeout << setw(15) << area(vpos,ci,L,nv,szList);
			}
			shapeout << endl;

			// print position info
			cout << "\t** PRINTING POSITIONS TO FILE... " << endl;
			printPos(posout, vpos, vrad, a0, calA0, L, cij, nv, szList, phi, NCELLS);

			// update last phi when printed, for next time
			lastPrintPhi = phi;
			printF++;
		}

		// grow or shrink particles by scale factor
		phi = 0.0;
		for (ci=0; ci<NCELLS; ci++){
			// scale preferred lengths
			l0[ci] *= scaleFactor;
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
			phi += a0[ci] + 0.25*PI*pow(l0[ci]*del,2.0)*(0.5*nv[ci] - 1);
		}
		phi /= L[0]*L[1];
		scaleFactor = sqrt((phi + dphiGrow)/phi);
	}
	if (k == kmax){
		cout << "	** ERROR: IN 2d cell jamming finding, k reached kmax without finding jamming. Ending." << endl;
		return 1;
	}

	// print shape info
	cout << "\t** PRINTING SHAPE INFO TO FILE... " << endl;
	shapeout << setw(6) << printF;
	shapeout << setw(15) << fireit;
	shapeout << setw(15) << phi;
	shapeout << setw(15) << pcheck;
	shapeout << setw(15) << U;
	for (ci=0; ci<NCELLS; ci++){
		shapeout << setw(15) << perimeter(vpos,ci,L,nv,szList);
		shapeout << setw(15) << area(vpos,ci,L,nv,szList);
	}
	shapeout << endl;

	// print position info
	cout << "\t** PRINTING POSITIONS TO FILE... " << endl;
	printPos(posout, vpos, vrad, a0, calA0, L, cij, nv, szList, phi, NCELLS);

	// close open objects
	shapeout.close();
	posout.close();

	return 0;
}



/* 

	&&&&&&&&&&&&&&&&&&&&&&& FUNCTION DEFINITIONS &&&&&&&&&&&&&&&&&&&&&

	FUNCTIONS DEFINED

	gindex 			: returns global vertex index (gi) given cell (ci) and local vertex index (vi)
	cindex 			: returns cell index (ci) given global vertex index (gi)

	area 			: returns area of cell ci
	perimeter 		: returns perimeter of cell ci

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










// -- PRINT TO FILE


// print cell positions
void printPos(ofstream& posout, vector<double>& vpos, vector<double>& vrad, vector<double>& a0, vector<double>& calA0, vector<double>& L, vector<int>& cij, vector<int>& nv, vector<int>& szList, double phi, int NCELLS){
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
		posout << setw(wnum) << right << calA0.at(ci);
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
		posout << setw(wnum) << setprecision(pnum) << right << vrad.at(gi);
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
			posout << setw(wnum) << setprecision(pnum) << right << vrad.at(gi);
			posout << setw(wnum) << setprecision(pnum) << right << xi;
			posout << setw(wnum) << setprecision(pnum) << right << yi;
			posout << endl;
		}
	}

	// print end frame
	posout << setw(w) << left << "ENDFR" << " " << endl;
}




