
/* 

	MAIN FILE FOR TUMOR INVASION INTO
		A MODEL ADIPOSE TISSUE

	** 	FIRST tN cells are tumor cells
		NEXT NCELLS - tN cells are adipocites
		Total number of cells = NCELLS

	Jack Treado
	10/20/2020, in the time of covid

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
const double timeStepMag 	= 0.01;
const double sizeRatio 		= 3.0;
const double phiT 			= 0.95;
const double dphi0 			= 5e-3;

// FIRE constants for initial minimizations (SP + DP)
const double alpha0      	= 0.2;
const double finc        	= 1.1;
const double fdec        	= 0.5;
const double falpha      	= 0.99;
const double Ftol0 			= 1e-6;

const int NSKIP 			= 1e3;
const int NMIN        		= 10;
const int NNEGMAX     		= 1000;
const int NDELAY      		= 50;
const int itmax       		= 5e7;

const int NACTIVESKIP 		= 2e3;

// DP force constants
const double ka 			= 1.0;			// area spring (should be = 1)
const double kl 			= 0.01;			// contractility
const double eint 			= 1.0;			// interaction energy 
const double del 			= 1.0;			// radius of vertices in units of l0

// active tumor cell
const double v0 			= 0.1;			// tumor velocity velocity (in units of sqrt(a0)/tau)
const double Ds 			= 0.1;			// spread of velocity coupling along tumor cell boundary
const double vmin 			= 1e-2*v0;		// minimum velocity of any given vertex along tc boundary 


// FUNCTION PROTOTYPES
int gindex(int ci, int vi, vector<int>& szList);
void cindices(int& ci, int& vi, int gi, int NCELLS, vector<int>& szList);

double area(vector<double>& dpos, int ci, vector<double>& L, vector<int>& nv, vector<int>& szList);
double perimeter(vector<double>& dpos, int ci, vector<double>& L, vector<int>& nv, vector<int>& szList);

// print to file
void printPos(ofstream& posout, vector<double>& vpos, vector<double>& a0, vector<double>& l0, vector<double>& L, vector<int>& nv, vector<int>& szList, double phi, int NCELLS); 


// MAIN
int main(int argc, char const *argv[]){
	// variables for indexing loops
	int i, ci, cj, vi, vj, gi, gj, d;

	// parameters to be read in 
	int NCELLS, NT, NV, tN, aN, NVTOT, cellDOF, vertDOF, seed;
	double NT_dbl, calA0, tCalA0, aCalA0, gam, Dr, phi0, Ftol;

	// set initial force tolerance
	Ftol = Ftol0;

	// read in parameters from command line input
	string NCELLS_str 			= argv[1];
	string NV_str 				= argv[2];
	string tN_str  				= argv[3];
	string tumorCalA0_str 		= argv[4];
	string adiposeCalA0_str 	= argv[5];
	string gam_str 				= argv[6];
	string phi0_str 			= argv[7];
	string Dr_str 				= argv[8];
	string NT_str 				= argv[9];
	string seed_str 			= argv[10];
	string positionFile 		= argv[11];

	stringstream NCELLSss(NCELLS_str);
	stringstream NVss(NV_str);
	stringstream tNss(tN_str);
	stringstream tCalA0ss(tumorCalA0_str);
	stringstream aCalA0ss(adiposeCalA0_str);
	stringstream gamss(gam_str);
	stringstream phi0ss(phi0_str);
	stringstream Drss(Dr_str);
	stringstream NTss(NT_str);
	stringstream seedss(seed_str);

	NCELLSss >> NCELLS;
	NVss >> NV;
	tNss >> tN;
	tCalA0ss >> tCalA0;
	aCalA0ss >> aCalA0;
	gamss >> gam;
	phi0ss >> phi0;
	Drss >> Dr;
	NTss >> NT_dbl;
	seedss >> seed;

	// cast input NT_dbl to integer
	NT = (int)NT_dbl;

	// number of tumor cells and adipose cells
	aN = NCELLS - tN;

	// open xyz file
	ofstream posout;
	posout.open(positionFile.c_str());
	if (!posout.is_open()){
		cout << "	** ERROR: position file " << positionFile << " could not be opened, ending." << endl;
		return 1;
	}

	// total number of vertices
	NVTOT = NCELLS*NV;

	// szList and nv (keep track of global vertex indices)
	vector<int> szList(NCELLS,0);
	vector<int> nv(NCELLS,0);
	nv.at(0) = NV;
	for (ci=1; ci<NCELLS; ci++){
		nv.at(ci) = NV;
		szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);
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
	tCalA0 *= NV*tan(PI/NV)/PI;
	aCalA0 *= NV*tan(PI/NV)/PI;

	// fundamental MD time units
	double dtMD, dt0, dt;

	dtMD 	= 1.0;
	dt0 	= timeStepMag*dtMD;
	dt 		= dt0;

	// output opening statement to console
	cout << "=======================================================" << endl << endl;

	cout << "		adipose2D.cpp 									" << endl;
	cout << "		Jack Treado, 2020   							" << endl;
	cout << "		Single active tumor cell 						" << endl;
	cout << "			in background of adipose tissue 			" << endl << endl;

	cout << "		NCELLS 		= " << NCELLS << "					" << endl;
	cout << "		# tcells 	= " << tN << "						" << endl;
	cout << "		# acells 	= " << aN << "						" << endl << endl;

	cout << "       NV (both) 	= " << NV << "						" << endl;
	cout << "		NVTOT 		= " << NVTOT << "					" << endl << endl;

	cout << "		tumor calA0 = " << tCalA0 << "					" << endl;
	cout << "		adi. calA0 	= " << aCalA0 << "					" << endl << endl;

	cout << "		phi0 		= " << phi0 << " 					" << endl;
	cout << "		v0 			= " << v0 << " 						" << endl;
	cout << "		ka 			= " << ka << "						" << endl;
	cout << "		kl 			= " << kl << "						" << endl;
	cout << "		gam 		= " << gam << " 					" << endl;
	cout << "		seed 		= " << seed << "					" << endl << endl;

	cout << "		pos file 	= " << positionFile << "			" << endl << endl;;
	
	cout << "=======================================================" << endl << endl;

	// seed random number generator
	srand48(seed);


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

	// initialize effective disk radius (for minimization), and l0 parameter
	areaSum = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		// set initial area
		if (ci < tN){
			lenscale = 1.0/sizeRatio;
			a0tmp = lenscale*lenscale;
			calA0tmp = tCalA0;
			nvtmp = NV;
		}
		else{
			lenscale = 1.0;
			a0tmp = 1.0;
			calA0tmp = aCalA0;
			nvtmp = NV;
		}

		// store preferred area
		a0.at(ci) 		= a0tmp;

		// set disk radius
		drad.at(ci) 	= 1.05*sqrt((2.0*a0tmp)/(nvtmp*sin(2.0*PI/nvtmp)));

		// set l0, vector radius
		l0.at(ci) 	= 2.0*lenscale*sqrt(PI*calA0tmp)/nvtmp;
		gi 			= szList.at(ci);
		for (vi=0; vi<nvtmp; vi++)
			vrad.at(gi+vi)	= 0.5*l0.at(ci)*del;

		// add to sum of particle areas (including contribution from vertices)
		areaSum 		+= a0tmp + 0.25*PI*pow(l0.at(ci)*del,2.0)*(0.5*nvtmp - 1);
		cout << "drad = " << drad.at(ci) << ", disk area = " << PI*pow(drad.at(ci),2.0) << ", l0 = " << l0.at(ci) << ", a0 = " << a0tmp << endl;
	}

	// determine box lengths from particle sizes and input packing fraction
	vector<double> L(NDIM,1.0);
	for (d=0; d<NDIM; d++)
		L.at(d) = sqrt(areaSum/phi0);

	// // initialize tumor cells in center of packing
	// for (ci=0; ci<tN; ci++){
	// 	dpos.at(NDIM*ci) = 0.1*L[0]*drand48() + 0.45*L[0];
	// 	dpos.at(NDIM*ci + 1) = 0.1*L[0]*drand48() + 0.45*L[0];
	// }

	// // initialize WAT cell centers randomly (but outside of center)
	// double xtmp, ytmp; 
	// int centerIt;
	// int centerItMax = 1e3;
	// bool inCenter = false;
	// for (ci=tN; ci<NCELLS; ci++){
	// 	centerIt = 0;
	// 	do {
	// 		// generate test points
	// 		xtmp = L[0]*drand48();
	// 		ytmp = L[1]*drand48();

	// 		// check to see if in center box
	// 		inCenter = xtmp < 0.55*L[0] && xtmp > 0.45*L[0];
	// 		inCenter = inCenter && ytmp < 0.55*L[1] && ytmp > 0.45*L[1];

	// 	} 
	// 	while (inCenter && centerIt < centerItMax);
	// 	if (centerIt == centerItMax){
	// 		cout << "\t ERROR: WAT cells could not be placed correctly. Ending.\n\n";
	// 		return 1;
	// 	}
	// 	else{
	// 		dpos.at(NDIM*ci) = xtmp;
	// 		dpos.at(NDIM*ci + 1) = ytmp;
	// 	}	
	// }

	// tumor at right
	for (ci=0; ci<tN; ci++){
		dpos.at(NDIM*ci) = (L[0] - 0.55*L[0])*drand48() + 0.55*L[0];
		dpos.at(NDIM*ci + 1) = L[1]*drand48();
	}

	// adipose at left
	for (ci=tN; ci<NCELLS; ci++){
		dpos.at(NDIM*ci) = 0.45*L[0]*drand48();
		dpos.at(NDIM*ci + 1) = L[1]*drand48();
	}



	/* * * * * * * * * * * * * * * * * * 

			INITIALIZE

				BOX-LINKED-LIST

	 * * * * * * * * * * * * * * * * * */


	// Cell-linked-list variables

	// box lengths in each direction
	vector<int> sb(NDIM,0);
	vector<double> lb(NDIM,0.0);
	int NBX = 1;
	for (d=0; d<NDIM; d++){
		// determine number of cells along given dimension by rmax
		sb[d] = round(L[d]/(2.5*l0.at(NCELLS-1)));

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

					MINIMZATION(S)

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


	// ----------------------------
	//
	// D P  M I N I M I Z A T I O N
	// 
	// ----------------------------



	// initialize disk velocity and force vectors
	vector<double> vvel(vertDOF,0.0);
	vector<double> vF(vertDOF,0.0);
	vector<double> vFold(vertDOF,0.0);

	// RESET FIRE VARIABLES
	P  			= 0;	
	fnorm 		= 0;
	vnorm 		= 0;
	alpha   	= alpha0;

	dtmax   	= 10*dt0;
	dtmin   	= 1e-4*dt0;
	dt 			= dt0;

	npPos      	= 0;
	npNeg      	= 0;
	npPMin      = 0;

	fireit    	= 0;
	fcheck  	= 10*Ftol;

	// linked list variables
	int boxid, bi, bj, pi, pj, sbtmp;
	int d0, dend;
	double U = 0.0;

	// shape force variables
	double Kl, l0tmp, atmp, li, lim1, cx, cy;
	double da, dli, dlim1;
	double lim2x, lim2y, lim1x, lim1y, lix, liy, lip1x, lip1y;
	double rim2x, rim2y, rim1x, rim1y, rix, riy, rip1x, rip1y, rip2x, rip2y;

	// loop until force relaxes
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

		// reset linked list 
		for (gi=0; gi<NVTOT+1; gi++)
			list[gi] = 0;

		// reset linked list head
		for (i=0; i<NBX; i++){
			head[i] = 0;
			last[i] = 0;
		}

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

		// FORCE UPDATE

		// interaction forces (USE BOX LINKED LIST)
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
								ftmp 				= eint*(1 - (rij/sij))/sij;
								fx 					= ftmp*(dx/rij);
								fy 					= ftmp*(dy/rij);

								// add to forces
								vF[NDIM*gi] 		-= fx;
								vF[NDIM*gi + 1] 	-= fy;

								vF[NDIM*gj] 		+= fx;
								vF[NDIM*gj + 1] 	+= fy;

								// increae potential energy
								U += 0.5*eint*pow((1 - (rij/sij)),2.0);
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
									ftmp 				= eint*(1 - (rij/sij))/sij;
									fx 					= ftmp*(dx/rij);
									fy 					= ftmp*(dy/rij);

									// add to forces
									vF[NDIM*gi] 		-= fx;
									vF[NDIM*gi + 1] 	-= fy;

									vF[NDIM*gj] 		+= fx;
									vF[NDIM*gj + 1] 	+= fy;

									// increae potential energy
									U += 0.5*eint*pow((1 - (rij/sij)),2.0);
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

					if (ci >= tN)
						l0tmp -= gam/(kl*nvtmp);

					// compute area deviation
					atmp = area(vpos,ci,L,nv,szList);
					da = atmp - a0tmp;

					// shape force parameters
					Kl = nvtmp*l0tmp*kl;

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

			// add to forces
			vF[NDIM*gi] 		+= ka*0.5*da*(rim1y - rip1y);
			vF[NDIM*gi + 1] 	+= ka*0.5*da*(rip1x - rim1x);


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
			vF[NDIM*gi] 		+= Kl*(dli*(lix/li) - dlim1*(lim1x/lim1));
			vF[NDIM*gi + 1] 	+= Kl*(dli*(liy/li) - dlim1*(lim1y/lim1));

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
		fcheck = fnorm/(NDIM*NCELLS);

		// update npPMin
		if (fcheck < Ftol)
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
			cout << "	** Uint = " << U << endl;
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
		cout << endl << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===========================================" << endl;
		cout << " 	F I R E 						" << endl;
		cout << "		M I N I M I Z A T I O N 	" << endl;
		cout << "	C O N V E R G E D! 				" << endl << endl;

		cout << "	(for initial DP minimization) " << endl;
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
		cout << "	** Uint = " << U << endl;
	}
























	/* * * * * * * * * * * * * * * * * * 

		COMPRESS TO TARGET

			PACKING FRACTION

	 * * * * * * * * * * * * * * * * * */


	// take equal sized steps in dphi
	int k, kmax, xind, yind;
	double dphi, pcheck, scaleFactor;

	// initial packing fraction
	phi0 = 0.0;
	for (ci=0; ci<NCELLS; ci++)
		phi0 += area(vpos,ci,L,nv,szList) + 0.25*PI*pow(l0[ci]*del,2.0)*(0.5*nv[ci] - 1.0);
	phi0 /= L[0]*L[1];

	// determine number of steps
	kmax = ((phiT - phi0)/dphi0);
	dphi = (phiT - phi0)/kmax;

	// compress to jamming, relax U and F using FIRE
	for (k=0; k<kmax+1; k++){

		// reduce Ftol
		if (phi0 > 0.9*phiT)
			Ftol = 1e-3*Ftol0;

		// grow or shrink particles by scale factor
		scaleFactor = sqrt((phi0 + dphi)/phi0);
		phi0 = 0.0;
		for (ci=0; ci<NCELLS; ci++){
			// scale preferred lengths
			l0[ci] *= scaleFactor;
			a0[ci] *= pow(scaleFactor,NDIM);

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
			phi0 += area(vpos,ci,L,nv,szList) + 0.25*PI*pow(l0[ci]*del,2.0)*(0.5*nv[ci] - 1.0);
		}
		phi0 /= L[0]*L[1];

		// RESET FIRE VARIABLES
		P  			= 0;	
		fnorm 		= 0;
		vnorm 		= 0;
		alpha   	= alpha0;

		dtmax   	= 10*dt0;
		dtmin   	= 1e-6*dt0;
		dt 			= dt0;

		npPos      	= 0;
		npNeg      	= 0;
		npPMin      = 0;

		fireit    	= 0;
		fcheck  	= 10*Ftol;

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

			// reset linked list 
			for (gi=0; gi<NVTOT+1; gi++)
				list[gi] = 0;

			// reset linked list head
			for (i=0; i<NBX; i++){
				head[i] = 0;
				last[i] = 0;
			}

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

			// FORCE UPDATE

			// interaction forces (USE BOX LINKED LIST)
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

			// normalize pressure by box area and number of particles
			pcheck /= NCELLS*L[0]*L[1];



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

						if (ci >= tN)
							l0tmp -= gam/(kl*nvtmp);

						// compute area deviation
						atmp = area(vpos,ci,L,nv,szList);
						da = atmp - a0tmp;

						// shape force parameters
						Kl = nvtmp*l0tmp*kl;

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

				// add to forces
				vF[NDIM*gi] 		+= ka*0.5*da*(rim1y - rip1y);
				vF[NDIM*gi + 1] 	+= ka*0.5*da*(rip1x - rim1x);


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
				vF[NDIM*gi] 		+= Kl*(dli*(lix/li) - dlim1*(lim1x/lim1));
				vF[NDIM*gi + 1] 	+= Kl*(dli*(liy/li) - dlim1*(lim1y/lim1));

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
			fcheck = fnorm/(NDIM*NCELLS);

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
				cout << "	** pcheck 	= " << pcheck << endl << endl;

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
			cout << endl << endl << endl;
			cout << "===========================================" << endl;
			cout << " 	F I R E 						" << endl;
			cout << "		M I N I M I Z A T I O N 	" << endl;
			cout << "	C O N V E R G E D! 				" << endl;
			cout << "===========================================" << endl;
			cout << endl;
			cout << "	** fireit 	= " << fireit << endl;
			cout << "	** fcheck 	= " << fcheck << endl;
			cout << "	** pcheck 	= " << pcheck << endl << endl;

			cout << "	** vnorm 	= " << vnorm << endl;
			cout << "	** dt 		= " << dt << endl;
			cout << "	** P 		= " << P << endl;
			cout << "	** Pdir 	= " << P/(fnorm*vnorm) << endl;
			cout << "	** alpha 	= " << alpha << endl << endl;

			cout << "	** current phi = " << phi0 << endl;
			cout << "	** k = " << k << "/" << kmax << endl << endl;

			// print vertex positions
			cout << "\t** PRINTING POSITIONS TO FILE... " << endl << endl << endl;
			printPos(posout, vpos, a0, l0, L, nv, szList, phi0, NCELLS);
		}
	}

	// print vertex positions
	cout << "\t** PRINTING POSITIONS TO FILE... " << endl << endl << endl;
	printPos(posout, vpos, a0, l0, L, nv, szList, phi0, NCELLS);































	/* * * * * * * * * * * * * * * * * * 

		OVERDAMPED DYNAMICS

			WITH ACTIVE TUMOR CELL

	 * * * * * * * * * * * * * * * * * */

	// DYNAMICS VARIABLES
	int tt, tcells; 
	double r1, r2, grv, ux, uy, rnorm, psitmp, dpsi, v0tmp;

	// initial directors
	vector<double> psi(tN,0.0);
	for (ci=0; ci<tN; ci++)
		psi.at(ci) = PI;

	// reset for active dynamocs
	dt = dt0;

	// draw draw 
	cout << endl;
	cout << "--  BEGINNING ACTIVE TUMOR CELL DYNAMICS" << endl;

	cout << "\t** Looping over time for NT = " << NT << " time steps with dt = " << dt << endl;
	for(tt=0; tt<NT; tt++){

		// pbcs and reset forces
		for (i=0; i<vertDOF; i++){
			// recenter in box
			if (vpos[i] > L[i % NDIM])
				vpos[i] -= L[i % NDIM];
			else if (vpos[i] < 0)
				vpos[i] += L[i % NDIM];

			// reset forces
			vF[i] = 0.0;
		}

		// reset linked list 
		for (gi=0; gi<(NVTOT+1); gi++)
			list[gi] = 0;

		// reset linked list head
		for (i=0; i<NBX; i++){
			head[i] = 0;
			last[i] = 0;
		}

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

		// FORCE UPDATE

		// interaction forces (USE BOX LINKED LIST)
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

		// normalize pressure by box area and number of particles
		pcheck /= NCELLS*L[0]*L[1];



		// shape forces (loop over global vertex labels)
		ci = 0;
		tcells = 1;
		for (gi=0; gi<NVTOT; gi++){

			// -- Area force (and get cell index ci)
			if (ci < NCELLS){
				if (gi == szList[ci]){
					// compute shape parameter
					nvtmp = nv[ci];
					a0tmp = a0[ci];
					l0tmp = l0[ci];

					if (ci >= tN){
						tcells = 0;
						l0tmp -= gam/(kl*nvtmp);
					}
					else
						tcells = 1;

					// compute area deviation
					atmp = area(vpos,ci,L,nv,szList);
					da = atmp - a0tmp;

					// shape force parameters
					Kl = nvtmp*l0tmp*kl;

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

			// add to forces
			vF[NDIM*gi] 		+= ka*0.5*da*(rim1y - rip1y);
			vF[NDIM*gi + 1] 	+= ka*0.5*da*(rip1x - rim1x);


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
			vF[NDIM*gi] 		+= Kl*(dli*(lix/li) - dlim1*(lim1x/lim1));
			vF[NDIM*gi + 1] 	+= Kl*(dli*(liy/li) - dlim1*(lim1y/lim1));

			// active force (only if tumor cell vertices)
			if (tcells == 1){
				// get angular distance from psi
				psitmp = atan2(riy,rix);
				dpsi = psitmp - psi[ci-1];
				dpsi -= 2.0*PI*round(dpsi/(2.0*PI));

				// get velocity scale
				v0tmp = vmin + (v0-vmin)*exp(-pow(dpsi,2.0)/(2.0*Ds));

				// get unit vectors
				rnorm = sqrt(rix*rix + riy*riy);
				ux = rix/rnorm;
				uy = riy/rnorm;

				// add to forces
				vF[NDIM*gi] += v0tmp*ux;
				vF[NDIM*gi + 1] += v0tmp*uy;
			}

			// update old coordinates
			rim2x = rim1x;
			rim1x = rix;
			rix = rip1x;

			rim2y = rim1y;
			rim1y = riy;
			riy = rip1y;
		}

		// update psi directions using box muller transform
		for (ci=0; ci<tN; ci++){
			r1 = drand48();
			r2 = drand48();
			grv = sqrt(-2.0*log(r1))*cos(2.0*PI*r2);

			psi[ci] += dt*2.0*Dr*grv;
		}
		


		// BROWNIAN MOTION VIA EULER UPDATE
		for (i=0; i<vertDOF; i++)
			vpos[i] += dt*vF[i];



		// print message console, print energy to file
		if (tt % NACTIVESKIP == 0){
			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << "			active cell 					" << endl;
			cout << "===========================================" << endl;
			cout << endl;
			cout << "	** tt = " << tt << endl;
			cout << "	** p = " << pcheck << endl;
			cout << "	** psi(0) = " << psi[0] << endl;

			// print vertex positions to check placement
			cout << "\t** PRINTING POSITIONS TO FILE... " << endl;
			printPos(posout, vpos, a0, l0, L, nv, szList, phi0, NCELLS);
		}
	}


	// close open objects
	posout.close();


	// print to console, return
	cout << "	** FINISHED MAIN FOR adipose2D.cpp, ENDING." << endl << endl << endl;
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
void printPos(ofstream& posout, vector<double>& vpos, vector<double>& a0, vector<double>& l0, vector<double>& L, vector<int>& nv, vector<int>& szList, double phi, int NCELLS){
	// local variables
	int ci, cj, vi, gi, ctmp;
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

		// cell information
		posout << setw(w) << left << "CINFO";
		posout << setw(w) << right << nv.at(ci);
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






















