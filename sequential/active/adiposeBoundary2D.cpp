/* 

	MAIN FILE FOR TUMOR INVASION INTO
		A MODEL ADIPOSE TISSUE

	** 	fixed number of adipose and tumor cells, 
		initialized at interface with partial periodic
		boundaries

	Jack Treado
	01/20/2021, in the time of covid (and joe biden!)

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
const double PI 			= 4.0*atan(1);
const int w 				= 10;
const int wnum 				= 25;
const int pnum 				= 8;

// simulation constants
const double timeStepMag 	= 0.0025;
const double polyd 			= 0.1;
const int nvmin 			= 12;
const double phi0 			= 0.6;
const double phiT 			= 0.9;
const double dphi 			= 0.005;

// FIRE constants for initial minimizations (SP + DP)
const double alpha0      	= 0.25;
const double finc        	= 1.1;
const double fdec        	= 0.5;
const double falpha      	= 0.99;
const double Ftol 			= 1e-8;
const double Ptol 			= 1e-8;
const double minDelA 		= 1e-2;

const int NSKIP 			= 2e4;
const int NMIN        		= 10;
const int NNEGMAX     		= 1000;
const int NDELAY      		= 50;
const int itmax       		= 5e7;

// DP force constants
const double ka 			= 1.0;			// area spring (should be = 1)
const double kl 			= 0.5;			// contractility
const double kb 			= 0;			// bending modulus
const double eint 			= 0.5;			// interaction energy scale
const double del 			= 1.0;			// radius of vertices in units of l0
const double aCalA0 		= 0.9;			// adipocyte deformability

// tumor invasion variables
const double Ds 			= 0.2;			// spread of velocity coupling along tumor cell boundary
const double Dr 			= 0.2;			// bare angular diffusion
const double Drmin 			= 1e-3;			// min angular diffusion, mimics aligning to collagen
const double pinbreak 		= 1.0; 			// fraction of rho0 that breaks a WAT pin spring
const double kpin 			= 0.1;			// pinning spring stiffness

// FUNCTION PROTOTYPES
int gindex(int ci, int vi, vector<int>& szList);
void cindices(int& ci, int& vi, int gi, int NCELLS, vector<int>& szList);

double area(vector<double>& dpos, int ci, vector<double>& L, vector<int>& nv, vector<int>& szList);
double perimeter(vector<double>& dpos, int ci, vector<double>& L, vector<int>& nv, vector<int>& szList);

// print to file
void printPos(ofstream& posout, vector<double>& vpos, vector<double>& vvel, vector<double>& a0, vector<double>& l0, vector<double>& L, vector<int>& nv, vector<int>& szList, double phi, int NCELLS, int tN); 


// MAIN
int main(int argc, char const *argv[]){
	// variables for indexing loops
	int i, ci, cj, vi, gi, gj, d;

	// parameters to be read in 
	int aN, tN, NCELLS, NT, NACTIVESKIP, NV, NVTOT, cellDOF, vertDOF, seed;
	double NT_dbl, NACTIVESKIP_dbl, areaRatio, tCalA0, l1, l2, kint, v0, vmin, Dr, dDr, phi;

	// read in parameters from command line input
	// test g++ -O3 sequential/active/adiposeBoundary2D.cpp -o tumor.o
	// test: ./tumor.o 8 2 24 1.10 0.01 0.011 0.05 0.1 1e5 2e3 1 pos.tes
	string aN_str 				= argv[1];
	string areaRatio_str 		= argv[2];
	string NV_str 				= argv[3];
	string tumorCalA0_str 		= argv[4];
	string l1_str 				= argv[5];
	string l2_str 				= argv[6];
	string v0_str 				= argv[7];
	string dDr_str 				= argv[8];
	string NT_str 				= argv[9];
	string NACTIVESKIP_str 		= argv[10];
	string seed_str 			= argv[11];
	string positionFile 		= argv[12];

	stringstream aNss(aN_str);
	stringstream areaRatioss(areaRatio_str);
	stringstream NVss(NV_str);
	stringstream tCalA0ss(tumorCalA0_str);
	stringstream l1ss(l1_str);
	stringstream l2ss(l2_str);
	stringstream v0ss(v0_str);
	stringstream dDrss(dDr_str);
	stringstream NTss(NT_str);
	stringstream NAss(NACTIVESKIP_str);
	stringstream seedss(seed_str);

	aNss >> aN;
	areaRatioss >> areaRatio;
	NVss >> NV;
	tCalA0ss >> tCalA0;
	l1ss >> l1;
	l2ss >> l2;
	v0ss >> v0;
	dDrss >> dDr;
	NTss >> NT_dbl;
	NAss >> NACTIVESKIP_dbl;
	seedss >> seed;


	// interaction parameter based on l1 and l2
	if (l1 > l2){
		cout << "	** ERROR: l1 = " << l1 << " which is > l2 = " << l2 << ". Ending." << endl;
		return 1;
	}
	kint = (eint*l1)/(l2 - l1);


	// seed random number generator
	srand48(seed);

	// cast input time scales to integer
	NT 				= (int)NT_dbl;
	NACTIVESKIP 	= (int)NACTIVESKIP_dbl;

	// minimum velocity
	vmin = 1e-2*v0;

	// open xyz file
	ofstream posout;
	posout.open(positionFile.c_str());
	if (!posout.is_open()){
		cout << "	** ERROR: position file " << positionFile << " could not be opened, ending." << endl;
		return 1;
	}

	// determine number of tumor cells by area ratio
	tN = round(aN * areaRatio);
	if (tN < 1){
		cout << "	** ERROR: tN = " << tN << ", which is less than 1. Ending." << endl;
		return 1;
	}

	// total number of cells
	NCELLS = tN + aN;

	// total number of vertices
	
	// szList and nv (keep track of global vertex indices)
	int nvtmp, NVtumor, NVadipose;
	double r1, r2, grv;
	vector<int> szList(NCELLS,0);
	vector<int> nv(NCELLS,0);
	nv.at(0) = NV;
	NVtumor = NV;
	NVadipose = 0;
	NVTOT = NV;

	// draw random number of vertices for each particle, add to NVTOT
	int imin, imax, rmin, rmax;
	rmin = NV;
	imin = 0;
	rmax = NV;
	imax = 0;
	for (ci=1; ci<NCELLS; ci++){
		// use Box-Muller to generate polydisperse sample
		r1 = drand48();
		r2 = drand48();
		grv = sqrt(-2.0*log(r1))*cos(2.0*PI*r2);
		nvtmp = floor(polyd*NV*grv + NV);
		if (nvtmp < nvmin)
			nvtmp = nvmin;

		if (nvtmp < rmin){
			rmin = nvtmp;
			imin = ci;
		}
		else if(nvtmp > rmax){
			rmax = nvtmp;
			imax = ci;
		}

		// store size of cell ci
		nv.at(ci) = nvtmp;
		szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);

		// add to total NV count
		NVTOT += nvtmp;

		// add to cell-specific count
		if (ci < tN)
			NVtumor += nvtmp;
	}

	// number of adipose vertices
	NVadipose = NVTOT - NVtumor;

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

	cout << "		adiposeBoundary2D.cpp 							" << endl;
	cout << "		Jack Treado, 2020   							" << endl;

	cout << "		NCELLS 		= " << NCELLS << "					" << endl;
	cout << "		# tcells 	= " << tN << "						" << endl;
	cout << "		# acells 	= " << aN << "						" << endl << endl;

	cout << "		NCELLS 		= " << NCELLS << "					" << endl;
	cout << "		NV (mean)	= " << NV << "						" << endl;
	cout << "		NVTOT 		= " << NVTOT << "					" << endl;
	cout << "		NVtumor 	= " << NVtumor << "					" << endl;
	cout << "		NVadipose 	= " << NVadipose << "					" << endl;
	cout << "		max NV 		= " << rmax << " at ci 	= " << imax << endl;
	cout << "		min NV 		= " << rmin << " at ci 	= " << imin << endl;
	cout << endl;

	cout << "		tumor calA0 = " << tCalA0 << "					" << endl;
	cout << "		adi. calA0 	= " << aCalA0 << "					" << endl << endl;

	cout << "		phi0 		= " << phi0 << " 					" << endl;
	cout << "		ka 			= " << ka << "						" << endl;
	cout << "		kl 			= " << kl << "						" << endl;
	cout << "		v0 			= " << v0 << " 						" << endl;
	cout << "		Dr 			= " << Dr << " 						" << endl;
	cout << "		seed 		= " << seed << "					" << endl << endl;

	cout << "		pos file 	= " << positionFile << "			" << endl << endl;;
	
	cout << "=======================================================" << endl << endl;



	/* * * * * * * * * * * * * * * * * * 

			INITIALIZE

				PARTICLES

	 * * * * * * * * * * * * * * * * * */

	// initialization variables
	double a0tmp, l0max, lenscale, calA0tmp, areaSum = 0.0;

	// initialize vectors for storing coordinates, shape information
	vector<double> vrad(NVTOT,1.0);
	vector<double> drad(NCELLS,1.0);

	vector<double> vpos(vertDOF,0.0);
	vector<double> dpos(cellDOF,0.0);

	vector<double> a0(NCELLS,1.0);
	vector<double> aMin(NCELLS,1.0);
	vector<double> l0(NCELLS,1.0);

	// initialize effective disk radius (for minimization), and l0 parameter
	areaSum = 0.0;
	l0max = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		// number of vertices
		nvtmp 			= nv.at(ci);

		// set initial areas
		lenscale = (double)nvtmp/NV;
		if (ci < tN)
			calA0tmp = tCalA0*(nvtmp*tan(PI/nvtmp)/PI);
		else{
			lenscale *= sqrt(areaRatio);
			calA0tmp = aCalA0*(nvtmp*tan(PI/nvtmp)/PI);
		}

		// store preferred area
		a0tmp 			= lenscale*lenscale;
		a0.at(ci) 		= a0tmp;

		// set disk radius
		drad.at(ci) 	= 1.05*sqrt((2.0*a0tmp)/(nvtmp*sin(2.0*PI/nvtmp)));

		// set l0, vector radius
		l0.at(ci) 		= 2.0*lenscale*sqrt(PI*calA0tmp)/nvtmp;
		gi 				= szList.at(ci);
		for (vi=0; vi<nvtmp; vi++)
			vrad.at(gi+vi)	= 0.5*l0.at(ci)*del;

		// check if l0max for linked list
		if (l0.at(ci) > l0max)
			l0max = l0.at(ci);

		// add to sum of tumor cell areas (including contribution from vertices)
		if (ci < tN)
			areaSum += a0tmp + 0.25*PI*pow(l0.at(ci)*del,2.0)*(0.5*nvtmp - 1);

		// output info
		cout << "\t ** ci = " << ci << ";  nv = " << nvtmp << ";  drad = " << drad.at(ci) << ";  disk area = " << PI*pow(drad.at(ci),2.0) << ", l0 = " << l0.at(ci) << ", a0 = " << a0tmp << endl;
	}

	// determine box lengths from particle sizes and input packing fraction
	vector<double> L(NDIM,0.0);
	L.at(1) = sqrt(areaSum/phi0);
	L.at(0) = 2.0*L[1];

	// initial packing fraction
	phi = areaSum/(L[1]*L[1]);

	// initialize tumor cell centers in left side of the box
	for (ci=0; ci<tN; ci++){
		dpos.at(NDIM*ci) 		= (L[1] - 2.0*drad[ci])*drand48() + drad[ci];
		dpos.at(NDIM*ci + 1) 	= L[1]*drand48();
	}

	// initialize WAT cell centers to the right
	for (ci=tN; ci<NCELLS; ci++){
		dpos.at(NDIM*ci) 		= (L[0] - L[1] - 2.0*drad[ci])*drand48() + L[1] + drad[ci];
		dpos.at(NDIM*ci + 1) 	= L[1]*drand48();
	}


	/* * * * * * * * * * * * * * * * * * 

			INITIAL SP

				FIRE MINIMZATION

	 * * * * * * * * * * * * * * * * * */

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

				// true distance (no PBCs in X)
				dx = xj - xi;
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

			// x boundary forces
			if (ci < tN) {
				if (xi < drad[ci])
					dF[NDIM*ci] += eint*(1.0 - (xi/drad[ci]))/drad[ci];
				else if (xi > L[1] - drad[ci])
					dF[NDIM*ci] -= eint*(1.0 - ((L[1] - xi)/drad[ci]))/drad[ci];
			}
			else {
				if (xi < L[1] + drad[ci])
					dF[NDIM*ci] += eint*(1.0 - ((xi - L[1])/drad[ci]))/drad[ci];
				else if (xi > L[0] - drad[ci])
					dF[NDIM*ci] -= eint*(1.0 - ((L[0] - xi)/drad[ci]))/drad[ci];
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
			vpos.at(NDIM*gi) 		= lenscale*cos((2.0*PI*vi)/nv.at(ci)) + dpos.at(NDIM*ci) + 1e-2*l0.at(ci)*drand48();
			vpos.at(NDIM*gi + 1)	= lenscale*sin((2.0*PI*vi)/nv.at(ci)) + dpos.at(NDIM*ci + 1) + 1e-2*l0.at(ci)*drand48();
		}
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
		sb[d] = round(L[d]/(1.5*l0max));

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

		// right-hand bc (no periodicity in x direction)
		if ((i+1) % scx == 0){
			nn[i][0] = -1;
			nn[i][2] = -1;
			nn[i][3] = -1;
		}
	} 

	// linked-list variables
	vector<int> head(NBX,0);
	vector<int> last(NBX,0);
	vector<int> list(NVTOT+1,0);














	/* * * * * * * * * * * * * * * * * * 

		COMPRESS TO TARGET

			PACKING FRACTION

				IN SEPARATE CHAMBERS

	 * * * * * * * * * * * * * * * * * */


	// initialize disk velocity and force vectors
	vector<double> vvel(vertDOF,0.0);
	vector<double> vF(vertDOF,0.0);
	vector<double> vFold(vertDOF,0.0);

	// linked list variables
	int boxid, bi, bj, pi, pj, sbtmp;

	// shape force variables
	double rho0, fa, fl, fb, l0tmp, atmp, ri, li, lim1, cx, cy;
	double da, dli, dlim1;
	double lim2x, lim2y, lim1x, lim1y, lix, liy, lip1x, lip1y;
	double rim2x, rim2y, rim1x, rim1y, rix, riy, rip1x, rip1y, rip2x, rip2y;

	// packing fractions for different cell types
	double tumorPhi, adiposePhi, phiMin;
	tumorPhi = 0.0;
	adiposePhi = 0.0;
	for (ci=0; ci<tN; ci++)
		tumorPhi += area(vpos,ci,L,nv,szList) + 0.25*PI*pow(l0[ci]*del,2.0)*(0.5*nv[ci] - 1.0);
	for (ci=tN+1; ci<NCELLS; ci++)
		adiposePhi += area(vpos,ci,L,nv,szList) + 0.25*PI*pow(l0[ci]*del,2.0)*(0.5*nv[ci] - 1.0);
	tumorPhi /= L[1]*L[1];
	adiposePhi /= L[1]*L[1];

	if (tumorPhi < adiposePhi)
		phiMin = tumorPhi;
	else
		phiMin = adiposePhi;

	// // area stress parameters
	// double minTumorDelA, minAdiposeDelA, minTotalDelA;
	// minTumorDelA = 0.0;
	// minAdiposeDelA = 0.0;
	// minTotalDelA = 0.0;

	// size scaling variables
	int k, kmax, xind, yind;
	double pcheck, scaleFactor;

	// determine number of steps
	k = 0;
	kmax = 1e4;

	// compress to jamming, relax U and F using FIRE
	cout << endl << endl << endl;
	cout << "	** Compressing to target packing fraction phiT = " << phiT << endl;
	// while (minTotalDelA < minDelA && k < kmax){
	while (phiMin < phiT && k < kmax){

		// increment counter
		k++;

		// RESET FIRE VARIABLES
		P  			= 0;	
		fnorm 		= 0;
		vnorm 		= 0;
		alpha   	= alpha0;

		dtmax   	= 10.0*dt0;
		dtmin   	= 1e-6*dt0;
		dt 			= dt0;

		npPos      	= 0;
		npNeg      	= 0;
		npPMin      = 0;

		fireit    	= 0;
		fcheck  	= 10.0*Ftol;

		// set length constant (variable due to particle growth)
		rho0 = sqrt(a0.at(0));

		// reset forces and velocities
		fill(vF.begin(), vF.end(), 0.0);
		fill(vFold.begin(), vFold.end(), 0.0);
		fill(vvel.begin(), vvel.end(), 0.0);

		// RELAX FORCES USING FIRE
		while ((fcheck > Ftol || npPMin < NMIN) && fireit < itmax){
			// VV POSITION UPDATE
			for (i=0; i<vertDOF; i++){
				// update position
				vpos[i] += dt*vvel[i] + 0.5*dt*dt*vF[i];

				// recenter in box (only if y)
				if (i % NDIM == 1){
					if (vpos[i] > L[1])
						vpos[i] -= L[1];
					else if (vpos[i] < 0)
						vpos[i] += L[1];
				}

				// reset forces
				vF[i] = 0.0;
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

						// particle distance (box boundaries in x direction)
						dx = vpos[NDIM*gj] - vpos[NDIM*gi];
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
						// if on boundary, don't check
						if (nn[bi][bj] == -1)
							continue;

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

							// particle distance (box boundaries in x direction)
							dx = vpos[NDIM*gj] - vpos[NDIM*gi];
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

						/// compute area deviation
						atmp = area(vpos,ci,L,nv,szList);
						da = (atmp/a0tmp) - 1.0;

						// shape force parameters (kl and kl are unitless energy ratios)
						fa = ka*da*(rho0/a0tmp);
						fl = kl*(rho0/l0tmp);
						fb = kb*(rho0/(l0tmp*l0tmp));

						// compute cell center of mass
						xi = vpos[NDIM*gi];
						yi = vpos[NDIM*gi + 1];
						cx = xi; 
						cy = yi;
						for (vi=1; vi<nvtmp; vi++){
							dx = vpos[NDIM*(gi+vi)] - xi;

							dy = vpos[NDIM*(gi+vi) + 1] - yi;
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

						rim1y = vpos[NDIM*im1[gi] + 1] - cy;
						rim1y -= L[1]*round(rim1y/L[1]);

						rim2x = vpos[NDIM*im1[im1[gi]]] - cx;

						rim2y = vpos[NDIM*im1[im1[gi]] + 1] - cy;
						rim2y -= L[1]*round(rim2y/L[1]);

						// boundary forces
						for (vi=0; vi<nvtmp; vi++){
							// positions using global indexing
							xi = vpos[NDIM*(gi+vi)];
							ri = vrad[gi+vi];

							// if near a wall, add to force

							// tumors
							if (ci < tN){
								if (xi < ri)
									vF[NDIM*(gi+vi)] += eint*(1.0 - (xi/ri))/ri;
								else if (xi > L[1] - ri)
									vF[NDIM*(gi+vi)] -= eint*(1.0 - ((L[1] - xi)/ri))/ri;
							}
							// adipocytes
							else{
								if (xi < L[1] + ri)
									vF[NDIM*(gi+vi)] += eint*(1.0 - ((xi - L[1])/ri))/ri;
								else if (xi > L[0] - ri)
									vF[NDIM*(gi+vi)] -= eint*(1.0 - ((L[0] - xi)/ri))/ri;
							}
						}

						// increment cell index
						ci++;
					}
				}


				// get next adjacent vertices
				rip1x = vpos[NDIM*ip1[gi]] - cx;

				rip1y = vpos[NDIM*ip1[gi] + 1] - cy;
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
				if ((ci-1) >= tN){
					// segment vectors for ip2
					rip2x = vpos[NDIM*ip1[ip1[gi]]] - cx;

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
			fcheck = fnorm/sqrt(NDIM*NVTOT);

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
		cout << "	* k 				= " << k << endl;
		cout << "	* dr 				= " << scaleFactor << endl;
		cout << "	* tumorPhi 			= " << tumorPhi << endl;
		cout << "	* adiposePhi 		= " << adiposePhi << endl;
		cout << "	* fcheck 			= " << fcheck << endl;
		cout << "	* pcheck 			= " << pcheck << endl;
		// cout << "	* minTotalDelA 		= " << minTotalDelA << endl;
		cout << endl;

		// // if first min, get aMin for each cell
		// if (k == 1){
		// 	for (ci=0; ci<NCELLS; ci++)
		// 		aMin[ci] = area(vpos,ci,L,nv,szList);
		// }

		
		// grow or shrink particles by scale factor
		tumorPhi = 0.0;
		adiposePhi = 0.0;
		for (ci=0; ci<NCELLS; ci++){
			if (ci < tN)
				tumorPhi += a0[ci] + 0.25*PI*pow(l0[ci]*del,2.0)*(0.5*nv[ci] - 1);
			else
				adiposePhi += area(vpos,ci,L,nv,szList) + 0.25*PI*pow(l0[ci]*del,2.0)*(0.5*nv[ci] - 1);
		}
		tumorPhi /= L[1]*L[1];
		adiposePhi /= L[1]*L[1];

		if (tumorPhi < adiposePhi)
			phiMin = tumorPhi;
		else
			phiMin = adiposePhi;

		// // get minimum fractional area strains
		// minTumorDelA = abs((area(vpos,0,L,nv,szList)/aMin[0]) - 1.0);
		// minAdiposeDelA = abs((area(vpos,tN,L,nv,szList)/aMin[tN]) - 1.0);
		// for (ci=0; ci<NCELLS; ci++){
		// 	da = abs((area(vpos,ci,L,nv,szList)/aMin[ci]) - 1.0);
		// 	if (ci < tN){
		// 		if (da < minTumorDelA)
		// 			minTumorDelA = da;
		// 	}
		// 	else{
		// 		if (da < minAdiposeDelA)
		// 			minAdiposeDelA = da;
		// 	}
		// }

		// // check total min delA
		// if (minTumorDelA < minAdiposeDelA)
		// 	minTotalDelA = minTumorDelA;
		// else
		// 	minTotalDelA = minAdiposeDelA;

		// // break if satisfies
		// if (minTotalDelA > minDelA){
		// 	for (ci=0; ci<NCELLS; ci++){
		// 		da = abs((area(vpos,ci,L,nv,szList)/aMin[ci]) - 1.0);
		// 		cout << "	** ci = " << ci << ", da = " << da << endl;
		// 	}
		// 	cout << "	** minTotalDelA = " << minTotalDelA << ", less than minDelA = " << minDelA << ", calling confluent and ending. " << endl;
		// 	break;
		// }

		// grow particles by scale factor
		for (ci=0; ci<NCELLS; ci++){
			// check whether or not to scale particles
			if (ci < tN){
				if (tumorPhi > phiT)
					continue;
				else
					scaleFactor = sqrt((tumorPhi + dphi)/tumorPhi);
			}
			else {
				if (adiposePhi > phiT)
					continue;
				else
					scaleFactor = sqrt((adiposePhi + dphi)/adiposePhi);
			}

			// scale preferred lengths
			l0[ci] *= scaleFactor;
			a0[ci] *= scaleFactor*scaleFactor;
			// aMin[ci] *= scaleFactor*scaleFactor;

			// first global index for ci
			gi = szList.at(ci);

			// compute cell center of mass
			xi = vpos[NDIM*gi];
			yi = vpos[NDIM*gi + 1];
			cx = xi; 
			cy = yi;
			for (vi=1; vi<nv.at(ci); vi++){
				dx = vpos.at(NDIM*(gi+vi)) - xi;

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
				dx = vpos.at(xind) - cx;

				dy = vpos.at(yind) - cy;
				dy -= L[1]*round(dy/L[1]);

				// update vertex positions
				vpos.at(xind)		+= (scaleFactor - 1.0)*dx;
				vpos.at(yind) 		+= (scaleFactor - 1.0)*dy;

				// scale vertex radii
				vrad.at(gi+vi) 		*= scaleFactor;
			}
		}
	}
	if (k == kmax){
		cout << "	** ERROR: In adipose/tumor compression step, k reached kmax without getting to phiT. Ending." << endl;
		return 1;
	}




























	/* * * * * * * * * * * * * * * * * * 

			RESET AND RESIZE

				BOX-LINKED-LIST

	 * * * * * * * * * * * * * * * * * */


	// re-determine l0max
	l0max = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		if (l0[ci] > l0max)
			l0max = l0[ci];
	}

	// recount cell box dimensions, total number of cells needed
	NBX = 1;
	for (d=0; d<NDIM; d++){
		// determine number of cells along given dimension by rmax
		sb[d] = round(L[d]/(2.0*l0max));

		// just in case, if < 3, change to 3 so box neighbor checking will work
		if (sb[d] < 3)
			sb[d] = 3;

		// determine box length by number of cells
		lb[d] = L[d]/sb[d];

		// recount total number of cells
		NBX *= sb[d];
	}

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

		// right-hand bc (no periodicity in x direction)
		if ((i+1) % scx == 0){
			nn[i][0] = -1;
			nn[i][2] = -1;
			nn[i][3] = -1;
		}
	} 

	// resize linked-list variables
	head.resize(NBX);
	last.resize(NBX);





































	/* * * * * * * * * * * * * * * * * * 

		OVERDAMPED DYNAMICS

			WITH ACTIVE TUMOR CELL

	 * * * * * * * * * * * * * * * * * */

	// print vertex positions to check placement
	cout << "\t** PRINTING INITIAL POSITIONS TO FILE... " << endl;
	printPos(posout, vpos, vF, a0, l0, L, nv, szList, phi, NCELLS, tN);

	// DYNAMICS VARIABLES
	int tt, tcells; 
	double ux, uy, rnorm, psitmp, dpsi, v0tmp, psiMean, psiStd;

	// attractive shell for tumor vertices
	double shellij, cutij;

	// WAT pin locations
	int pinCellInd = 0;
	double dpin, dpinx, dpiny = 0.0;
	vector<double> watpins(NDIM*aN,0.0);
	vector<int> pinattached(aN,1);
	for (ci=tN; ci<NCELLS; ci++){
		// first global index for ci
		gi = szList.at(ci);

		// compute cell center of mass
		xi = vpos[NDIM*gi];
		yi = vpos[NDIM*gi + 1];
		cx = xi; 
		cy = yi;
		for (vi=1; vi<nv.at(ci); vi++){
			dx = vpos.at(NDIM*(gi+vi)) - xi;

			dy = vpos.at(NDIM*(gi+vi) + 1) - yi;
			dy -= L[1]*round(dy/L[1]);

			xi += dx;
			yi += dy;

			cx += xi;
			cy += yi;
		}
		cx /= nv.at(ci);
		cy /= nv.at(ci);

		// save pin locations
		watpins.at(NDIM*(ci-tN)) = cx;
		watpins.at(NDIM*(ci-tN) + 1) = cy;
	}

	// initialize vertex-vertex contact network
	int NVVCTCS = 0.5*NVTOT*(NVTOT-1);
	double zta = 0.0;
	vector<bool> gij(NVVCTCS,0);

	// initialize directors
	double Drtmp = Dr;
	vector<double> DrList(tN,Dr);
	vector<double> psi(tN,0.0);

	// reset for active dynamocs
	dt = dt0;

	// draw draw 
	cout << endl;
	cout << "--  BEGINNING ACTIVE TUMOR CELL DYNAMICS" << endl;
	cout << "--  Looping over time for NT = " << NT << " time steps with dt = " << dt << endl;
	for(tt=0; tt<NT; tt++){

		// pbcs and reset forces
		for (i=0; i<vertDOF; i++){
			// recenter in box (only if y)
			if (i % NDIM == 1){
				if (vpos[i] > L[1])
					vpos[i] -= L[1];
				else if (vpos[i] < 0)
					vpos[i] += L[1];
			}

			// reset forces
			vF[i] = 0.0;
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

					// attraction distances
					shellij = (1.0 + l2)*sij;
					cutij = (1.0 + l1)*sij;

					// particle distance
					dx = vpos[NDIM*gj] - vpos[NDIM*gi];
					if (dx < shellij){
						dy = vpos[NDIM*gj + 1] - vpos[NDIM*gi + 1];
						dy -= L[1]*round(dy/L[1]);
						if (dy < shellij){
							rij = sqrt(dx*dx + dy*dy);
							// if two tumor cells, check for attraction
							if (rij < shellij && rij > cutij && gi < NVtumor && gj < NVtumor){
								// force scale
								ftmp 				= kint*((rij/sij) - 1.0 - l2)/sij;
								fx 					= ftmp*(dx/rij);
								fy 					= ftmp*(dy/rij);

								// add to forces
								vF[NDIM*gi] 		-= fx;
								vF[NDIM*gi + 1] 	-= fy;

								vF[NDIM*gj] 		+= fx;
								vF[NDIM*gj + 1] 	+= fy;

								// add to virial expression for pressure
								pcheck += dx*fx + dy*fy;

								// add to contact network
								if (gj > gi)
									gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 1;
								else
									gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2] = 1;
							}
							else if (rij < cutij && rij > sij && gi < NVtumor && gj < NVtumor){
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

								// add to contact network
								if (gj > gi)
									gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 1;
								else
									gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2] = 1;
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

								// add to contact network
								if (gj > gi)
									gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 1;
								else
									gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2] = 1;
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
						shellij = (1.0 + l2)*sij;
						cutij = (1.0 + l1)*sij;

						// particle distance
						dx = vpos[NDIM*gj] - vpos[NDIM*gi];
						if (dx < shellij){
							dy = vpos[NDIM*gj + 1] - vpos[NDIM*gi + 1];
							dy -= L[1]*round(dy/L[1]);
							if (dy < shellij){
								rij = sqrt(dx*dx + dy*dy);
								// if two tumor cells, check for attraction
								if (rij < shellij && rij > cutij && gi < NVtumor && gj < NVtumor){
									// force scale
									ftmp 				= kint*((rij/sij) - 1.0 - l2)/sij;
									fx 					= ftmp*(dx/rij);
									fy 					= ftmp*(dy/rij);

									// add to forces
									vF[NDIM*gi] 		-= fx;
									vF[NDIM*gi + 1] 	-= fy;

									vF[NDIM*gj] 		+= fx;
									vF[NDIM*gj + 1] 	+= fy;

									// add to virial expression for pressure
									pcheck += dx*fx + dy*fy;

									// add to contact network
									if (gj > gi)
										gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 1;
									else
										gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2] = 1;
								}
								else if (rij < cutij && rij > sij && gi < NVtumor && gj < NVtumor){
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

									// add to contact network
									if (gj > gi)
										gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 1;
									else
										gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2] = 1;
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

									// add to contact network
									if (gj > gi)
										gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 1;
									else
										gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2] = 1;
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
		psiMean = 0.0;
		psiStd = 0.0;
		for (gi=0; gi<NVTOT; gi++){

			// -- Area force (and get cell index ci)
			if (ci < NCELLS){
				if (gi == szList[ci]){
					// compute shape parameter
					nvtmp = nv[ci];
					a0tmp = a0[ci];
					l0tmp = l0[ci];

					/// compute area deviation
					atmp = area(vpos,ci,L,nv,szList);
					da = (atmp/a0tmp) - 1.0;

					// shape force parameters (kl and kl are unitless energy ratios)
					fa = ka*da*(rho0/a0tmp);
					fl = kl*(rho0/l0tmp);
					fb = kb*(rho0/(l0tmp*l0tmp));

					// compute cell center of mass
					xi = vpos[NDIM*gi];
					yi = vpos[NDIM*gi + 1];
					cx = xi; 
					cy = yi;
					for (vi=1; vi<nvtmp; vi++){
						dx = vpos[NDIM*(gi+vi)] - xi;

						dy = vpos[NDIM*(gi+vi) + 1] - yi;
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

					rim1y = vpos[NDIM*im1[gi] + 1] - cy;
					rim1y -= L[1]*round(rim1y/L[1]);

					rim2x = vpos[NDIM*im1[im1[gi]]] - cx;

					rim2y = vpos[NDIM*im1[im1[gi]] + 1] - cy;
					rim2y -= L[1]*round(rim2y/L[1]);

					// if adipocyte, get distance to pin
					if (ci >= tN){
						// get pin indices
						pinCellInd = ci - tN;

						// get distance to pin if attached
						if (pinattached[pinCellInd] == 1){
							dpinx = watpins[NDIM*pinCellInd] - cx;
							dpiny = watpins[NDIM*pinCellInd + 1] - cy;
							dpiny -= L[1]*round(dpiny/L[1]);
							dpin = sqrt(pow(dpinx,2.0) + pow(dpiny,2.0));
						}
					}
					// otherwise, update director
					else{
						// generate random variable
						r1 = drand48();
						r2 = drand48();
						grv = sqrt(-2.0*log(r1))*cos(2.0*PI*r2);

						// update director
						psi[ci] += sqrt(dt*2.0*DrList[ci])*grv;

						// add to mean and std dev of directors
						psiMean += psi[ci];
						psiStd += psi[ci]*psi[ci];
					}

					// boundary forces
					for (vi=0; vi<nvtmp; vi++){
						// positions using global indexing
						xind = NDIM*(gi+vi);
						yind = xind + 1;

						xi = vpos[xind];
						ri = vrad[gi+vi];

						// if near a wall, add to force
						if (xi < ri)
							vF[xind] += eint*(1.0 - (xi/ri))/ri;
						else if (xi > L[0] - ri)
							vF[xind] -= eint*(1.0 - ((L[0] - xi)/ri))/ri;

						// if an adipocyte and pin is intact, compute force due to pinning spring
						if (ci >= tN && pinattached[pinCellInd] == 1){
							if (dpin < pinbreak*sqrt(a0tmp)){
								vF[xind] 	+= (kpin/nvtmp)*dpinx;
								vF[yind] 	+= (kpin/nvtmp)*dpiny;
							}
							else
								pinattached[pinCellInd] = 0;
						}
					}

					// increment cell index
					ci++;
				}
			}


			// get next adjacent vertices
			rip1x = vpos[NDIM*ip1[gi]] - cx;

			rip1y = vpos[NDIM*ip1[gi] + 1] - cy;
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

			// bending force for adipocytes
			if (gi >= NVtumor){
				// segment vectors for ip2
				rip2x = vpos[NDIM*ip1[ip1[gi]]] - cx;

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
			// active force (only if tumor cell vertices)
			else {
				// get angular distance from psi
				psitmp = atan2(riy,rix);
				dpsi = psitmp - psi[ci-1];
				dpsi -= 2.0*PI*round(dpsi/(2.0*PI));

				// get velocity scale
				v0tmp = vmin + (v0-vmin)*exp(-pow(dpsi,2.0)/(2.0*Ds*Ds));

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

		// compute psi stats for output
		psiMean /= tN;
		psiStd /= tN;
		psiStd -= psiMean*psiMean;
		psiStd = sqrt(psiStd);

		// BROWNIAN MOTION VIA EULER UPDATE
		for (i=0; i<vertDOF; i++)
			vpos[i] += dt*vF[i];

		// update Dr based on contact with adipocytes
		gi = 0;
		for (ci=0; ci<tN; ci++){
			// compute number of tumor-adipocyte contacts
			zta = 0;
			for (vi=0; vi<nv[ci]; vi++){
				for (gj=NVtumor; gj<NVTOT; gj++){
					if (gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2])
						zta += 1.0;
				}
				gi++;
			}

			// scale Dr based on tumor-adipocyte contacts
			Drtmp = Dr*(1 - (zta/nv[ci])*dDr);
			if (Drtmp > Drmin)
				DrList[ci] = Drtmp;
			else
				DrList[ci] = Drmin;
		}

		// print message console, print position to file
		if (tt % NACTIVESKIP == 0){
			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << "			active tumor cells 				" << endl;
			cout << "===========================================" << endl;
			cout << endl;
			cout << "	** tt 			= " << tt << endl;
			cout << "	** p 			= " << pcheck << endl;
			cout << "	** psiMean 		= " << psiMean << endl;
			cout << "	** psiStd 		= " << psiStd << endl;

			// compute instantaneous packing fraction
			phi = 0.0;
			for (ci=0; ci<NCELLS; ci++)
				phi += area(vpos,ci,L,nv,szList) + 0.25*PI*pow(l0[ci]*del,2.0)*(0.5*nv[ci] - 1);
			phi /= L[0]*L[1];

			// print vertex positions to check placement
			cout << "\t** PRINTING POSITIONS TO FILE... " << endl;
			printPos(posout, vpos, vF, a0, l0, L, nv, szList, phi, NCELLS, tN);
		}
	}





	// close file objects
	posout.close();

	// end
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
	for (vi=0; vi<nv[ci]; vi++){
		// next vertex
		vip1 = (vi + 1) % nv[ci];
		gip1 = gindex(ci,vip1,szList);

		// get positions (check minimum images)
		dx 		= vpos[NDIM*gip1] - xi;
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
	for (vi=0; vi<nv[ci]; vi++){
		// next vertex
		vip1 = (vi + 1) % nv[ci];
		gip1 = gindex(ci,vip1,szList);

		// get positions (check minimum images)
		dx 		= vpos[NDIM*gip1] - xi;
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
void printPos(ofstream& posout, vector<double>& vpos, vector<double>& vvel, vector<double>& a0, vector<double>& l0, vector<double>& L, vector<int>& nv, vector<int>& szList, double phi, int NCELLS, int tN){
	// local variables
	int ci, vi, gi;
	double xi, yi, vxi, vyi, dx, dy, Lx, Ly;

	// save box sizes
	Lx = L.at(0);
	Ly = L.at(1);

	// print information starting information
	posout << setw(w) << left << "NEWFR" << " " << endl;
	posout << setw(w) << left << "NUMCL" << setw(w) << right << NCELLS << setw(w) << right << tN << endl;
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

		// get vertex velocities
		vxi = vvel.at(NDIM*gi);
		vyi = vvel.at(NDIM*gi + 1);

		// place back in box center
		xi = fmod(xi,Lx);
		yi = fmod(yi,Ly);

		posout << setw(w) << left << "VINFO";
		posout << setw(w) << left << ci;
		posout << setw(w) << left << 0;

		// output initial vertex information
		posout << setw(wnum) << setprecision(pnum) << right << xi;
		posout << setw(wnum) << setprecision(pnum) << right << yi;
		posout << setw(wnum) << setprecision(pnum) << right << vxi;
		posout << setw(wnum) << setprecision(pnum) << right << vyi;
		posout << endl;

		// vertex information for next vertices
		for (vi=1; vi<nv.at(ci); vi++){
			// get global vertex index for next vertex
			gi++;

			// get next vertex positions (use MIC)
			dx = vpos.at(NDIM*gi) - xi;
			xi += dx;

			dy = vpos.at(NDIM*gi + 1) - yi;
			dy -= Ly*round(dy/Ly);
			yi += dy;

			// get vertex velocities
			vxi = vvel.at(NDIM*gi);
			vyi = vvel.at(NDIM*gi + 1);

			// Print indexing information
			posout << setw(w) << left << "VINFO";
			posout << setw(w) << left << ci;
			posout << setw(w) << left << vi;

			// output vertex information
			posout << setw(wnum) << setprecision(pnum) << right << xi;
			posout << setw(wnum) << setprecision(pnum) << right << yi;
			posout << setw(wnum) << setprecision(pnum) << right << vxi;
			posout << setw(wnum) << setprecision(pnum) << right << vyi;
			posout << endl;
		}
	}

	// print end frame
	posout << setw(w) << left << "ENDFR" << " " << endl;
}



