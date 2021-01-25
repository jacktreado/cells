/* 

	MAIN FILE FOR 2D MESOPHYLL TISSUE

	NCELLS polydisperse disperse DPb particles

	Features:
		-- contact dependent adhesion
		-- crosslinking with binding kinetics
		-- perimeter aging
		-- quasistatic decompression

	ADDED 

	01/13/21
		-- RIGIDIFICATION: cells all have preferred curvature that lags
			behind instantaneous curvature, will rigidify shape
			as decompression progresses

	01/14/21
		-- LOCALIZED CONTRACTILITY: each perimeter spring will have an l0
			localized to it, rather than an l0 for each cell
		-- CONTACT-DEPENDENT AGING: The perimeter will age twice as fast
			on contacts that are engaged vs void contacts
		-- SHAPE LIMIT: There is an upper limit to cell preferred shape



	NOTE 01/24/21
		-- Look into effect of contactScale and voidScale, deformation vs void size scaling??
		-- Should these be variable inputs?

	Jack Treado
	11/03/2020, in the time of covid

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
const double timeStepMag 	= 0.01;
const double phiInit 		= 0.2;
const double dphiGrow 		= 0.01;
const double dphiShrink 	= 0.001;
const double dphiPrint 		= 0.01;

// FIRE constants for initial minimizations (SP + DP)
const double alpha0      	= 0.2;
const double finc        	= 1.1;
const double fdec        	= 0.5;
const double falpha      	= 0.99;
const double Ftol 			= 1e-12;

const int NSKIP 			= 1e4;
const int NMIN        		= 10;
const int NNEGMAX     		= 1000;
const int NDELAY      		= 20;
const int itmax       		= 1e7;

// DP force constants
const double ka 			= 1.0;			// area spring (should be = 1)
const double eint 			= 1.0;			// baseline interaction energy 
const double del 			= 1.0;			// radius of vertices in units of l0

// displacement magnitude for MC bonds
const double bondDisp 		= 1.0;

// shape aging constants
const double contactScale 	= 2.0;			// rate of growth rel. to lambdal of contact vertices
const double voidScale 		= 1.0;			// rate of growth rel. to lambdal of void-facing vertices
const double calA0Thresh 	= 2.0;			// max preferred shape parameter allowable

// FUNCTION PROTOTYPES

// indexing
int gindex(int ci, int vi, vector<int>& szList);
void cindices(int& ci, int& vi, int gi, int NCELLS, vector<int>& szList);

// particle geometry
double area(vector<double>& dpos, int ci, vector<double>& L, vector<int>& nv, vector<int>& szList);
double perimeter(vector<double>& dpos, int ci, vector<double>& L, vector<int>& nv, vector<int>& szList);

// system potential energy without spring network (for Metropolis choice)
double potentialEnergyNoNetwork(vector<double>& vpos, 
	vector<double>& vrad, 
	vector<double>& a0, 
	vector<double>& l0, 
	vector<double>& delta0,
	vector<double>& s0,
	vector<double>& L, 
	vector<int>& nv, 
	vector<int>& szList, 
	vector<int> im1, 
	vector<int> ip1, 
	double kl, 
	double kb, 
	int NCELLS);

// remove rattlers from contact network, return rattler number
int removeRattlers(vector<int>& cij);

// print to file
void printPos(ofstream& posout, vector<double>& vpos, vector<double>& vrad, vector<double>& a0, vector<double>& l0, vector<double>& L, vector<int>& cij, vector<int>& nv, vector<int>& szList, double phi, int NCELLS); 



// MAIN
int main(int argc, char const *argv[]){
	// variables for indexing loops
	int i, ci, cj, vi, vj, gi, gj, d;

	// parameters to be read in 
	int NCELLS, NV, NVTOT, cellDOF, vertDOF, seed;
	double polyd, calA0Input, phi, phiMax, phiMin, kl, kb, espring, lambdaL, lambdaB, betaEff;

	// set spring energy
	espring = eint;

	// read in parameters from command line input
	// test: g++ -O3 sequential/mesoRigidify2D.cpp -o meso.cpp
	// test: ./meso.o 16 32 1.01 0.1 1.05 0.4 1.0 0.01 0.01 0.01 1.0 1 pos.test
	string NCELLS_str 			= argv[1];
	string NV_str 				= argv[2];
	string calA0_str 			= argv[3];
	string polyd_str 			= argv[4];
	string phiMax_str  			= argv[5];
	string phiMin_str 			= argv[6];
	string kl_str 				= argv[7];
	string kb_str 				= argv[8];
	string lambdaL_str 			= argv[9];
	string lambdaB_str			= argv[10];
	string betaEff_str 			= argv[11];
	string seed_str 			= argv[12];
	string positionFile 		= argv[13];

	stringstream NCELLSss(NCELLS_str);
	stringstream NVss(NV_str);
	stringstream calA0ss(calA0_str);
	stringstream polydss(polyd_str);
	stringstream phiMaxss(phiMax_str);
	stringstream phiMinss(phiMin_str);
	stringstream klss(kl_str);
	stringstream kbss(kb_str);
	stringstream lambdaLss(lambdaL_str);
	stringstream lambdaBss(lambdaB_str);
	stringstream betaEffss(betaEff_str);
	stringstream seedss(seed_str);

	NCELLSss >> NCELLS;
	NVss >> NV;
	calA0ss >> calA0Input;
	polydss >> polyd;
	phiMaxss >> phiMax;
	phiMinss >> phiMin;
	klss >> kl;
	kbss >> kb;
	lambdaLss >> lambdaL;
	lambdaBss >> lambdaB;
	betaEffss >> betaEff;
	seedss >> seed;

	// open xyz file
	ofstream posout;
	posout.open(positionFile.c_str());
	if (!posout.is_open()){
		cout << "	** ERROR: position file " << positionFile << " could not be opened, ending." << endl;
		return 1;
	}

	// total number of vertices
	

	// szList and nv (keep track of global vertex indices)
	int nvtmp;
	double r1, r2, grv;
	vector<int> szList(NCELLS,0);
	vector<int> nv(NCELLS,0);
	nv.at(0) = NV;
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

		// print vertex info
		cout << "ci = " << ci << ";  nvtmp = " << nvtmp << endl;

		// store size of cell ci
		nv.at(ci) = nvtmp;
		szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);

		// add to total NV count
		NVTOT += nvtmp;
	}
	cout << "** minimum: at i = " << imin << ", rmin = " << rmin << endl;
	cout << "** maximum: at i = " << imax << ", rmax = " << rmax << endl;

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

	return 0;

	// fundamental MD time units
	double dtMD, dt0, dt;

	dtMD 	= 1.0;
	dt0 	= timeStepMag*dtMD;
	dt 		= dt0;

	// output opening statement to console
	cout << "=======================================================" << endl << endl;

	cout << "		mesophyll2D.cpp 								" << endl;
	cout << "		Jack Treado, 2020   							" << endl;

	cout << "		NCELLS 		= " << NCELLS << "					" << endl;
	cout << "       NV (all) 	= " << NV << "						" << endl;
	cout << "		NVTOT 		= " << NVTOT << "					" << endl << endl;

	cout << "		calA0 		= " << calA0Input << "				" << endl;

	cout << "		phiInit 	= " << phiInit << " 				" << endl;
	cout << "		phiMax 		= " << phiMax << " 					" << endl;
	cout << "		phiMin 		= " << phiMin << "					" << endl;
	cout << "		kl 			= " << kl << "						" << endl;
	cout << "		kb 			= " << kb << "						" << endl;
	cout << "		lambdaL 	= " << lambdaL << " 				" << endl;
	cout << "		lambdaB 	= " << lambdaB << " 				" << endl;
	cout << "		betaEff 	= " << betaEff << " 				" << endl;
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
	double lenscale, calA0tmp, a0tmp, areaSum = 0.0;

	// initialize vectors for storing coordinates, shape information
	vector<double> vrad(NVTOT,1.0);
	vector<double> drad(NCELLS,1.0);

	vector<double> vpos(vertDOF,0.0);
	vector<double> dpos(cellDOF,0.0);

	vector<double> a0(NCELLS,1.0);
	vector<double> l0(NCELLS,1.0);
	vector<double> calA0(NCELLS,1.0);

	// initialize effective disk radius (for minimization), and l0 parameter
	areaSum = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		// nv for cell ci
		nvtmp = nv.at(ci);

		// store preferred area
		lenscale 		= ((double)nvtmp/NV);
		a0tmp 			= lenscale*lenscale;
		a0.at(ci) 		= a0tmp;

		// scale calA0 tmp by calAv
		calA0tmp 		= calA0Input*(nvtmp*tan(PI/nvtmp)/PI);
		calA0.at(ci) 	= calA0tmp;

		// set disk radius
		drad.at(ci) 	= 1.05*sqrt((2.0*a0tmp)/(nvtmp*sin(2.0*PI/nvtmp)));

		// set l0 based on a0, calA0, vector radius
		l0.at(ci) 	= 2.0*sqrt(PI*calA0tmp*a0tmp)/nvtmp;
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
		L.at(d) = sqrt(areaSum/phiInit);
	phi = phiInit;

	// initialize tumor cells in center of packing
	for (ci=0; ci<NCELLS; ci++){
		dpos.at(NDIM*ci) = L[0]*drand48();
		dpos.at(NDIM*ci + 1) = L[1]*drand48();
	}


	// initialize cc contact network
	int NCTCS = 0.5*NCELLS*(NCELLS-1);
	vector<int> cij(NCTCS,0);


	// initialize vertex-vertex spring network
	int NVVCTCS = 0.5*NVTOT*(NVTOT-1);
	vector<bool> gij(NVVCTCS,0);


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
		sb[d] = round(L[d]/(3.0*l0.at(NCELLS-1)));

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

	// RIGIDIFY SIM: initialize preferred curvatures (unit of length)
	vector<double> s0(NVTOT,0.0);

	// jamming check variables
	int k, kmax, xind, yind;
	double pcheck;

	// max number of jamming iteractions
	k = 0;
	kmax = 1e4;

	// compute scale factors
	double scaleFactor = sqrt((phiInit + dphiGrow)/phiInit);

	// linked list variables
	int boxid, bi, bj, pi, pj, sbtmp;

	// temporary tolerance (to speed initial compression)
	double Ftoltmp = 1e3*Ftol;

	// total potential energy
	double U = 0.0;

	// length unit variable
	double rho0 = 0.0;

	// shape force variables
	double fa, fl, fb, l0tmp, atmp, li, lim1, cx, cy;
	double da, dli, dlim1;
	double lim2x, lim2y, lim1x, lim1y, lix, liy, lip1x, lip1y;
	double rim2x, rim2y, rim1x, rim1y, rix, riy, rip1x, rip1y, rip2x, rip2y;

	// contact variables
	int Nvv, Ncc;

	// compress to jamming, relax U and F using FIRE
	while (phi < phiMax && k < kmax){
		// update iterator
		k++;

		// update tolerance
		if (phi > 0.9*phiMax)
			Ftoltmp = Ftol;	

		// RESET FIRE VARIABLES
		P  			= 0;	
		fnorm 		= 0;
		vnorm 		= 0;
		alpha   	= alpha0;

		dtmax   	= 10*dt0;
		dtmin   	= 1e-8*dt0;
		dt 			= dt0;

		npPos      	= 0;
		npNeg      	= 0;
		npPMin      = 0;

		fireit    	= 0;
		fcheck  	= 10*Ftoltmp;

		// reset forces
		for (i=0; i<vertDOF; i++){
			vF[i] = 0.0;
			vFold[i] = 0.0;
			vvel[i] = 0.0;
		}

		// set length constant (variable due to particle growth)
		rho0 = sqrt(a0.at(0));

		// RELAX FORCES USING FIRE
		while ((fcheck > Ftoltmp || npPMin < NMIN) && fireit < itmax){
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
			fill(gij.begin(), gij.end(), false);

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

									if (ci != cj){
										if (ci > cj)
											cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
										else
											cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++;
										
										if (gi > gj)
											gij.at(NVTOT*gj + gi - (gj+1)*(gj+2)/2) = 1;
										else
											gij.at(NVTOT*gi + gj - (gi+1)*(gi+2)/2) = 1; 
									}

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

										// add to vertex-vertex contact network (ONLY IF DIFFERENT CELLS)
										if (ci != cj){
											if (ci > cj)
												cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
											else
												cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++;

											if (gi > gj)
												gij.at(NVTOT*gj + gi - (gj+1)*(gj+2)/2) = 1;
											else
												gij.at(NVTOT*gi + gj - (gi+1)*(gi+2)/2) = 1; 
										}
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

						// shape force parameters (kl and kl are unitless energy ratios)
						fa = ka*da*(rho0/a0tmp);		// derivation from the fact that rho0^2 does not necessarily cancel a0tmp
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
			if (fcheck < Ftoltmp)
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

		// update number of contacts
		Nvv = 0;
		Ncc = 0;
		for (i=0; i<NCTCS; i++){
			Nvv += cij[i];
			if (cij[i] > 0)
				Ncc++;
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
		cout << "	* Nvv 			= " << Nvv << endl;
		cout << "	* Ncc 			= " << Ncc << endl;
		cout << endl;

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


	// print initial dense state
	cout << "\t** PRINTING POSITIONS TO FILE... " << endl << endl << endl;
	printPos(posout, vpos, vrad, a0, calA0, L, cij, nv, szList, phi, NCELLS);










	/* * * * * * * * * * * * * * * * * * 

			DECOMPRESS TO

				FORM CELL GEL

	 * * * * * * * * * * * * * * * * * */

	// aging parameters
	double sip1, sip1x, sip1y, si, six, siy, sim1, sim1x, sim1y;
	double s0tmp, d0tmp, p0tmp, meanSegLength, meanPrefCurv;
	int Nb, NbRmv;

	// vector of perimeter segment length
	vector<double> linst(NVTOT,0.0);
	vector<int> contactVert(NVTOT,0);
	vector<double> delta0(NVTOT,1.0);

	// vector of instantaneous curvatures
	vector<double> sinst(NVTOT,0.0);

	// new scale factor
	scaleFactor = sqrt((phi - dphiShrink)/phi);

	// set Ftoltmp to Ftol
	Ftoltmp = Ftol;

	// # of interparticle contacts for each cell
	int ctmp;
	bool gtmp;
	double zij;
	vector<int> z(NCELLS,0);
	for (ci=0; ci<NCELLS; ci++){
		for (cj=0; cj<NCELLS; cj++){
			if (ci != cj){
				if (ci > cj)
					ctmp = cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2];
				else
					ctmp = cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2];
			}
			else
				ctmp = 0;

			// add to vector of per-cell vv contacts
			z.at(ci) += ctmp;
		}
	}

	// MC variables
	double rdraw, dU, dUtot, poff;
	vector<double> altpos(vertDOF,0.0);

	// loop until phi is below phiMin
	double lastPrintPhi = phi;
	k = 0;
	while (phi > phiMin && k < kmax){
		// update iterator
		k++;

		// RESET FIRE VARIABLES
		P  			= 0;	
		fnorm 		= 0;
		vnorm 		= 0;
		alpha   	= alpha0;

		dtmax   	= 10*dt0;
		dtmin   	= 2e-2*dt0;
		dt 			= dt0;

		npPos      	= 0;
		npNeg      	= 0;
		npPMin      = 0;

		fireit    	= 0;
		fcheck  	= 10*Ftoltmp;

		// reset forces
		fill(vF.begin(), vF.end(), 0.0);
		fill(vFold.begin(), vFold.end(), 0.0);
		fill(vvel.begin(), vvel.end(), 0.0);

		// set length constant (variable due to particle growth)
		rho0 = sqrt(a0.at(0));

		// RELAX FORCES USING FIRE
		while ((fcheck > Ftoltmp || npPMin < NMIN) && fireit < itmax){
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




			// FORCE UPDATE

			// reset contact network
			fill(cij.begin(), cij.end(), 0);

			// REPULSIVE interaction forces (USE BOX LINKED LIST)
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
									else
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
										else
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


			// SPRING NETWORK FORCES
			fill(contactVert.begin(), contactVert.end(), 0);
			for (gi=0; gi<NVTOT; gi++){
					
				// loop over other vertices and check connections					
				for (gj=gi+1; gj<NVTOT; gj++){

					// connections
					gtmp = gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2];

					// if connected, compute spring force based on z_\mu\nu
					if (gtmp){

						// label gi and gj as vertices in contact
						contactVert[gi] = 1;
						contactVert[gj] = 1;

						// contact distance
						sij = vrad[gi] + vrad[gj];

						// get vertex-vertex distances
						dx = vpos[NDIM*gj] - vpos[NDIM*gi];
						dx -= L[0]*round(dx/L[0]);

						dy = vpos[NDIM*gj + 1] - vpos[NDIM*gi + 1];
						dy -= L[1]*round(dy/L[1]);

						// true distance
						rij = sqrt(dx*dx + dy*dy);

						// only compute forces if spring is extended
						if (rij > sij){
							// get cell indices for vertex gi
							cindices(ci, vi, gi, NCELLS, szList);

							// get cell indices for vertex gj
							cindices(cj, vj, gj, NCELLS, szList);

							// get zij
							zij = 0.5*(z[ci] + z[cj]) + 1.0;

							// force scale
							ftmp 				= (espring/zij)*(1 - (rij/sij))*(rho0/sij);
							fx 					= ftmp*(dx/rij);
							fy 					= ftmp*(dy/rij);

							// add to forces
							vF[NDIM*gi] 		-= fx;
							vF[NDIM*gi + 1] 	-= fy;

							vF[NDIM*gj] 		+= fx;
							vF[NDIM*gj + 1] 	+= fy;

							// increae potential energy
							U += 0.5*(espring/zij)*pow((1 - (rij/sij)),2.0);

							// add to virial expression for pressure
							pcheck += dx*fx + dy*fy;
						}
					}
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

						// update potential energy
						U += 0.5*ka*da*da;

						// shape force parameters (kl and kl are unitless energy ratios)
						fa = ka*da*(rho0/a0tmp);		// derivation from the fact that rho0^2 does not necessarily cancel a0tmp
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
				lim1x 		= rix - rim1x;
				lim1y 		= riy - rim1y;

				lix 		= rip1x - rix;
				liy 		= rip1y - riy;

				// segment lengths
				lim1 		= sqrt(lim1x*lim1x + lim1y*lim1y);
				li 			= sqrt(lix*lix + liy*liy);

				// update instantaneous segment length
				linst[gi] 	= li;

				// segment deviations
				dlim1  		= (lim1/l0tmp) - delta0[im1[gi]];
				dli 		= (li/l0tmp) - delta0[gi];

				// add to forces
				vF[NDIM*gi] 		+= fl*(dli*(lix/li) - dlim1*(lim1x/lim1));
				vF[NDIM*gi + 1] 	+= fl*(dli*(liy/li) - dlim1*(lim1y/lim1));

				// update potential energy
				U += 0.5*kl*dli*dli;


				// -- Bending force

				// update instantaneous curvature
				six 		= lix - lim1x;
				siy 		= liy - lim1y;
				si 			= sqrt(six*six + siy*siy);
				sinst[gi] 	= si;

				// segment vectors for ip2
				rip2x = vpos[NDIM*ip1[ip1[gi]]] - cx;
				rip2x -= L[0]*round(rip2x/L[0]);

				rip2y = vpos[NDIM*ip1[ip1[gi]] + 1] - cy;
				rip2y -= L[1]*round(rip2y/L[1]);

				lip1x = rip2x - rip1x;
				lip1y = rip2y - rip1y;

				lim2x = rim1x - rim2x;
				lim2y = rim1y - rim2y;

				// compute instantaneous curvatures
				sim1x = lim1x - lim2x;
				sim1y = lim1y - lim2y;
				sim1 = sqrt(sim1x*sim1x + sim1y*sim1y);

				sip1x = lip1x - lix;
				sip1y = lip1y - liy;
				sip1 = sqrt(sip1x*sip1x + sip1y*sip1y);

				// add to force
				vF[NDIM*gi] 		+= fb*((sim1x/sim1)*(s0[im1[gi]] - sim1) - 2.0*(six/si)*(s0[gi] - si) + (sip1x/sip1)*(s0[ip1[gi]] - sip1));
				vF[NDIM*gi + 1] 	+= fb*((sim1y/sim1)*(s0[im1[gi]] - sim1) - 2.0*(siy/si)*(s0[gi] - si) + (sip1y/sip1)*(s0[ip1[gi]] - sip1));

				// add to potential energy
				U += 0.5*(kb/(l0tmp*l0tmp))*pow(si - s0[gi],2.0);

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
			if (fcheck < Ftoltmp)
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
				cout << "	** total U 	= " << U << endl;

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

				// decrease time step if past initial delay
				if (fireit > NDELAY){
					// decrease time step 
					if (dt*fdec > dtmin)
						dt *= fdec;

					// reset alpha
					alpha = alpha0;
				}

				// take half step backwards, reset velocities
				for (i=0; i<vertDOF; i++){
					// take half step backwards
					vpos[i] -= 0.5*dt*vvel[i] + 0.25*dt*dt*vF[i];

					// reset vertex velocities
					vvel[i] = 0.0;
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

		// BOND MC AND CONTACT NETWORK UPDATE
		Nb 		= 0;
		NbRmv 	= 0;
		dUtot 	= 0.0;
		for (gi=0; gi<NVTOT; gi++){	
			cindices(ci, vi, gi, NCELLS, szList);				
			for (gj=gi+1; gj<NVTOT; gj++){
				// get cell index of gj
				cindices(cj, vj, gj, NCELLS, szList);
				if (ci != cj){
					// contact distance
					sij = vrad[gi] + vrad[gj];

					// get vertex-vertex distances
					dx = vpos[NDIM*gj] - vpos[NDIM*gi];
					dx -= L[0]*round(dx/L[0]);

					dy = vpos[NDIM*gj + 1] - vpos[NDIM*gi + 1];
					dy -= L[1]*round(dy/L[1]);

					// true distance
					rij = sqrt(dx*dx + dy*dy);

					// mean contact number
					zij = 0.5*(z[ci] + z[cj]) + 1.0;

					// get current contact status
					gtmp = gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2];

					// reconnect bond if vertices come into contact
					if (rij <= sij && !gtmp){
						// update number of bonds
						Nb++;

						// add connection
						gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 1;

						// add bond to global contact network
						if (ci > cj)
							cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
						else
							cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++; 
					}
					// if extended and bond exists, detach if distance exceeds threshold
					else if (rij > sij && gtmp){
						// construct proposed configuration with displaced vertices
						altpos = vpos;

						// displace vertices in alternative configuration
						altpos[NDIM*gi] 		-= (bondDisp*sij)*(dx/rij);
						altpos[NDIM*gi + 1] 	-= (bondDisp*sij)*(dy/rij);

						altpos[NDIM*gj] 		+= (bondDisp*sij)*(dx/rij);
						altpos[NDIM*gj + 1] 	+= (bondDisp*sij)*(dy/rij);

						// compute change in potential energy
						dU = potentialEnergyNoNetwork(altpos, vrad, a0, l0, delta0, s0, L, nv, szList, im1, ip1, kl, kb, NCELLS);
						dU -= U;

						// add in change to potential energy based on bond displacement 
						// 	-- note in prefactor, sij cancels out based on def of bondDisp vs delta in notes
						dU += ((2.0*bondDisp*espring)/zij)*(((bondDisp*sij + rij)/sij) - 1.0);

						// add to total possible change in U
						dUtot += dU;

						// remove if bond detaching decreases energy
						if (dU < 0){
							gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 0;
							NbRmv++;
						}
						// else, remove conditionally
						else{
							poff = exp(-betaEff*dU);
							rdraw = drand48();

							// detach
							if (poff > rdraw){
								gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 0;
								NbRmv++;
							}
							// else, keep and add bond to global contact network
							else{
								if (ci > cj)
									cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
								else
									cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++; 
							}
						}
					}
				}
			}
		}

		// update number of contacts
		Nvv = 0;
		Ncc = 0;
		for (i=0; i<NCTCS; i++){
			Nvv += cij[i];
			if (cij[i] > 0)
				Ncc++;
		}

		for (ci=0; ci<NCELLS; ci++){
			z[ci] = 0;
			for (cj=0; cj<NCELLS; cj++){
				if (ci != cj){
					if (ci > cj)
						ctmp = cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2];
					else
						ctmp = cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2];
				}
				else
					ctmp = 0;

				// add to vector of per-cell vv contacts
				if (ctmp > 0)
					z[ci]++;
			}
		}


		// SHAPE AGING
		fill(calA0.begin(), calA0.end(), 0.0);
		meanPrefCurv = 0.0;
		meanSegLength = 0.0;
		gi = 0;
		for (ci=0; ci<NCELLS; ci++){

			// loop over vertices, age curvature and perimeter rest values
			p0tmp = 0.0;
			for (vi=0; vi<nv[ci]; vi++){
				// update s0 (preferred curvature)
				s0tmp = s0[gi];
				s0[gi] += lambdaB*(sinst[gi] - s0tmp);

				// update mean preferred curvature
				meanPrefCurv += s0[gi];

				// perimeter tmp variables
				l0tmp = l0[ci];
				d0tmp = delta0[gi];

				// only grow if calA0 below threshold
				if (calA0[ci] < calA0Thresh && l0tmp*d0tmp < linst[gi]){
					// update delta0 (preferred segment length) based on if between contact verts
					// otherwise, treat as void contact
					if (contactVert[gi] == 1 && contactVert[ip1[gi]] == 1)					
						delta0[gi] += contactScale*lambdaL*((linst[gi]/l0tmp) - d0tmp);
					else
						delta0[gi] += voidScale*lambdaL*((linst[gi]/l0tmp) - d0tmp);
				}

				// update total preferred perimeter
				p0tmp += delta0[gi]*l0tmp;

				// update mean preferred segment length
				meanSegLength += delta0[gi]*l0tmp;

				// update global vertex index
				gi++;
			}

			// update preferred shape parameter calA0
			calA0[ci] = pow(p0tmp,2.0)/(4.0*PI*a0[ci]);
		}
		meanPrefCurv /= NVTOT;
		meanSegLength /= NVTOT;


		// output to console
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===============================================" << endl << endl;
		cout << " 	Q U A S I S T A T I C  				" << endl;
		cout << " 	  	I S O T R O P I C 				" << endl;
		cout << "			G E L A T I O N 			" << endl << endl;
		cout << "===============================================" << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* dr 			= " << scaleFactor << endl;
		cout << "	* phi 			= " << phi << endl;
		cout << "	* fcheck 		= " << fcheck << endl;
		cout << "	* pcheck 		= " << pcheck << endl;
		cout << "	* U 		 	= " << U << endl << endl;
		cout << "	* Contacts:" << endl;
		cout << "	* Nvv 			= " << Nvv << endl;
		cout << "	* Ncc 			= " << Ncc << endl << endl;
		cout << "	* Spring network: " << endl;
		cout << "	* dUtot 		= " << dUtot << endl;
		cout << "	* Nb 			= " << Nb << endl;
		cout << "	* NbCurr 		= " << Nb - NbRmv << endl;
		cout << "	* NbRmv 		= " << NbRmv << endl << endl;
		cout << "	* Shape aging: " << endl;
		cout << "	* meanPrefCurv 	= " << meanPrefCurv << endl;
		cout << "	* meanSegLength = " << meanSegLength << endl;
		cout << endl;

		if ((lastPrintPhi - phi) > dphiPrint){
			// print positions
			cout << "\t** PRINTING POSITIONS TO FILE... " << endl << endl << endl;
			printPos(posout, vpos, vrad, a0, calA0, L, cij, nv, szList, phi, NCELLS);

			// update last phi when printed, for next time
			lastPrintPhi = phi;
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
		scaleFactor = sqrt((phi - dphiShrink)/phi);
	}

	// print final time
	cout << "\t** PRINTING POSITIONS TO FILE... " << endl << endl << endl;
	printPos(posout, vpos, vrad, a0, calA0, L, cij, nv, szList, phi, NCELLS);

	// close open objects
	posout.close();

	// end simulation
	cout << "\n\n\nFINISHED MAIN FOR mesophyll2D.cpp, ENDING." << endl << endl << endl;
	return 0;
}










/* 

	&&&&&&&&&&&&&&&&&&&&&&& FUNCTION DEFINITIONS &&&&&&&&&&&&&&&&&&&&&

	FUNCTIONS DEFINED

	gindex 				: returns global vertex index (gi) given cell (ci) and local vertex index (vi)
	cindex 				: returns cell index (ci) given global vertex index (gi)

	area 				: returns area of cell ci
	perimeter 			: returns perimeter of cell ci

	potentialEnergy 	: returns total potential energy of entire system given any configuration

	removeRattlers		: remove all rattlers from a contact network

	printPos 			: output vertex positions to .pos file for processing and visualization

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




// TOTAL SYSTEM POTENTIAL ENERGY WITHOUT SPRING NETWORK
double potentialEnergyNoNetwork(vector<double>& vpos, 
	vector<double>& vrad, 
	vector<double>& a0, 
	vector<double>& l0, 
	vector<double>& delta0,
	vector<double>& s0,
	vector<double>& L, 
	vector<int>& nv, 
	vector<int>& szList, 
	vector<int> im1, 
	vector<int> ip1, 
	double kl, 
	double kb, 
	int NCELLS){

	// local variables
	int ci, vi, gi, gj, nvtmp;
	double atmp, a0tmp, l0tmp;
	double dx, dy, rij, ri, rj, sij;
	double xi, yi, cx, cy, da, dli;
	double rix, riy, rim1x, rim1y, rip1x, rip1y, lim1x, lim1y, lix, liy, li, six, siy, si;
	double U = 0.0;

	// compute variables from input
	int NVTOT = 0;
	for (ci=0; ci<NCELLS; ci++)
		NVTOT += nv[ci];

	// add up potential energy contributions
	ci = 0;
	for (gi=0; gi<NVTOT; gi++){
		// radius of i
		ri = vrad[gi];

		// interactions with other repulsive vertices
		for (gj=gi+1; gj<NVTOT; gj++){
			// only if not adjacent
			if (gj != im1[gi] && gj != ip1[gi]){
				// radius of j
				rj = vrad[gj];

				// contact distance
				sij = ri + rj;

				// vertex-vertex distance
				dx = vpos[NDIM*gj] - vpos[NDIM*gi];
				dx -= L[0]*round(dx/L[0]);
				if (dx < sij){
					dy = vpos[NDIM*gj + 1] - vpos[NDIM*gi + 1];
					dy -= L[1]*round(dy/L[1]);
					if (dy < sij){
						rij = sqrt(dx*dx + dy*dy);

						// if overlaps, add to potential
						if (rij < sij)
							U += 0.5*eint*pow((1 - (rij/sij)),2.0);
					}
				}
			}
		}

		// shape potential energy
		if (ci < NCELLS){
			if (gi == szList[ci]){
				// compute shape parameter
				nvtmp = nv[ci];
				l0tmp = l0[ci];
				a0tmp = a0[ci];

				// compute area deviation
				atmp = area(vpos,ci,L,nv,szList);
				da = (atmp/a0tmp) - 1.0;

				// add area energy
				U += 0.5*ka*da*da;
				
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

				// increment cell index
				ci++;
			}

			// get next adjacent vertices
			rip1x = vpos.at(NDIM*ip1[gi]) - cx;
			rip1x -= L[0]*round(rip1x/L[0]);

			rip1y = vpos.at(NDIM*ip1[gi] + 1) - cy;
			rip1y -= L[1]*round(rip1y/L[1]);

			// -- Perimeter force

			// segments
			lim1x 		= rix - rim1x;
			lim1y 		= riy - rim1y;

			lix 		= rip1x - rix;
			liy 		= rip1y - riy;

			// segment length (just li for energy)
			li 			= sqrt(lix*lix + liy*liy);
			dli 		= (li/l0tmp) - delta0[gi];

			// add to potential energy
			U += 0.5*kl*dli*dli;

			// -- Bending force

			// get instantaneous curvature
			six 		= lix - lim1x;
			siy 		= liy - lim1y;
			si 			= sqrt(six*six + siy*siy);

			// add to potential energy
			U += 0.5*(kb/(l0tmp*l0tmp))*pow(si - s0[gi],2.0);

			// update old coordinates for next vertex
			rim1x = rix;
			rix = rip1x;

			rim1y = riy;
			riy = rip1y;
		}
	}

	// return potential energy
	return U;
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


