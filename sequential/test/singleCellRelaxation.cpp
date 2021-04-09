
/* 

	MAIN FILE FOR FIRE RELAXATION 
		OF SINGLE DP, DPb and DPbb cells

	Jack Treado
	10/28/2020, in the time of covid

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

// namespaces
using namespace Eigen;
using namespace std;

// GLOBAL CONSTANTS
const int NCELLS 			= 1;
const double PI 			= 4*atan(1);
const int w 				= 10;
const int wnum 				= 25;
const int pnum 				= 14;

// simulation constants
const double timeStepMag 	= 0.01;
const double pertScale 		= 0.01;
const double phi0 			= 0.01;

// FIRE constants for initial minimizations (SP + DP)
const double alpha0      	= 0.2;
const double finc        	= 1.1;
const double fdec        	= 0.5;
const double falpha      	= 0.99;
const double Ftol 			= 1e-14;

const int NSKIP 			= 1e4;
const int NMIN        		= 10;
const int NNEGMAX     		= 1000;
const int NDELAY      		= 50;
const int itmax       		= 1e7;


// DP force constants
const double ka 			= 1.0;			// area spring (should be = 1)
const double del 			= 1.0;			// radius of vertices in units of l0



// FUNCTION PROTOTYPES
int gindex(int ci, int vi, vector<int>& szList);
void cindices(int& ci, int& vi, int gi, int NCELLS, vector<int>& szList);

double area(vector<double>& dpos, int ci, vector<double>& L, vector<int>& nv, vector<int>& szList);
double perimeter(vector<double>& dpos, int ci, vector<double>& L, vector<int>& nv, vector<int>& szList);

void printPos(ofstream& posout, vector<double>& vpos, vector<double>& a0, vector<double>& l0, vector<double>& L, vector<int>& nv, vector<int>& szList, double phi, int NCELLS);

// MAIN
int main(int argc, char const *argv[]){
	// variables for indexing loops
	int i, ci, cj, vi, vj, gi, gj, d;

	// parameters to be read in 
	int NV, NT, vertDOF, seed;
	double kl, kb, kbb, calA0;

	// read in parameters from command line input
	string NV_str 			= argv[1];
	string calA0_str 		= argv[2];
	string kl_str 			= argv[3];
	string kb_str 			= argv[4];
	string kbb_str 			= argv[5];
	string seed_str 		= argv[6];
	string positionFile 	= argv[7];
	string vdosFile 		= argv[8];

	stringstream NVss(NV_str);
	stringstream calA0ss(calA0_str);
	stringstream klss(kl_str);
	stringstream kbss(kb_str);
	stringstream kbbss(kbb_str);
	stringstream seedss(seed_str);

	NVss >> NV;
	calA0ss >> calA0;
	klss >> kl;
	kbss >> kb;
	kbbss >> kbb;
	seedss >> seed;

	// open xyz file
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

	// szList and nv (keep track of global vertex indices)
	cout << "makin vectors" << endl;
	vector<int> szList(1,0);
	vector<int> nv(1,0);
	nv.at(0) = NV;

	// degree of freedom counts
	vertDOF = NDIM*NV;

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
	calA0 *= NV*tan(PI/NV)/PI;

	cout << "getting times" << endl;

	// fundamental MD time units
	double dtMD, dt0, dt;

	dtMD 	= 1.0;
	dt0 	= timeStepMag*dtMD;
	dt 		= dt0;

	// output opening statement to console
	cout << "=======================================================" << endl << endl;

	cout << "		singleCellRelaxation.cpp 						" << endl;
	cout << "		Jack Treado, 2020   							" << endl;
	cout << "		Single cell relaxation using FIRE 				" << endl << endl;

	cout << "       NV 			= " << NV << "						" << endl; 
	cout << "		calA0 		= " << calA0 << "					" << endl << endl;

	cout << "		kl 			= " << kl << "						" << endl;
	cout << "		kb 			= " << kb << "						" << endl;
	cout << "		kbb 		= " << kbb << " 					" << endl;
	cout << "		seed 		= " << seed << "					" << endl << endl;

	cout << "		pos file 	= " << positionFile << "			" << endl;
	cout << "		vdos file 	= " << vdosFile << " 				" << endl << endl;
	
	cout << "=======================================================" << endl << endl;

	// seed random number generator
	srand48(seed);


	/* * * * * * * * * * * * * * * * * * 

				INITIALIZE

					PARTICLE

	 * * * * * * * * * * * * * * * * * */

	// initialize vectors for storing coordinates, shape information
	vector<double> vrad(NV,1.0);
	vector<double> vpos(vertDOF,0.0);
	double a0, l0;

	// store preferred area
	a0 			= 1.0;

	// set l0, vector radius
	l0 	= 2.0*sqrt(PI*calA0)/NV;
	for (vi=0; vi<NV; vi++)
		vrad.at(vi)	= 0.5*l0.at(ci)*del;

	// determine box lengths from particle sizes and input packing fraction
	vector<double> L(NDIM,1.0);
	for (d=0; d<NDIM; d++)
		L.at(d) = sqrt(a0/phi0);


	// initialize particle positions
	for (vi=0; vi<NV; vi++){
		// get distance from cell center to vertex if reg poly
		lenscale = sqrt((2.0*a0)/(NV*sin((2.0*PI)/NV)));

		// set vertex positions
		vpos.at(NDIM*vi) 		= lenscale*cos((2.0*PI*vi)/NV) + 0.5*L[0] + pertScale*l0*drand48();
		vpos.at(NDIM*vi + 1)	= lenscale*sin((2.0*PI*vi)/NV) + 0.5*L[1] + pertScale*l0*drand48();
	}

	// define length unit
	double rho0 = sqrt(a0);


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
	dtmin   	= 1e-6*dt0;

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
	double Kl, Kb, Kbb, l0, atmp, li, lim1, kappai, cx, cy;
	double da, dli, dlim1;
	double lim2x, lim2y, lim1x, lim1y, lix, liy, lip1x, lip1y;
	double rim2x, rim2y, rim1x, rim1y, rix, riy, rip1x, rip1y, rip2x, rip2y;

	// compute Dc for belt spring
	double Dc = l0/sin(PI/NV);

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


		// FORCE UPDATE

		// shape forces (loop over global vertex labels)
		for (vi=0; vi<NV; vi++){

			// -- Area force (and get cell index ci)
			if (vi == 0){
				// compute area deviation
				atmp = area(vpos,0,L,nv,szList);
				da = (atmp/a0) - 1.0;

				// compute dL for belt spring (use floor so, if NV odd, one vertex does not get force, but all forces are balanced)
				dL = 0.0;
				for (vj=0; vi<floor(0.5*NV); vi++){
					vop = vj + floor(0.5*NV);
					dx = vpos.at(NDIM*vop) - vpos.at(NDIM*vj);
					dy = vpos.at(NDIM*vop + 1) - vpos.at(NDIM*vj + 1);
					Ltmp = sqrt(dx*dx + dy*dy);
					dL += Ltmp/floor(0.5*NV);
				}
				dL = (dL/Dc) - 1.0;

				// shape force parameters (kl and kl are unitless energy ratios)
				fa = da*(rho0/a0);		// derivation from the fact that rho0^2 does not necessarily cancel a0
				fl = NV*kl*(rho0/l0);
				fb = NV*kb*(rho0/(l0*l0));
				fbb = kbb*dL*(2.0*rho0/(NV*Dc));
				
				// compute cell center of mass
				xi = vpos[NDIM*vi];
				yi = vpos[NDIM*vi + 1];
				cx = xi; 
				cy = yi;
				for (vj=1; vj<NV; vj++){
					dx = vpos.at(NDIM*vi) - xi;
					dx -= L[0]*round(dx/L[0]);

					dy = vpos.at(NDIM*vi + 1) - yi;
					dy -= L[1]*round(dy/L[1]);

					xi += dx;
					yi += dy;

					cx += xi;
					cy += yi;
				}
				cx /= NV;
				cy /= NV;

				// get coordinates relative to center of mass
				rix = vpos[NDIM*vi] - cx;
				riy = vpos[NDIM*vi + 1] - cy;

				// get (prior) adjacent vertices
				rim1x = vpos[NDIM*im1[vi]] - cx;
				rim1x -= L[0]*round(rim1x/L[0]);

				rim1y = vpos[NDIM*im1[vi] + 1] - cy;
				rim1y -= L[1]*round(rim1y/L[1]);

				rim2x = vpos[NDIM*im1[im1[vi]]] - cx;
				rim2x -= L[0]*round(rim2x/L[0]);

				rim2y = vpos[NDIM*im1[im1[vi]] + 1] - cy;
				rim2y -= L[1]*round(rim2y/L[1]);
			}


			// get next adjacent vertices
			rip1x = vpos.at(NDIM*ip1[vi]) - cx;
			rip1x -= L[0]*round(rip1x/L[0]);

			rip1y = vpos.at(NDIM*ip1[vi] + 1) - cy;
			rip1y -= L[1]*round(rip1y/L[1]);



			// -- Area force
			vF[NDIM*vi] 		+= 0.5*fa*(rim1y - rip1y);
			vF[NDIM*vi + 1] 	+= 0.5*fa*(rip1x - rim1x);


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
			dlim1  	= (lim1/l0) - 1.0;
			dli 	= (li/l0) - 1.0;

			// add to forces
			vF[NDIM*vi] 		+= fl*(dli*(lix/li) - dlim1*(lim1x/lim1));
			vF[NDIM*vi + 1] 	+= fl*(dli*(liy/li) - dlim1*(lim1y/lim1));


			// -- Bending force
			if (kb > 0){
				// segment vectors for ip2
				rip2x = vpos[NDIM*ip1[ip1[vi]]] - cx;
				rip2x -= L[0]*round(rip2x/L[0]);

				rip2y = vpos[NDIM*ip1[ip1[vi]] + 1] - cy;
				rip2y -= L[1]*round(rip2y/L[1]);

				lip1x = rip2x - rip1x;
				lip1y = rip2y - rip1y;

				lim2x = rim1x - rim2x;
				lim2y = rim1y - rim2y;

				// add to force
				vF[NDIM*vi] 		+= fb*(3.0*(lix - lim1x) + lim2x - lip1x);
				vF[NDIM*vi + 1] 	+= fb*(3.0*(liy - lim1y) + lim2y - lip1y);
			}


			// -- Belt spring force
			if (kbb > 0 && vi < floor(0.5*NV)){
				// get index of opposite vertex
				vop 				= vi + floor(0.5*NV);

				// get distances to opposite vertices
				dx 					= vpos[NDIM*vop] - vpos[NDIM*vj];
				dy 					= vpos[NDIM*vop + 1] - vpos[NDIM*vj + 1];
				Ltmp 				= sqrt(dx*dx + dy*dy);

				// add to forces on vi
				vF[NDIM*vi]			+= fbb*(dx/Ltmp);
				vF[NDIM*vi + 1] 	+= fbb*(dy/Ltmp);

				// add to forces on opposite vertex
				vF[NDIM*vop]		-= fbb*(dx/Ltmp);
				vF[NDIM*vop + 1] 	-= fbb*(dy/Ltmp);
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


	// print vertex positions to check placement
	cout << "\t** PRINTING POSITIONS TO FILE... " << endl;
	printPos(posout, vpos, a0, l0, L, nv, szList, phi0, NCELLS);

	


























	/* * * * * * * * * * * * * * * * * * 

			   COMPUTE VDOS

	 * * * * * * * * * * * * * * * * * */


	// LOCAL VARIABLES

	// integers
	int l, lxm1, lym1, lxp1, lyp1;
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
	Eigen::MatrixXd Hvv(vertDOF,vertDOF);		// stiffness matrix for belt spring
	Eigen::MatrixXd Sbb(vertDOF,vertDOF);		// stress matrix for belt spring
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
			Hvv(k,l) = 0.0;
			Svv(k,l) = 0.0;
			S(k,l) = 0.0;
			H(k,l) = 0.0;
			M(k,l) = 0.0;
		}
	}

	// Loop over cells, compute shape forces for each individual cell and contributions from
	// vertex-vertex interactions
	rho0 = sqrt(a0.at(0));

	// print statement
	cout << "		-- Computing dynamical matrix elements for single cell" << endl;
		

	// ------------------------------------------
	// 
	// 				SHAPE
	// 					CONTRIBUTIONS
	//
	// ------------------------------------------

	// area deviations
	delA = (area(vpos,0,L,nv,szList)/a0) - 1.0;

	// dimensionless stiffness constants
	Kl1 = kl*(rho0*rho0)/l0;		// units = L
	Kl2 = Kl1/l0;					// units = 1

	Kb1 = kb*(rho0*rho0);			// units = L^2
	Kb2 = Kb1/(l0*l0);				// units = 1



	// loop over vertices, compute each DM element
	for (vi=0; vi<NV; vi++){

		// wrap vertices
		vim2 		= (vi - 2 + NV) % NV;
		vim1 		= (vi - 1 + NV) % NV;
		vip1 		= (vi + 1) % NV;
		vip2 		= (vi + 2) % NV;			

		// vertex elements
		kxm2 		= NDIM*(vim2);
		kym2 		= NDIM*(vim2) + 1;

		kxm1 		= NDIM*(vim1);
		kym1 		= NDIM*(vim1) + 1;

		kx 			= NDIM*(vi);
		ky 			= NDIM*(vi) + 1;

		kxp1 		= NDIM*(vip1);
		kyp1 		= NDIM*(vip1) + 1;

		kxp2 		= NDIM*(vip2);
		kyp2 		= NDIM*(vip2) + 1;

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
		deli 		= (li/l0) - 1.0;
		delim1 		= (lim1/l0) - 1.0;


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
		kapim1 			= sqrt(pow(lim1x - lim2x,2.0) + pow(lim1y - lim2y,2.0))/l0;
		kapi 			= sqrt(pow(lix - lim1x,2.0) + pow(liy - lim1y,2.0))/l0;
		kapip1 			= sqrt(pow(lip1x - lix,2.0) + pow(lip1y - liy,2.0))/l0;

		// curvature derivatives

		// derivatives of kapim1
		dkapim1_dxi 	= (lim1x - lim2x)/(kapim1*l0*l0);
		dkapim1_dyi 	= (lim1y - lim2y)/(kapim1*l0*l0);

		// derivatives of kapi
		dkapi_dxip1 	= (lix - lim1x)/(kapi*l0*l0);
		dkapi_dyip1 	= (liy - lim1y)/(kapi*l0*l0);
		dkapi_dxi 		= -2.0*dkapi_dxip1;
		dkapi_dyi 		= -2.0*dkapi_dyip1;	

		// derivatives of kapip1
		dkapip1_dxi 	= (lip1x - lix)/(kapip1*l0*l0);
		dkapip1_dyi 	= (lip1y - liy)/(kapip1*l0*l0);
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
	    Sb(kx,kx)       = Kb2*(6.0 - (l0*dkapim1_dxi)*(l0*dkapim1_dxi) - (l0*dkapi_dxi)*(l0*dkapi_dxi) - (l0*dkapip1_dxi)*(l0*dkapip1_dxi));
	    Sb(ky,ky)       = Kb2*(6.0 - (l0*dkapim1_dyi)*(l0*dkapim1_dyi) - (l0*dkapi_dyi)*(l0*dkapi_dyi) - (l0*dkapip1_dyi)*(l0*dkapip1_dyi));
	    
	    Sb(kx,ky)       = -Kb1*(dkapim1_dxi*dkapim1_dyi + dkapi_dxi*dkapi_dyi + dkapip1_dxi*dkapip1_dyi);
	    Sb(ky,kx)       = Sb(kx,ky);
	    
	    // 1off block diagonal
	    Sb(kx,kxp1)     = -2.0*Kb2*(2.0 - (l0*dkapi_dxip1)*(l0*dkapi_dxip1) - (l0*dkapip1_dxi)*(l0*dkapip1_dxi));
	    Sb(ky,kyp1)     = -2.0*Kb2*(2.0 - (l0*dkapi_dyip1)*(l0*dkapi_dyip1) - (l0*dkapip1_dyi)*(l0*dkapip1_dyi));
	    
	    Sb(kx,kyp1)     = -Kb1*(dkapi_dxi*dkapi_dyip1 + dkapip1_dxi*dkapip1_dyip1);
	    Sb(ky,kxp1)     = -Kb1*(dkapi_dyi*dkapi_dxip1 + dkapip1_dyi*dkapip1_dxip1);
	    
	    // 2off block diagonal
	    Sb(kx,kxp2)     = Kb2*(1.0 - (l0*dkapip1_dxi)*(l0*dkapip1_dxi));
	    Sb(ky,kyp2)     = Kb2*(1.0 - (l0*dkapip1_dyi)*(l0*dkapip1_dyi));
	    
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
	    if (kbb > 0 && vi < floor(0.5*NV)){
	    	// get index of opposite vertex
			vop 				= vi + floor(0.5*NV);

			// get matrix element indices for opposite
			kopx 				= NDIM*vop;
			kopy 				= kopx + 1;

			// get distances to opposite vertices
			dx 					= vpos[NDIM*vop] - vpos[NDIM*vj];
			dy 					= vpos[NDIM*vop + 1] - vpos[NDIM*vj + 1];
			Ltmp 				= sqrt(dx*dx + dy*dy);

			uLx 				= dx/Ltmp;
			uLy 				= dy/Ltmp;



			// STIFFNESS MATRIX

			// block diagonal
			Hbb(kx,kx) 			= ((4.0*kbb*a0)/(pow(NV*Dc,2.0)))*(dx/Ltmp)*uLx*uLx;
			Hbb(ky,ky) 			= ((4.0*kbb*a0)/(pow(NV*Dc,2.0)))*(dx/Ltmp)*uLy*uLy;
			Hbb(kx,ky) 			= ((4.0*kbb*a0)/(pow(NV*Dc,2.0)))*(dx/Ltmp)*uLx*uLy;
			Hbb(ky,kx) 			= Hbb(kx,ky);
			
			// off diagonals
	        Hbb(kx,kopx)       = -Hbb(kx,kx);
        	Hbb(ky,kopy)       = -Hbb(ky,ky);
        	Hbb(kx,kopy)       = -Hbb(kx,ky);
        	Hbb(ky,kopx)       = -Hbb(ky,kx);
        
        	Hbb(kopx,kx)       = Hbb(kx,kopx);
        	Hbb(kopy,ky)       = Hbb(ky,kopy);
        	Hbb(kopx,ky)       = Hbb(ky,kopx);
        	Hbb(kopy,kx)       = Hbb(kx,kopy);
        
        	// block diagonal (vop)
        	Hbb(kopx,kopx)    = Hbb(kx,kx);
        	Hbb(kopy,kopy)    = Hbb(ky,ky);
        	Hbb(kopx,kopy)    = Hbb(kx,ky);
        	Hbb(kopy,kopx)    = Hbb(ky,kx);


        	// -- stress elements
        
	        // block diagonal
	        Sbb(kx,kx)          = ((2.0*kbb*dL*a0)/(NV*Dc*Ltmp))*uLy*uLy;
	        Sbb(ky,ky)          = ((2.0*kbb*dL*a0)/(NV*Dc*Ltmp))*uLx*uLx;
	        Sbb(kx,ky)          = -(()2.0*kbb*dL*a0)/(NV*Dc*Ltmp))*uLx*uLy;
	        Sbb(ky,kx)          = Sbb(kx,ky);
	        
	        // off diagonals
	        Sbb(kx,kopx)       = -Sbb(kx,kx);
	        Sbb(ky,kopy)       = -Sbb(ky,ky);
	        Sbb(kx,kopy)       = -Sbb(kx,ky);
	        Sbb(ky,kopx)       = -Sbb(ky,kx);
	        
	        Sbb(kopx,kx)       = Sbb(kx,kopx);
	        Sbb(kopy,ky)       = Sbb(ky,kopy);
	        Sbb(kopx,ky)       = Sbb(ky,kopx);
	        Sbb(kopy,kx)       = Sbb(kx,kopy);
	        
	        // block diagonal (vop)
	        Sbb(kopx,kopx)    = Sbb(kx,kx);
	        Sbb(kopy,kopy)    = Sbb(ky,ky);
	        Sbb(kopx,kopy)    = Sbb(kx,ky);
	        Sbb(kopy,kopx)    = Sbb(ky,kx);
	    }














		

	    // -- AREA SPRING (stress matrix)
	    Sa(kx,kyp1) = 0.5*delA*((rho0*rho0)/a0);
		Sa(ky,kxp1) = -0.5*delA*((rho0*rho0)/a0);

		Sa(kyp1,kx) = Sa(kx,kyp1);
		Sa(kxp1,ky) = Sa(ky,kxp1);

		// area derivatives (for stiffness matrix)
		da_dxi      = 0.5*(liy + lim1y);
		da_dyi      = -0.5*(lim1x + lix);

		// loop over other vertices, for area elasticity stiffness matrix
		for (vj=vi; vj<NV; vj++){

			// wrap jp1 and jm1
			vjp1 		= (vj + 1) % NV;
			vjm1 		= (vj - 1 + NV) % NV;

			// dof elements
			lxm1 		= NDIM*vjm1;
			lym1 		= lxm1 + 1;

			lx 			= NDIM*vj;
			ly 			= lx + 1;

			lxp1		= NDIM*vjp1;
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
			Ha(kx,lx) = da_dxi*da_dxj*((rho0*rho0)/pow(a0,2.0));
	        Ha(kx,ly) = da_dxi*da_dyj*((rho0*rho0)/pow(a0,2.0));
	        
	        Ha(ky,lx) = da_dyi*da_dxj*((rho0*rho0)/pow(a0,2.0));
	        Ha(ky,ly) = da_dyi*da_dyj*((rho0*rho0)/pow(a0,2.0));
	        
	        Ha(lx,kx) = Ha(kx,lx);
	        Ha(ly,kx) = Ha(kx,ly);
	        
	        Ha(lx,ky) = Ha(ky,lx);
	        Ha(ly,ky) = Ha(ky,ly);
		}
	}

	// compute D from sum of other dynamical matrices
	// initialize all matrices to be 0 initially
	for (k=0; k<vertDOF; k++){
		for (l=0; l<vertDOF; l++){
			H(k,l) = Ha(k,l) + Hl(k,l) + Hb(k,l) + Hvv(k,l);
			S(k,l) = -Sa(k,l) - Sl(k,l) - Sb(k,l) - Svv(k,l);
			M(k,l) = H(k,l) - S(k,l);
		}
	}

	// compute eigenvalues
	cout << "\t** Computing eigenvalues and eigenvectors of M, H and S matrices" << endl;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> allModes(M);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> hModes(H);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> sModes(S);

	// define eigenvector matrix
	Eigen::MatrixXd evecs = allModes.eigenvectors();

	// Use Graham-Schmidt Orthogonalization to rotate modes

	// print eigenvalues to vdos file
	cout << "\t** Printing evals to file" << endl;
	vdosout << vertDOF << endl;
	vdosout << allModes.eigenvalues() << endl;
	vdosout << hModes.eigenvalues() << endl;
	vdosout << sModes.eigenvalues() << endl;
	vdosout << evecs << endl;


	// close open objects
	posout.close();
	vdosout.close();


	// print to console, return
	cout << "\n\n\nFINISHED MAIN FOR singleCellRelaxation.cpp, ENDING." << endl << endl << endl;
	return 0;
}