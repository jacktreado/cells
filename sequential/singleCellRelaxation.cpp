
/* 

	MAIN FILE FOR THERMAL NVE SIMULATION OF NCELLS
		USING LINKED-LIST SPEED UP

	Jack Treado
	08/05/2020, in the time of covid

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
const double sizeRatio 		= 1.4;
const double sizeFraction 	= 0.5;

// FIRE constants for initial minimizations (SP + DP)
const double alpha0      	= 0.2;
const double finc        	= 1.1;
const double fdec        	= 0.5;
const double falpha      	= 0.99;
const double Ftol 			= 1e-10;

const int NSKIP 			= 1e3;
const int NMIN        		= 100;
const int NNEGMAX     		= 10000;
const int NDELAY      		= 1000;
const int itmax       		= 1e7;


// DP force constants
const double ka 			= 1.0;			// area spring (should be = 1)
const double eint 			= 1.0;			// interaction energy 
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
	int NCELLS, NT, smallNV, largeNV, smallN, largeN, NVTOT, NVSMALL, cellDOF, vertDOF, seed;
	double NT_dbl, phi0, T0, kl, kb, calA0, smallCalA0, largeCalA0;

	// read in parameters from command line input
	string NCELLS_str 		= argv[1];
	string smallNV_str 		= argv[2];
	string calA0_str 		= argv[3];
	string phi0_str 		= argv[4];
	string T0_str 			= argv[5];
	string NT_str 			= argv[6];
	string kl_str 			= argv[7];
	string kb_str 			= argv[8];
	string seed_str 		= argv[9];
	string positionFile 	= argv[10];
	string energyFile 		= argv[11];

	stringstream NCELLSss(NCELLS_str);
	stringstream smallNVss(smallNV_str);
	stringstream calA0ss(calA0_str);
	stringstream phi0ss(phi0_str);
	stringstream T0ss(T0_str);
	stringstream NTss(NT_str);
	stringstream klss(kl_str);
	stringstream kbss(kb_str);
	stringstream seedss(seed_str);

	NCELLSss >> NCELLS;
	smallNVss >> smallNV;
	calA0ss >> calA0;
	phi0ss >> phi0;
	T0ss >> T0;
	NTss >> NT_dbl;
	klss >> kl;
	kbss >> kb;
	seedss >> seed;

	// cast input NT_dbl to integer
	NT = (int)NT_dbl;

	// open xyz file
	cout << "opening files" << endl;
	ofstream posout;
	posout.open(positionFile.c_str());
	if (!posout.is_open()){
		cout << "	** ERROR: position file " << positionFile << " could not be opened, ending." << endl;
		return 1;
	}

	// open contact matrix file
	ofstream enout;
	enout.open(energyFile.c_str());
	if (!enout.is_open()){
		cout << "	** ERROR: energy file " << energyFile << " could not be opened, ending." << endl;
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

	cout << "		phi0 		= " << phi0 << " 					" << endl;
	cout << "		T0 			= " << T0 << " 						" << endl;
	cout << "		kl 			= " << kl << "						" << endl;
	cout << "		kb 			= " << kb << "						" << endl;
	cout << "		seed 		= " << seed << "					" << endl << endl;

	cout << "		pos file 	= " << positionFile << "			" << endl;
	cout << "		en file 	= " << energyFile << " 				" << endl << endl;
	
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

	// initialize cell centers randomly
	for (ci=0; ci<cellDOF; ci += 2)
		dpos.at(ci) = L[ci % 2]*drand48();
	for (ci=cellDOF-1; ci>0; ci -= 2)
		dpos.at(ci) = L[ci % 2]*drand48();




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
	dtmin   	= 1e-1*dt0;

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
	double Kl, Kb, l0tmp, atmp, li, lim1, kappai, cx, cy;
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

					// compute area deviation
					atmp = area(vpos,ci,L,nv,szList);
					da = atmp - a0tmp;

					// shape force parameters
					Kl = nvtmp*l0tmp*kl;
					Kb = kb/(nvtmp*pow(l0tmp,2.0));

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
				vF[NDIM*gi] 		+= Kb*(3.0*(lix - lim1x) + lim2x - lip1x);
				vF[NDIM*gi + 1] 	+= Kb*(3.0*(liy - lim1y) + lim2y - lip1y);
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


	// print vertex positions to check placement
	cout << "\t** PRINTING POSITIONS TO FILE... " << endl;
	printPos(posout, vpos, a0, l0, L, nv, szList, phi0, NCELLS);







































	return 0;