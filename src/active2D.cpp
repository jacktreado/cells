/*

	Methods file for 2D active particles

*/


// include file
#include "deformableParticles2D.h"
#include "cellPacking2D.h"

// namespace
using namespace std;

// constant PI
const double PI = 4*atan(1);

// channel width (in particle diameters)
const double w0 = 1.5;			


// FUNCTION TO INITALIZE PARTICLE POSITIONS AS IF THEY WERE ACTIVE SOFT PARTICLES
// WITH DIAMETER SIGMA
// 	** In Pipe flow geometry, to simulate ABP flow for zebrafish study
// 
// 	NOTE: radii are input in units of sigma, the mean particle diameter
// 
// 	** ALSO, use fireMinimizeSP and spAttractiveForces, two functions NOT in this file
void cellPacking2D::initializeActiveStickySP(vector<double>& radii, int NV, double phiDisk, double sizeDispersion, double Lscale){
	// local variables
	int ci, vi, d, nvtmp;
	double r1, r2, g1, radsum;
	double xpos, ypos;
	double xmin, xmax, ymin, ymax;
	double calA;
	double rtmp, l0tmp, a0tmp;
	double delval = 1.0;

	// minimum number of vertices
	const int nvmin = 12;

	// check inputs to stick SP initialization
	if (radii.size() < NCELLS){
		cout << "	** ERROR: in initializing sticky SP, input radii vector size = " << radii.size() << ", which is != NCELLS (= " << NCELLS << "). ending." << endl;
		exit(1);
	}

	// output to console
	cout << "		-- In active stickySP initialization, initializing active SP particles" << endl;

	// initialize length scales as gaussian random variables (will becomes area square roots)
	radsum = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		// generate random numbers
		r1 = drand48();
		r2 = drand48();

		// calculate gaussian random variable using Box-Muller transform
		g1 = sqrt(-2.0*log(r1))*cos(2*PI*r2);

		// get radius
		radii.at(ci) = 0.5*(g1*sizeDispersion + 1.0);

		// add to lenscales sum for boundary size
		radsum += radii.at(ci)*radii.at(ci);
	}

	// determine box length from particle sizes and input packing fraction
	L.at(0) = sqrt(Lscale*PI*radsum/phiDisk);
	L.at(1) = sqrt(PI*radsum/(Lscale*phiDisk));

	// set phi to input
	phi = phiDisk;

	// reseed rng
	srand48(56835698*seed);

	// initialize cell information
	cout << "		-- Ininitializing cell objects" << endl;
	for (ci=0; ci<NCELLS; ci++){
		// boundary information
		for (d=0; d<NDIM; d++)
			cell(ci).setL(d,L.at(d));

		// x-direction is periodic, y-direction is not
		cell(ci).setpbc(0,1);
		cell(ci).setpbc(1,0);

		// number of vertices ( SIGMA SETS # OF VERTS )
		nvtmp = round(2.0*radii.at(ci)*NV);
		if (nvtmp > nvmin)
 			cell(ci).setNV(nvtmp);
		else
			cell(ci).setNV(nvmin);

		// array information
		cell(ci).initializeVertices();
		cell(ci).initializeCell();

		// initial length of polygon side
		l0tmp = 2.0*radii.at(ci)*sin(PI/nvtmp);

		// use rtmp slightly smaller than lenscale, so no overlaps at end
		rtmp = radii.at(ci) - 0.25*delval*l0tmp;

		// calculate a0 and l0 based on fact that they are regular polygons
		a0tmp = 0.5*nvtmp*pow(rtmp,2.0)*sin(2.0*PI/nvtmp);
		l0tmp = 2.0*rtmp*sin(PI/nvtmp);

		// set preferred area and length 
		cell(ci).seta0(a0tmp);
		cell(ci).setl0(l0tmp);
		cell(ci).setdel(delval);
	}

	// initialize particle positions
	cout << "		-- Ininitializing cell positions" << endl;
	for (ci=0; ci<NCELLS; ci++){
		// set min and max values of positions
		xmin = radii.at(ci);
		xmax = L.at(0) - radii.at(ci);
		ymin = radii.at(ci);
		ymax = L.at(1) - radii.at(ci);

		// get random location in pipe
		xpos = (xmax-xmin)*drand48() + xmin;
		ypos = (ymax-ymin)*drand48() + ymin;

		// set as initial position of com
		cell(ci).setCPos(0,xpos);
		cell(ci).setCPos(1,ypos);

		// initialize vertices as a regular polygon
		cell(ci).regularPolygon();
	}

	// initial time scales (t_0^2 = mass*sigma/f_0, mass = 0.25*PI*sigma^2)
	cout << "		-- Ininitializing time scale" << endl;
	dt = 0.05*sqrt(0.25*PI);
	dt0 = dt;

	// use FIRE in hopper geometry to relax overlaps
	cout << "		-- Using FIRE to relax overlaps..." << endl;
	fireMinimizeSP(radii,0.0);
}



// FUNCTION FOR SINGLE ACTIVE PARTICLE WITH BOTH CENTRAL AND VERTEX ACTIVITIES
void cellPacking2D::singleActiveCell(int NV, double phiInit, double calA0, double Dc, double vari, double v0){
	// local variables
	int vi, d, t;
	double rtmp, l0tmp, a0tmp, veltmp, postmp;
	double r1, r2, g1, g2, rvtmp, da;

	// set diameter of initial polygonal shapes to be 1
	rtmp = 0.5;

	// get boundary sizes (set to 1)
	for (d=0; d<NDIM; d++){
		L.at(d) = sqrt(PI*rtmp*rtmp/phiInit);
		cell(0).setL(d,L.at(d));
		cell(0).setpbc(d,1);
	}

	// number of vertices
	cell(0).setNV(NV);

	// array information
	cell(0).initializeVertices();
	cell(0).initializeCell();

	// initial length of polygon side
	l0tmp = 2.0*rtmp*sin(PI/NV);
	a0tmp = 0.5*NV*pow(rtmp,2.0)*sin(2.0*PI/NV);

	// set preferred area and length 
	cell(0).seta0(a0tmp);
	cell(0).setl0(l0tmp);
	cell(0).setdel(1.0);

	// set com in box center
	for (d=0; d<NDIM; d++)
		cell(0).setCPos(d,0.5*L.at(d));

	// initialize vertices as a regular polygon
	cell(0).regularPolygon();

	// set asphericity
	cell(0).setAsphericityConstA(calA0);

	// director for cell center
	double psi = 0.0;
	double dpsi = 0.0;

	// angular position of given site
	double psiVi = 0.0;

	// normal w.r.t. cell center
	double ntmp = 0.0;

	// self-propulsion mag
	double v0tmp = 0.0;

	// directors for each individual vertex
	// vector<double> alpha(NV,0.0);
	// vector<double> theta(NV,0.0);
	// int vip1, vim1;

	// get initial directors for each
	// for (vi=0; vi<NV; vi++){
	// 	r1 = drand48();
	// 	alpha.at(vi) = 2.0*PI*(2.0*r1 - 1.0);
	// }

	// set time step
	dt = 0.01*sqrt(PI*cell(0).area());
	dt0 = dt;

	// loop over time, integrate overdamped eqn of motion with active motility
	for (t=0; t<NT; t++){

		// output some information to console
		if (t % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 		Single cell motility, t = " << t << endl << endl;
			cout << "===================================================" << endl;

			// print if object has been opened already
			if (packingPrintObject.is_open()){
				cout << "	* Printing single cell vertex positions to file" << endl;
				printSystemPositions();
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing single cell vertex energy to file" << endl;
				printSystemEnergy();
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << totalKineticEnergy() << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << endl << endl;
		}

		// update positions based on forces (EULER)
		for (vi=0; vi<NV; vi++){
			// angular position of site
			psiVi = atan2(cell(0).vrel(vi,1),cell(0).vrel(vi,0));

			// get angular distance
			dpsi = psiVi - psi;
			dpsi -= 2.0*PI*round(dpsi/(2.0*PI));

			// get magnitude of self-propulsion
			v0tmp = v0*exp(-pow(dpsi,2.0)/(2.0*vari));

			// loop over dimensions, update positions and reset forces for next time
			for (d=0; d<NDIM; d++){
				// component of random vel director (0 = x, 1 = y)
				// rvtmp = (1-d)*cos(alpha.at(vi)) + d*sin(alpha.at(vi));

				// get velocities (= forces in overdamped regime)
				// veltmp = cell(0).vvel(vi,d) + v0*rvtmp;

				// get normal component (point in opposite direction if psi pointed away from site)
				ntmp = cell(0).vrel(vi,d)/sqrt(cell(0).vrel(vi,0)*cell(0).vrel(vi,0) + cell(0).vrel(vi,1)*cell(0).vrel(vi,1));
				if (abs(dpsi) > 0.5*PI)
					ntmp *= -1;


				// add activity to velocity
				veltmp = cell(0).vvel(vi,d) + v0tmp*ntmp;

				// if new position in outflow region, place back in hopper
				postmp = cell(0).vpos(vi,d) + dt*veltmp;

				// update positions (EULER STEP)
				cell(0).setVPos(vi,d,postmp);

				// update velocities
				cell(0).setVVel(vi,d,veltmp);

				// reset forces and energies
				cell(0).setVForce(vi,d,0.0);
				cell(0).setUInt(vi,0.0);
			}

			// determine angle for the current velocity
			// theta.at(vi) = atan2(cell(0).vvel(vi,1),cell(0).vvel(vi,0));
		}


		/*
		// UPDATE ACTIVE DIRECTORS
		for (vi=0; vi<NV; vi++){
			// draw uniform random variables
			r1 = drand48();
			r2 = drand48();

			// use Box-Muller trnsfrm to get GRVs for active variable
			g1 = sqrt(-2.0*log(r1))*cos(2*PI*r2);

			// determine vip1, vim1
			vip1 = (vi+1) % NV;
			vim1 = (vi-1+NV) % NV;

			// update alpha based on euler scheme and velocity alignment
			da = sin(theta.at(vip1) - alpha.at(vi)) + ep*sin(psi - alpha.at(vi)) + 2.0*Dv*g1;
			alpha.at(vi) += dt*da;
			if (alpha.at(vi) > 2.0*PI)
				alpha.at(vi) -= 2.0*PI;
			else if (alpha.at(vi) < -2.0*PI)
				alpha.at(vi) += 2.0*PI;
		}
		*/

		// draw uniform random variables
		r1 = drand48();
		r2 = drand48();

		// use Box-Muller trnsfrm to get GRVs for active variable
		g1 = sqrt(-2.0*log(r1))*cos(2*PI*r2);

		// update psi
		psi += dt*2.0*Dc*g1;

		// DEBUG: PRINT OUT EVERY TIME WHEN NAN POPS UP
		/*
		if (t > 18000){
			cout << "Printing everything, nan coming up soon!" << endl;

			// print config and energy
			if (packingPrintObject.is_open()){
				cout << "	* Printing vetex positions to file" << endl;
				printSystemPositions();
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing cell energy to file" << endl;
				printSystemEnergy();
			}

			if (t > 18500){
				cout << "should have found a nan by now, ending" << endl;
				exit(1);
			}

			cout << endl;
			cout << "t = " << t << endl;
			cout << "phi = " << phi << endl;
			for (vi=0; vi<NV; vi++){
				cout << alpha.at(vi) << "; ";
				cout << cell(0).vpos(vi,0) << "; " << cell(0).vpos(vi,1) << "; ";
				cout << cell(0).vvel(vi,0) << "; " << cell(0).vvel(vi,1) << "; ";
				cout << cell(0).vforce(vi,0) << "; " << cell(0).vforce(vi,1);
				cout << endl;
			}
		}
		*/

		// update cpos
		cell(0).updateCPos();

		// calculate forces between vertices
		cell(0).shapeForces();

		// update velocities based on forces
		for (vi=0; vi<NV; vi++){
			for (d=0; d<NDIM; d++)
				cell(0).setVVel(vi,d,cell(0).vforce(vi,d));
		}
	}
}




// FUNCTION TO CALCULATE FORCES BETWEEN STICKY, ACTIVE SP PARTICLES
void cellPacking2D::spActiveForces(vector<double>& radii){
	// local variables
	int ci, cj, vi, d;
	double l1, l2;
	double meanAttract, meanDiam;
	double energyScale;
	double overlap = 0.0;
	double uv = 0.0;
	double ftmp, utmp;
	vector<double> distanceVec(NDIM,0.0);
	double contactDistance = 0.0; 
	double centerDistance = 0.0;

	// fraction of a for l1
	const double l1Scale = 0.5;

	// vector to store number of attractive contacts per cell
	vector<int> nac(NCELLS,0);

	// reset virial stresses to 0
	sigmaXX = 0.0;
	sigmaXY = 0.0;
	sigmaYX = 0.0;
	sigmaYY = 0.0;

	// reset forces
	for (ci=0; ci<NCELLS; ci++){
		for (d=0; d<NDIM; d++)
			cell(ci).setCForce(d,0.0);
	}

	// reset contacts
	resetContacts();

	// get number of attractive contacts
	// loop over cell pairs
	for (ci=0; ci<NCELLS; ci++){
		for (cj=0; cj<NCELLS; cj++){

			// determine if ci and cj engage in attractive bond
			meanAttract = 0.5*(cell(ci).geta()*radii.at(ci) + cell(cj).geta()*radii.at(cj));
			meanDiam = radii.at(ci) + radii.at(cj);

			// get cell-cell distance
			centerDistance = 0.0;
			for (d=0; d<NDIM; d++)
				centerDistance += pow(cell(ci).cellDistance(cell(cj),d),2.0);
			centerDistance = sqrt(centerDistance);

			if (centerDistance < meanAttract + meanDiam)
				addContact(ci,cj);
		}
	}

	// determine number of attractive contacts on each particle
	for (ci=0; ci<NCELLS; ci++){
		for (cj=ci+1; cj<NCELLS; cj++){
			nac.at(ci) += contacts(ci,cj);
			nac.at(cj) += contacts(ci,cj);
		}
	}

	// loop over cells and cell pairs, calculate shape and interaction forces
	for (ci=0; ci<NCELLS; ci++){

		// loop over pairs, add info to contact matrix
		for (cj=ci+1; cj<NCELLS; cj++){
			
			// contact distance
			contactDistance = radii.at(ci) + radii.at(cj);

			// attractive shell
			meanAttract = 0.5*(cell(ci).geta() + cell(cj).geta());

			// center-to-center distance
			centerDistance = 0.0;
			for (d=0; d<NDIM; d++){
				// vectorial quantity
				distanceVec.at(d) = cell(ci).cellDistance(cell(cj),d);

				// add to distance
				centerDistance += pow(distanceVec.at(d),2);
			}
			centerDistance = sqrt(centerDistance);

			// attraction quantities
			l2 = meanAttract;
			l1 = 0.5*l1Scale*l2*(1.0/nac.at(cj) + 1.0/nac.at(ci));

			// check if within interaction zone
			if (centerDistance < meanAttract + contactDistance){
				// overlap scale
				overlap = centerDistance/contactDistance;

				// energy scale
				energyScale = contactDistance;

				// IF in zone to use repulsive force (and, if attractiveParam > 0, bottom of attractive well)
				if (centerDistance < contactDistance*(1 + l1)){
					// interaction potential
					utmp = 0.5 * energyScale * (pow(1 - overlap,2) - l1*l2);

					// scalar part of force
					ftmp = 1.0 - overlap;
				}
				// IF attractiveParam > 0, in attractive well
				else if (centerDistance >= contactDistance*(1 + l1) && centerDistance < contactDistance*(1 + l2) && meanAttract > 0.0){
					// interaction potential
					utmp =  -(0.5*energyScale*l1/(l2 - l1)) * pow(overlap - 1 - l2,2);

					// scalar part of force
					ftmp = (l1/(l2 - l1)) * (overlap - 1.0 - l2);
				}

				// add to u and f based on utmp and ftmp

				// potential energy from utmp
				for (vi=0; vi<cell(ci).getNV(); vi++)
					cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());
				for (vi=0; vi<cell(cj).getNV(); vi++)
					cell(cj).setUInt(vi,cell(cj).uInt(vi) + utmp/cell(cj).getNV());

				// forces from ftmp
				for (d=0; d<NDIM; d++){
					// unit vector
					uv = distanceVec.at(d)/centerDistance;

					// add to forces (MIND FORCE DIRECTION; rij points from i -> j, so need extra minus sign)
					cell(ci).setCForce(d,cell(ci).cforce(d) - ftmp*uv);
					cell(cj).setCForce(d,cell(cj).cforce(d) + ftmp*uv);

					// add to virial stresses
					if (d == 0){
						sigmaXX += ftmp*uv*distanceVec.at(0);
						sigmaXY += ftmp*uv*distanceVec.at(1);
					}
					else{
						sigmaYX += ftmp*uv*distanceVec.at(0);
						sigmaYY += ftmp*uv*distanceVec.at(1);
					}
				}
			} 
		}
	}

	// normalize virial stresses by the area
	sigmaXX /= L.at(0)*L.at(1);
	sigmaXY /= L.at(0)*L.at(1);
	sigmaYX /= L.at(0)*L.at(1);
	sigmaYY /= L.at(0)*L.at(1);
}






// ZEBRAFISH HORSESHOE GEOMETRY SIMS
void cellPacking2D::initializeActiveZebrafish(vector<double>& radii, int NV, double phiDisk, double sizeDispersion, double R0){
	// local variables
	int ci, vi, d, nvtmp, areaint, Ncurr, Ntmp;
	double Ach, Ahs;
	double r1, r2, g1, radsum;
	double w, ltmp;
	double xpos, ypos;
	double xmin, xmax, ymin, ymax, ypackmax, ydiff;
	double calA;
	double rtmp, l0tmp, a0tmp;
	double delval = 1.0;

	// minimum number of vertices
	const int nvmin = 12;

	// check inputs to stick SP initialization
	if (radii.size() < NCELLS){
		cout << "	** ERROR: in initializing sticky SP, input radii vector size = " << radii.size() << ", which is != NCELLS (= " << NCELLS << "). ending." << endl;
		exit(1);
	}

	// output to console
	cout << "		-- In active stickySP initialization, initializing active SP particles" << endl;

	// get initial area of horseshoe and channels
	Ahs = 0.5*PI*R0*R0 + 2.0*w0*R0;
	Ach = R0*2.0*(R0-w0);

	// decide number of particles to place in sim box
	Ncurr = round(phiDisk*(Ahs + Ach)/(0.25*PI));
	if (Ncurr > NCELLS)
		Ncurr = NCELLS;

	// initialize length scales as gaussian random variables (will becomes area square roots)
	radsum = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		// generate random numbers
		r1 = drand48();
		r2 = drand48();

		// calculate gaussian random variable using Box-Muller transform
		g1 = sqrt(-2.0*log(r1))*cos(2*PI*r2);

		// get radius
		radii.at(ci) = 0.5*(g1*sizeDispersion + 1.0);

		// add to lenscales sum for boundary size
		if (ci < Ncurr)
			radsum += radii.at(ci)*radii.at(ci);
	}

	// determine horseshow radius from particle sizes and input packing fraction
	L.at(0) = R0;
	L.at(1) = 0.0;

	// determine main channel width
	w = R0 - w0;

	// set phi of initial particles
	phi = PI*radsum/(Ach + Ahs);
	cout << "		-- Initial phi for first " << Ncurr << "/" << NCELLS << " particles = " << phi << endl;

	// reseed rng
	srand48(56835698*seed);

	// initialize cell information
	cout << "		-- Ininitializing cell objects" << endl;
	for (ci=0; ci<NCELLS; ci++){
		// boundary information
		for (d=0; d<NDIM; d++)
			cell(ci).setL(d,L.at(d));

		// x- and y-directions are not periodic
		cell(ci).setpbc(0,0);
		cell(ci).setpbc(1,0);

		// number of vertices ( SIGMA SETS # OF VERTS )
		nvtmp = round(2.0*radii.at(ci)*NV);
		if (nvtmp > nvmin)
 			cell(ci).setNV(nvtmp);
		else
			cell(ci).setNV(nvmin);

		// array information
		cell(ci).initializeVertices();
		cell(ci).initializeCell();

		// initial length of polygon side
		l0tmp = 2.0*radii.at(ci)*sin(PI/nvtmp);

		// use rtmp slightly smaller than lenscale, so no overlaps at end
		rtmp = radii.at(ci) - 0.25*delval*l0tmp;

		// calculate a0 and l0 based on fact that they are regular polygons
		a0tmp = 0.5*nvtmp*pow(rtmp,2.0)*sin(2.0*PI/nvtmp);
		l0tmp = 2.0*rtmp*sin(PI/nvtmp);

		// set preferred area and length 
		cell(ci).seta0(a0tmp);
		cell(ci).setl0(l0tmp);
		cell(ci).setdel(delval);
	}

	// initialize particle positions
	cout << "		-- Ininitializing cell positions" << endl;
	for (ci=0; ci<NCELLS; ci++){
		// partices in main area
		if (ci < Ncurr){
			// draw random integer between 0 and 4, to decide where
			// to place particle
			areaint = round(4*drand48());

			// main channel
			if (areaint <= 1){
				xmin = -0.5*w;
				xmax = 0.5*w;
				ymin = -R0;
				ymax = 0.0;
			}
			// left side channel
			else if (areaint == 2){
				xmin = -R0;
				xmax = -0.5*(w + 2.0*w0);
				ymin = -R0;
				ymax = 0.0;
			}
			// right side channel
			else if (areaint == 3){
				xmin = 0.5*(w + 2.0*w0);
				xmax = R0;
				ymin = -R0;
				ymax = 0.0;
			}
			// horseshoe chamber
			else if (areaint == 4){
				xmin = -0.5*R0;
				xmax = 0.5*R0;
				ymin = 0;
				ymax = 2.0*w0 + 0.75*R0;
			}

			// get random location in pipe
			xpos = (xmax-xmin)*drand48() + xmin;
			ypos = (ymax-ymin)*drand48() + ymin;

			// set as initial position of com
			cell(ci).setCPos(0,xpos);
			cell(ci).setCPos(1,ypos);
		}
		// cells in reservoir
		else{
			// place well away from tailbid
			xpos = 0.0;
			ypos = -4.0*R0;

			// set as initial position of com
			cell(ci).setCPos(0,xpos);
			cell(ci).setCPos(1,ypos);
		}

		// initialize vertices as a regular polygon
		cell(ci).regularPolygon();
	}

	// initial time scales (t_0^2 = mass*sigma/f_0, mass = 0.25*PI*sigma^2)
	cout << "		-- Ininitializing time scale" << endl;
	dt = 0.05*sqrt(0.25*PI);
	dt0 = dt;

	// use FIRE in zebrafish horseshow geometry to relax overlaps
	cout << "		-- Using FIRE to relax overlaps for first Ncurr = " << Ncurr << "particles ..." << endl;

	// minimize energy only for first Ncurr particles
	Ntmp = NCELLS;
	NCELLS = Ncurr;
	fireMinimizeZebrafishSP(radii,0.0);
	NCELLS = Ntmp;
}


void cellPacking2D::fireMinimizeZebrafishSP(vector<double>& radii, double attractiveParam){
	// HARD CODE IN FIRE PARAMETERS
	const double alpha0 	= 0.25;
	const double finc 		= 1.01;
	const double fdec 		= 0.5;
	const double falpha 	= 0.99;
	const double dtmax 		= 10*dt0;
	const double dtmin 		= 0.02*dt0;
	const int NMIN 			= 200;
	const int NNEGMAX 		= 2000;
	const int NDELAY 		= 1000;
	int npPos				= 0;
	int npNeg 				= 0;
	int npPMIN				= 0;
	int closed 				= 1;
	double alpha 			= alpha0;
	double alphat 			= alpha;
	double t 				= 0.0;
	double Ptol 			= 1e-8;
	double Ktol 			= 1e-12;
	bool converged 			= false;

	// local variables
	int ci,vi,d,itr,itrMax;
	double P,vstarnrm,fstarnrm,vtmp,ftmp,ptmp;
	double Knew, Pvirial;
	double Kcheck, Pcheck;

	// reset time step
	dt = dt0;

	// reset velocities to 0
	for (ci=0; ci<NCELLS; ci++){
		for (d=0; d<NDIM; d++){
			cell(ci).setCVel(d,0.0);
			cell(ci).setCForce(d,0.0);
		}
	}

	// psi for printing
	vector<double> psi(NCELLS,0.0);

	// initial wall pressure to 0
	double wallPressure = 0.0;

	// initialize forces (neglect damping forces, only interactions)
	resetContacts();
	spActiveZebrafishWallForces(radii,wallPressure);

	// initialize virial pressure from pressure from last time
	Pvirial = 0.5*(sigmaXX + sigmaYY);

	// update kinetic energy based on com velocity
	Knew = 0.0;
	for (ci=0; ci<NCELLS; ci++)
		Knew += 0.5*(PI*pow(radii.at(ci),2))*sqrt(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2));

	// calc check variables
	Kcheck = Knew/(NDIM*NCELLS);
	Pcheck = Pvirial/(NDIM*NCELLS);

	// iterate through MD time until system converged
	itrMax = 1e6;
	for (itr=0; itr<itrMax; itr++){

		// output some information to console
		if (itr % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 	FIRE MINIMIZATION, itr = " << itr << endl << endl;
			cout << "===================================================" << endl;
			cout << "	* Run data:" << endl;
			cout << "	* NCELLS 	= " << NCELLS << endl;
			cout << "	* Kcheck 	= " << Kcheck << endl;
			cout << "	* Pcheck 	= " << Pcheck << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << "	* alpha 	= " << alpha << endl;
			cout << "	* P 		= " << P << endl;
			cout << endl;
			cout << "	* zebrafish data:" << endl;
			cout << "	* wallPressure = " << wallPressure << endl;
			cout << endl << endl;
		}

		// Step 1. calculate P and norms
		P = 0.0;
		vstarnrm = 0.0;
		fstarnrm = 0.0;
		for (ci=0; ci<NCELLS; ci++){
			for (d=0; d<NDIM; d++){
				// get tmp variables
				ftmp = cell(ci).cforce(d);
				vtmp = cell(ci).cvel(d);

				// calculate based on all vertices on all cells
				P += ftmp*vtmp;
				vstarnrm += vtmp*vtmp;
				fstarnrm += ftmp*ftmp;
			}
		}


		// get norms
		vstarnrm = sqrt(vstarnrm);
		fstarnrm = sqrt(fstarnrm);


		// Step 2. Adjust simulation based on net motion of system
		if (P > 0){
			// increment pos counter
			npPos++;

			// reset neg counter
			npNeg = 0;

			// update alpha_t for next time
			alphat = alpha;

			// alter sim if enough positive steps taken
			if (npPos > NMIN){
				// change time step
				if (dt*finc < dtmax)
					dt *= finc;
				else
					dt = dtmax;

				// decrease alpha
				alpha *= falpha;
			}
		}
		else{
			// reset pos counter
			npPos = 0;

			// rest neg counter
			npNeg++;

			// check for stuck sim
			if (npNeg > NNEGMAX){
				// // print if objects have been opened already
				// if (packingPrintObject.is_open()){
				// 	cout << "	* Printing zebrafish SP center positions to file" << endl;
				// 	printPositionsZebrafishSP(radii,psi);
				// }
				
				// if (energyPrintObject.is_open()){
				// 	cout << "	* Printing zebrafish SP energy to file" << endl;
				// 	printEnergyZebrafishSP(Knew);
				// }

				cout << "	** FIRE has stalled..." << endl;
				cout << "	** Kcheck = " << Kcheck << endl;
				cout << "	** Pcheck = " << Pcheck << endl;
				cout << "	** itr = " << itr << ", t = " << t << endl;
				cout << "	** Breaking out of FIRE protocol." << endl;
				break;
			}

			// decrease time step if past initial delay
			if (itr > NMIN){
				// decrease time step 
				if (dt*fdec > dtmin)
					dt *= fdec;
				else
					dt = dtmin;

				// change alpha
				alpha = alpha0;
				alphat = alpha;
			}

			// take half step backwards
			for (ci=0; ci<NCELLS; ci++){
				for (d=0; d<NDIM; d++)
					cell(ci).setCPos(d,cell(ci).cpos(d) - 0.5*dt*cell(ci).cvel(d));
			}

			// reset velocities to 0
			for (ci=0; ci<NCELLS; ci++){
				for (d=0; d<NDIM; d++)
					cell(ci).setCVel(d,0.0);
			}
		}

		// update velocities if forces are acting
		if (fstarnrm > 0){
			for (ci=0; ci<NCELLS; ci++){
				for (d=0; d<NDIM; d++){
					vtmp = (1 - alphat)*cell(ci).cvel(d) + alphat*(cell(ci).cforce(d)/fstarnrm)*vstarnrm;
					cell(ci).setCVel(d,vtmp);
				}
			}
		}

		// verlet position update
		spPosVerlet();

		// reset contacts before force calculation
		resetContacts();

		// calculate forces between disks
		spAttractiveForces(radii,attractiveParam);
		spActiveZebrafishWallForces(radii,wallPressure);

		// verlet velocity update
		spVelVerlet(radii);

		// update t
		t += dt;

		// update virial pressure
		Pvirial = 0.5*(sigmaXX + sigmaYY);

		// update kinetic energy based on com velocity
		Knew = 0.0;
		for (ci=0; ci<NCELLS; ci++)
			Knew += 0.5*(PI*pow(radii.at(ci),2))*sqrt(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2));

		// update if Pvirial under tol
		if (abs(Pvirial) < Ptol)
			npPMIN++;
		else
			npPMIN = 0;

		// calc check variables
		Kcheck = Knew/(NDIM*NCELLS);
		Pcheck = Pvirial/(NDIM*NCELLS);

		// check for convergence
		converged = (abs(Pcheck) < Ptol && npPMIN > NMIN);
		converged = (converged || (abs(Pcheck) > Ptol && Kcheck < Ktol));

		if (converged){
			// // print if objects have been opened already
			// if (packingPrintObject.is_open()){
			// 	cout << "	* Printing zebrafish SP center positions to file" << endl;
			// 	printPositionsZebrafishSP(radii,psi);
			// }
			
			// if (energyPrintObject.is_open()){
			// 	cout << "	* Printing zebrafish SP energy to file" << endl;
			// 	printEnergyZebrafishSP(Knew);
			// }

			cout << "	** FIRE has converged!" << endl;
			cout << "	** Kcheck = " << Kcheck << endl;
			cout << "	** Pcheck = " << Pcheck << endl;
			cout << "	** itr = " << itr << ", t = " << t << endl;
			cout << "	** Breaking out of FIRE protocol." << endl;
			break;
		}
	}

	// reset dt to be original value before ending function
	dt = dt0;

	// if no convergence, just stop
	if (itr == itrMax)
		cout << "	** FIRE not converged in itrMax = " << itr << " force evaluations" << endl;
}


// calculate wall forces due to horseshoe
void cellPacking2D::spActiveZebrafishWallForces(vector<double>& radii, double& wallPressure){
	// local variables
	int ci, vi, d; 							// indices
	int intwall, wall;						// wall check variables
	double R0, w, h;						// zebrafish boundary variables
	double x, y;							// particle positions (IN UNITS OF SIGMA)
	double sigma; 							// particle diameter
	double lw, lwx, lwy;					// elements of vector pointing from wall to particle
	double overlap;							// overlap of particle with wall
	double ftmp, utmp, uv;					// force/energy of particle overlap with walls
	double capDist;							// distance of particles to spherical caps
	double edgeDist;						// distance of particles to horseshoe edge
	double centerDist; 						// distance of particles to horseshoe center
	vector<double> distVec(NDIM,0.0);		// vector for distances between particles and wall

	// reset wall pressure to 0
	wallPressure = 0.0;

	// head radius and tail height
	R0 = L.at(0);
	h = L.at(1);

	// main channel width
	w = R0 - w0;

	// loop over cells
	for (ci=0; ci<NCELLS; ci++){
		// get sigma (2*radius)
		sigma = 2.0*radii.at(ci);

		// get particle positions
		x = cell(ci).cpos(0);
		y = cell(ci).cpos(1);

		// check horseshoe side walls
		if (y < h){
			// if particle in main channel
			if (x < 0.5*(w + w0) && x > -0.5*(w + w0)){
				// check variable for wall force
				wall = 0;

				// particle touching left wall, force vector points in +x dir
				if (x < 0.5*(sigma - w)){
					lwx = x + 0.5*w;
					wall = -1;
				}
				// particle touching right wall, force vector points in -x dir
				else if (x > 0.5*(w - sigma)){
					lwx = 0.5*w - x;
					wall = 1;
				}
				else 
					wall = 0;

				if (wall != 0){
					// overlap with wall
					overlap = 2.0*lwx/sigma;

					// add to x force ONLY
					ftmp = 1 - overlap;
					cell(ci).setCForce(0,cell(ci).cforce(0) - wall*ftmp);

					// add to energies
					utmp = 0.25*sigma*pow(1 - overlap,2);
					for (vi=0; vi<cell(ci).getNV(); vi++)
						cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());
				}				
			}
			// if particle is in a side channel
			else{
				// first check interior channel walls
				if ( (x > -0.5*(w + 2.0*w0 + sigma) && x < 0) || (x < 0.5*(w + 2.0*w0 + sigma) && x > 0) ){
					// check variable for interior wall force
					intwall = 0;

					// particle touching left interior wall
					if (x > -0.5*(w + 2.0*w0 + sigma) && x < 0){
						lwx = -0.5*(w + 2.0*w0) - x;
						intwall = -1;
					}
					// particle touching right interior wall
					else if (x < 0.5*(w + 2.0*w0 + sigma) && x > 0){
						lwx = x - 0.5*(w + 2.0*w0);
						intwall = 1;
					}

					// overlap with wall
					overlap = 2.0*lwx/sigma;

					// add to x force ONLY
					ftmp = 1 - overlap;
					cell(ci).setCForce(0,cell(ci).cforce(0) + intwall*ftmp);

					// add to energies
					utmp = 0.25*sigma*pow(1 - overlap,2);
					for (vi=0; vi<cell(ci).getNV(); vi++)
						cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());
				}
				// next check exterior walls
				if ( (x < -R0 + 0.5*sigma && x < 0) || (x > R0 - 0.5*sigma && x > 0) ) {
					// check variable for exterior wall force
					wall = 0;

					// particle touching left exterior wall
					if (x < -R0 + 0.5*sigma && x < 0){
						lwx = x + R0;
						wall = -1;
					}
					else if (x > R0 - 0.5*sigma && x > 0){
						lwx = R0 - x;
						wall = 1;
					}

					// overlap with wall
					overlap = 2.0*lwx/sigma;

					// add to x force ONLY
					ftmp = 1 - overlap;
					cell(ci).setCForce(0,cell(ci).cforce(0) - wall*ftmp);

					// add to energies
					utmp = 0.25*sigma*pow(1 - overlap,2);
					for (vi=0; vi<cell(ci).getNV(); vi++)
						cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());

				}
			}

			// bottom walls
			if (y < -R0 + 0.5*sigma){
				// get wall overlap
				lwy = y + R0;

				// overlap with wall
				overlap = 2.0*lwy/sigma;

				// add to y force ONLY
				ftmp = 1 - overlap;
				cell(ci).setCForce(1,cell(ci).cforce(1) + ftmp);

				// add to energies
				utmp = 0.25*sigma*pow(1 - overlap,2);
				for (vi=0; vi<cell(ci).getNV(); vi++)
					cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());
			}
		}
		// otherwise, in the horseshoe section
		else{
			// first check disk caps where channels end
			if (x < 0){
				// check left cap
				capDist = pow(x + 0.5*(w + w0),2.0) + pow(y - h,2.0);
				capDist = sqrt(capDist);

				// get distance from particle center to cap edge
				centerDist = capDist - 0.5*w0;

				// save normal of spherical cap to vector
				distVec.at(0) = (x + 0.5*(w + w0))/capDist;
				distVec.at(1) = (y - h)/capDist;
			}
			else{
				// check right cap
				capDist = pow(x - 0.5*(w + w0),2.0) + pow(y - h,2.0);
				capDist = sqrt(capDist);

				// get distance from particle center to cap edge
				centerDist = capDist - 0.5*w0;

				// save distances to vector
				distVec.at(0) = (x - 0.5*(w + w0))/capDist;
				distVec.at(1) = (y - h)/capDist;
			}

			// if overlap with cap, calculate ghost particle force
			if (centerDist < 0.5*sigma){
				// calculate overlap
				overlap = 2.0*centerDist/sigma;

				// scalar part of force
				ftmp = 1 - overlap;

				// vectorial parts of force
				for (d=0; d<NDIM; d++)
					cell(ci).setCForce(d,cell(ci).cforce(d) + ftmp*distVec.at(d));

				// add to energies
				utmp = 0.25*sigma*pow(1 - overlap,2);
				for (vi=0; vi<cell(ci).getNV(); vi++)
					cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());
			}

			// then check overlap with horseshoe edge: get radial distance from horseshoe center
			centerDist = sqrt(x*x + (y-h-2.0*w0)*(y-h-2.0*w0));

			// get distance to horseshoe edge
			edgeDist = R0 - centerDist;

			// get unit vector from horshoe to particle
			distVec.at(0) = -x/centerDist;
			distVec.at(1) = (h+2.0*w0-y)/centerDist;

			// if interacting with horseshoe edge, add forces
			if (edgeDist < 0.5*sigma){
				// calculate overlap
				overlap = 2.0*edgeDist/sigma;

				// scalar part of force
				ftmp = 1 - overlap;

				// vectorial parts of force
				for (d=0; d<NDIM; d++)
					cell(ci).setCForce(d,cell(ci).cforce(d) + ftmp*distVec.at(d));

				// add to wall pressure
				wallPressure += ftmp/(PI*R0);

				// add to energies
				utmp = 0.25*sigma*pow(1 - overlap,2);
				for (vi=0; vi<cell(ci).getNV(); vi++)
					cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());
			}
		}
	}
}


void cellPacking2D::spActiveZebrafishNVE(vector<double>& radii, double v0){
	// local variables
	int t;
	int ci = 0;
	double K = 0;
	double Pvirial = 0.0;
	double wallPressure = 0.0;

	// psi for printing
	vector<double> psi(NCELLS,0.0);

	// initialize particles in main channel
	for (ci=0; ci<NCELLS; ci++){
		cell(ci).setCVel(1,v0);
	}

	// loop over time, run NVE dynamics
	for (t=0; t<NT; t++){
		// update virial pressure
		Pvirial = 0.5*(sigmaXX + sigmaYY);

		// update kinetic energy based on com velocity
		K = 0.0;
		for (ci=0; ci<NCELLS; ci++)
			K += 0.5*(PI*pow(radii.at(ci),2))*(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2));

		// output some information to console
		if (t % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 	ZEBRAFISH NVE, t = " << t << endl << endl;
			cout << "===================================================" << endl;
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << K << endl;
			cout << "	* Pvirial 	= " << Pvirial << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << endl;
			cout << "	* zebrafish data:" << endl;
			cout << "	* wallPressure = " << wallPressure << endl;
			cout << endl << endl;

			// print if objects have been opened already
			if (packingPrintObject.is_open()){
				cout << "	* Printing zebrafish SP center positions to file" << endl;
				printPositionsZebrafishSP(radii,psi);
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing zebrafish SP energy to file" << endl;
				printEnergyZebrafishSP(K);
			}
		}

		// verlet position update
		spPosVerlet();

		// reset contacts before force calculation
		resetContacts();

		// calculate forces between disks
		spAttractiveForces(radii,0.1);
		spActiveZebrafishWallForces(radii,wallPressure);

		// verlet velocity update
		spVelVerlet(radii);
	}
}

// Active particles: can be active brownian particles if vtau is large enough
void cellPacking2D::spActiveZebrafishVicsek(vector<double>& radii, double attractionParam, double v0, double Dr, double vtau, double Pthresh, double dh){
	// local variables
	int t, d, ci, vi, Ncurr, Ntmp;
	double K = 0;
	double Pvirial = 0.0;
	double wallPressure = 0.0;
	double veltmp, postmp;
	double r1, r2, grv, dpsi, rvtmp;
	double ymin = 0.0;

	// wall values
	double R0 = L.at(0);
	double h = L.at(1);
	double w = R0 - w0;
	double hcum = 0.0;

	// determine number of particles in main channel
	Ncurr = 0;
	for (ci=0; ci<NCELLS; ci++){
		if (cell(ci).cpos(1) > -3.0*R0)
			Ncurr++;
	}

	// angular directors (all point to the right initially)
	vector<double> psi(NCELLS,0.0);

	// initialize particles in main channel
	for (ci=0; ci<Ncurr; ci++){
		cell(ci).setCVel(1,v0);
		psi.at(ci) = 2.0*PI*drand48() - PI;
	}

	// loop over time, run brownian dynamics with 
	for (t=0; t<NT; t++){
		// update virial pressure
		Pvirial = 0.5*(sigmaXX + sigmaYY);

		// update kinetic energy based on com velocity
		K = 0.0;
		for (ci=0; ci<Ncurr; ci++)
			K += 0.5*(PI*pow(radii.at(ci),2))*(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2));

		// output some information to console
		if (t % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 	ZEBRAFISH vicsek particles, t = " << t << endl << endl;
			cout << "===================================================" << endl;
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << K << endl;
			cout << "	* Pvirial 	= " << Pvirial << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << "	* Ncurr 	= " << Ncurr << endl;
			cout << endl;
			cout << "	* zebrafish data:" << endl;
			cout << "	* wallPressure = " << wallPressure << endl;
			cout << "	* h = " << h << endl;
			cout << " 	* dh = " << dh << endl;
			cout << "	* Pthresh = " << Pthresh << endl;
			cout << endl << endl;

			// print if objects have been opened already
			if (packingPrintObject.is_open()){
				cout << "	* Printing zebrafish SP center positions to file" << endl;
				printPositionsZebrafishSP(radii,psi);
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing zebrafish SP energy to file" << endl;
				printEnergyZebrafishSP(K);
			}
		}

		// update positions of particles in main channel based on Euler update
		for (ci=0; ci<Ncurr; ci++){

			// loop over dimensions, update positions and reset forces for next time
			for (d=0; d<NDIM; d++){
				// component of random vel director (0 = x, 1 = y)
				rvtmp = (1-d)*cos(psi.at(ci)) + d*sin(psi.at(ci));

				// get velocities (= forces in overdamped regime)
				veltmp = cell(ci).cvel(d) + v0*rvtmp;

				// if new position in outflow region, place back in hopper
				postmp = cell(ci).cpos(d) + dt*veltmp;

				// update positions (EULER STEP)
				cell(ci).setCPos(d,postmp);

				// update velocities
				cell(ci).setCVel(d,veltmp);

				// reset forces
				cell(ci).setCForce(d,0.0);
			}

			// set interaction energy to 0
			for (vi=0; vi<cell(ci).getNV(); vi++)
				cell(ci).setUInt(vi,0.0);

			//	UPDATE ACTIVE DIRECTOR
			r1 = drand48();
			if (vtau < 1e2)
				psi.at(ci) += dt*((1.0/vtau)*asin(cos(psi.at(ci))*cell(ci).cvel(1) - sin(psi.at(ci))*cell(ci).cvel(0)) + 4.0*Dr*PI*(r1 - 0.5));
			else
				psi.at(ci) += dt*4.0*Dr*PI*(r1 - 0.5);
		}

		// update wall position based on forces
		if (wallPressure > Pthresh){
			// update size
			h += dh;
			L.at(1) = h;
			hcum += dh;

			// if cumulative gain is large enough, then add a reserve particle
			if (2.0*hcum*R0 > 0.25*PI && Ncurr < NCELLS){
				// increment number
				Ncurr++;

				// reset cumulative h
				hcum = 0.0;

				// set new position, velocity of particle
				cell(Ncurr-1).setCPos(0,(w - 2.0*radii.at(Ncurr-1))*r1 - 0.5*w + radii.at(Ncurr-1));
				cell(Ncurr-1).setCPos(1,-R0 + radii.at(Ncurr-1));
				psi.at(Ncurr-1) = 0.5*PI;
			}
		}

		// update boundary
		// h += dh;
		// L.at(1) = h;
		// hcum += dh;

		// // if cumulative gain is large enough, then add a reserve particle
		// if (2.0*hcum*R0 > 0.25*PI && Ncurr < NCELLS){
		// 	// increment number
		// 	Ncurr++;

		// 	// reset cumulative h
		// 	hcum = 0.0;

		// 	// set new position, velocity of particle
		// 	cell(Ncurr-1).setCPos(0,(w - 2.0*radii.at(Ncurr-1))*r1 - 0.5*w + radii.at(Ncurr-1));
		// 	cell(Ncurr-1).setCPos(1,-R0 + radii.at(Ncurr-1));
		// 	psi.at(Ncurr-1) = 0.5*PI;
		// }

		// reset contacts before force calculation
		resetContacts();

		// calculate forces between disks (only check forces up to Ncurr)
		Ntmp = NCELLS;
		NCELLS = Ncurr;
		spAttractiveForces(radii,attractionParam);
		spActiveZebrafishWallForces(radii,wallPressure);
		NCELLS = Ntmp;

		// update velocities based on forces
		for (ci=0; ci<Ncurr; ci++){
			for (d=0; d<NDIM; d++)
				cell(ci).setCVel(d,cell(ci).cforce(d));
		}
	}
}






// PIPE FLOW GEOMETRY

void cellPacking2D::spActivePipeWallForces(vector<double>& radii){
	// local variables
	int ci, vi;
	double overlap, sigma, x, y, lwy;
	double ftmp, utmp;

	// loop over cells
	for (ci=0; ci<NCELLS; ci++){
		// get sigma (2*radius)
		sigma = 2*radii.at(ci);

		// get particle positions
		x = cell(ci).cpos(0);
		y = cell(ci).cpos(1);

		// if true, interacting with bottom wall
		if (y < 0.5*sigma){
			// vector from wall to particle
			lwy = y;

			// overlap with wall
			overlap = 2.0*lwy/sigma;

			// add to y force ONLY (points in positive y direction)
			ftmp = 1 - overlap;
			cell(ci).setCForce(1,cell(ci).cforce(1) + ftmp);

			// add to energies
			utmp = 0.25*sigma*pow(1 - overlap,2);
			for (vi=0; vi<cell(ci).getNV(); vi++)
				cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());

			// virial stress in YY direction
			sigmaYY += ftmp*lwy;
		}

		// if true, interacting with top wall
		if (y > L.at(1) - 0.5*sigma){
			// vector from particle to wall
			lwy = L.at(1) - y;

			// overlap with wall
			overlap = 2.0*lwy/sigma;

			// add to y force ONLY (points in negative y direction)
			ftmp = 1 - overlap;
			cell(ci).setCForce(1,cell(ci).cforce(1) - ftmp);

			// add to energies
			utmp = 0.25*sigma*pow(1 - overlap,2);
			for (vi=0; vi<cell(ci).getNV(); vi++)
				cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());

			// virial stress in YY direction 
			sigmaYY -= ftmp*lwy;
		}
	}
}

void cellPacking2D::spActivePipeNVE(vector<double>& radii, double T0){
	// local variables
	int t, ci, d;
	double Pvirial, K;
	int closed = 1;

	// check that NT has been set 
	if (NT <= 0){
		cout << "	** ERROR: in active pipe SP NVE, sim length NT = " << NT << ", which is <= 0. ending." << endl;
		exit(1);
	}

	// initialize velocities using Gaussian random variables
	vector<double> pmean(NDIM,0.0);
	double r1, r2, grv, vscale, mtmp;

	// loop over velocities, give them initial conditions
	for (ci=0; ci<NCELLS; ci++){
		for (d=0; d<NDIM; d++){
			// draw uniform random variables
			r1 = drand48();
			r2 = drand48();

			// use Box-Muller trnsfrm to get GRV
			grv = sqrt(-2.0*log(r1))*cos(2*PI*r2);

			// add to cell velocity
			cell(ci).setCVel(d,grv);		
		}
	}

	// get system momentum
	for (d=0; d<NDIM; d++)
		pmean.at(d) /= NCELLS;

	// subtract off mean
	K = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		for (d=0; d<NDIM; d++){
			// particle mass
			mtmp = PI*pow(radii.at(ci),2);

			// subtract of com motion
			cell(ci).setCVel(d,cell(ci).cvel(d) - pmean.at(d)/mtmp);

			// calc ek
			K += 0.5*mtmp*pow(cell(ci).cvel(d),2);
		}
	}

	// scale velocities so K to start is T0
	vscale = sqrt(T0/K);
	for (ci=0; ci<NCELLS; ci++){
    	for (d=0; d<NDIM; d++)
        	cell(ci).setCVel(d,cell(ci).cvel(d)*vscale);
    }

	// loop over time, run NVE dynamics
	for (t=0; t<NT; t++){
		// update virial pressure
		Pvirial = 0.5*(sigmaXX + sigmaYY);

		// update kinetic energy based on com velocity
		K = 0.0;
		for (ci=0; ci<NCELLS; ci++)
			K += 0.5*(PI*pow(radii.at(ci),2))*sqrt(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2));

		// output some information to console
		if (t % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 		Soft, active Particle NVE, t = " << t << endl << endl;
			cout << "===================================================" << endl;

			// print if object has been opened already
			if (packingPrintObject.is_open()){
				cout << "	* Printing SP center positions to file" << endl;
				printPositionsStickySP(radii);
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing SP energy to file" << endl;
				printEnergyStickySP();
			}

			if (statPrintObject.is_open()){
				cout << "	* Printing SP energy to file" << endl;
				printSystemContacts();
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << K << endl;
			cout << "	* virial P 	= " << Pvirial << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << endl << endl;
		}

		// verlet position update
		spPosVerlet();

		// reset contacts before force calculation
		resetContacts();

		// calculate forces between disks
		spAttractiveForces(radii,0.1);
		spActivePipeWallForces(radii);

		// verlet velocity update
		spVelVerlet(radii);
	}
}

void cellPacking2D::spActivePipeFlow(vector<double>& radii, double attractiveParam, double v0, double Dr){
	// local variables
	int t, ci, vi, d;
	double Pvirial, K, veltmp, postmp;
	double r1, r2, grv, dpsi, rvtmp;

	// angular directors (all point to the right initially)
	vector<double> psi(NCELLS,0.0);

	// loop over time, integrate overdamped eqn of motion with active motility
	for (t=0; t<NT; t++){
		// update virial pressure
		Pvirial = 0.5*(sigmaXX + sigmaYY);

		// update kinetic energy based on com velocity
		K = 0.0;
		for (ci=0; ci<NCELLS; ci++)
			K += 0.5*(PI*pow(radii.at(ci),2))*sqrt(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2));

		// output some information to console
		if (t % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 		Soft Particle FLOW, t = " << t << endl << endl;
			cout << "===================================================" << endl;

			// print if object has been opened already
			if (packingPrintObject.is_open()){
				cout << "	* Printing SP center positions to file" << endl;
				printPositionsStickySP(radii);
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing SP energy to file" << endl;
				printEnergyStickySP();
			}

			if (statPrintObject.is_open()){
				cout << "	* Printing SP energy to file" << endl;
				printSystemContacts();
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << K << endl;
			cout << "	* virial P 	= " << Pvirial << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << endl << endl;
		}

		// update positions based on forces (EULER)
		for (ci=0; ci<NCELLS; ci++){
			// loop over dimensions, update positions and reset forces for next time
			for (d=0; d<NDIM; d++){
				// component of random vel director (0 = x, 1 = y)
				rvtmp = (1-d)*cos(psi.at(ci)) + d*sin(psi.at(ci));

				// get velocities (= forces in overdamped regime)
				veltmp = cell(ci).cvel(d) + v0*rvtmp;

				// if new position in outflow region, place back in hopper
				postmp = cell(ci).cpos(d) + dt*veltmp;

				// update positions (EULER STEP)
				cell(ci).setCPos(d,postmp);

				// update velocities
				cell(ci).setCVel(d,veltmp);

				// reset forces
				cell(ci).setCForce(d,0.0);
			}

			// set interaction energy to 0
			for (vi=0; vi<cell(ci).getNV(); vi++)
				cell(ci).setUInt(vi,0.0);

			/*
				UPDATE ACTIVE DIRECTOR
			*/

			// draw uniform random variables
			r1 = drand48();
			r2 = drand48();

			// use Box-Muller trnsfrm to get GRV for active variable
			grv = sqrt(-2.0*log(r1))*cos(2*PI*r2);

			// update psi based on euler scheme
			psi.at(ci) += dt*2.0*Dr*grv;
		}

		// reset contacts before force calculation
		resetContacts();

		// calculate forces between disks
		spAttractiveForces(radii,attractiveParam);
		spActivePipeWallForces(radii);

		// update velocities based on forces
		for (ci=0; ci<NCELLS; ci++){
			for (d=0; d<NDIM; d++)
				cell(ci).setCVel(d,cell(ci).cforce(d));
		}
	}
}


// PRINT FUNCTIONS

void cellPacking2D::printPositionsZebrafishSP(vector<double>& radii, vector<double>& psi){
	// local variables
	int ci,d;
	int w1 = 25;
	int w2 = 15;

	// check to see if file is open
	if (!packingPrintObject.is_open()) {
		cout << "	ERROR: packingPrintObject is not open in spActiveZebrafishPosPrint(), ending." << endl;
		exit(1);
	}

	// print information starting information
	packingPrintObject << setw(w1) << left << "NEWFR" << " " << endl;
	packingPrintObject << setw(w1) << left << "NUMCL" << setw(w2) << right << NCELLS << endl;

	// print hopper information
	packingPrintObject << setw(w1) << left << "BNDRY";
	packingPrintObject << setw(w2) << right << L.at(0);
	packingPrintObject << setw(w2) << right << L.at(1);
	packingPrintObject << setw(w2) << setprecision(6) << right << w0;
	packingPrintObject << endl;

	// print stress information
	packingPrintObject << setw(w1) << left << "VRIAL";
	packingPrintObject << setw(w2) << right << sigmaXX;
	packingPrintObject << setw(w2) << right << sigmaXY;
	packingPrintObject << setw(w2) << right << sigmaYX;
	packingPrintObject << setw(w2) << right << sigmaYY;
	packingPrintObject << endl;

	// print header for information
	packingPrintObject << setw(w1) << left << "SINFO";
	packingPrintObject << setw(w2) << right << "id";
	packingPrintObject << setw(w2) << right << "r";
	packingPrintObject << setw(w2) << right << "x";
	packingPrintObject << setw(w2) << right << "y";
	packingPrintObject << setw(w2) << right << "vx";
	packingPrintObject << setw(w2) << right << "vy";
	packingPrintObject << setw(w2) << right << "fx";
	packingPrintObject << setw(w2) << right << "fy";
	packingPrintObject << endl;

	// loop over cells, print positions, forces, velocities
	for (ci=0; ci<NCELLS; ci++){
		// print row label
		packingPrintObject << setw(w1) << left << "SCELL";

		// print sp index
		packingPrintObject << setw(w2) << right << ci;

		// print radius
		packingPrintObject << setw(w2) << right << radii.at(ci);

		// print polarization
		packingPrintObject << setw(w2) << right << psi.at(ci);

		// print sp positions
		for (d=0; d<NDIM; d++)
			packingPrintObject << setw(w2) << right << cell(ci).cpos(d);

		// print sp velocities
		for (d=0; d<NDIM; d++)
			packingPrintObject << setw(w2) << right << cell(ci).cvel(d);

		// print new line
		packingPrintObject << endl;
	}

	// print end frame
	packingPrintObject << setw(w1) << left << "ENDFR" << " " << endl;
}

void cellPacking2D::printEnergyZebrafishSP(double K){
	// local variables
	int ci, d;
	double w = L.at(0) - w0;

	// check to see if file is open
	if (!energyPrintObject.is_open()) {
		cout << "	ERROR: energyPrintObject is not open in printSystemEnergy(), ending." << endl;
		exit(1);
	}

	// calculate zebrafish-specific states

	// fraction of cells on left side of boundary
	double fracOnLeft = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		if (cell(ci).cpos(0) < 0)
			fracOnLeft += 1.0;
	}
	fracOnLeft /= NCELLS;

	// polarization
	double polarization = 0.0;
	vector<double> meanUnitVec(NDIM,0.0);
	for (ci=0; ci<NCELLS; ci++){

		// get ci vel vector size
		polarization = 0.0;
		for (d=0; d<NDIM; d++)
			polarization += pow(cell(ci).cvel(d),2);
		polarization = sqrt(polarization);

		// update mean unit vect
		for (d=0; d<NDIM; d++)
			meanUnitVec.at(d) += cell(ci).cvel(d)/(polarization*NCELLS);
	}

	// calculate polarization from magnitude of mean velocity direction
	polarization = 0.0;
	for (d=0; d<NDIM; d++)
		polarization += pow(meanUnitVec.at(d),2);
	polarization = sqrt(polarization);

	// calculate mean ADM velocity
	double meanADMVel = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		if (cell(ci).cpos(1) < L.at(1) && cell(ci).cpos(0) < 0.5*(w + w0) && cell(ci).cpos(0) > -0.5*(w + w0))
			meanADMVel += sqrt(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2))/NCELLS;
	}

	// loop over particles, print cell energy
	energyPrintObject << setw(30) << setprecision(16) << right << interactionPotentialEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << K;
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaXX;
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaXY;
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaYX;
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaYY;
	energyPrintObject << setw(30) << setprecision(16) << right << fracOnLeft;
	energyPrintObject << setw(30) << setprecision(16) << right << polarization;
	energyPrintObject << setw(30) << setprecision(16) << right << meanADMVel;
	energyPrintObject << endl;
}










