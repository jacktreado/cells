/*

	Methods file for 2D hopper flows

*/


// include file
#include "deformableParticles2D.h"
#include "cellPacking2D.h"

// namespace
using namespace std;

// constants
const double PI = 4*atan(1);

// 2D HOPPER FLOW FUNCTIONS




// FUNCTION TO INITALIZE PARTICLE POSITIONS AS IF THEY WERE SOFT PARTICLES
// WITH DIAMETER SIGMA
// 
// assume following will already be initialized:
// 	** NCELLS, NPRINT, L (using w,w0,th)
// 	** SET PBCS TO 0
// 
// 	NOTE: radii should be input in units of sigma, the mean particle diameter
// 		** also will give the ratio of nv_i to NV, the number of vertices on the mean
void cellPacking2D::initializeHopperSP(vector<double>& radii, double w0, double w, double th, double Lmin, int NV){
	// local variables
	int ci, vi, nvtmp;
	double calA;
	double xpos, ypos;
	double xmin, xmax, ymin, ymax;

	// minimum number of vertices
	const int nvmin = 12;

	// check inputs to hopper initialization
	if (w0 < 0.0){
		cout << "	** ERROR: in initializing hopper, input w0 = " << w0 << ", which is < 0. ending." << endl;
		exit(1);
	}
	else if (w < 0.0){
		cout << "	** ERROR: in initializing hopper, input w = " << w << ", which is < 0. ending." << endl;
		exit(1);
	}
	else if (w > w0){
		cout << "	** ERROR: in initializing hopper, input w = " << w << ", which is > w0. ending." << endl;
		exit(1);
	}
	else if (th < 0){
		cout << "	** ERROR: in initializing hopper, input th = " << th << ", which is < 0. ending." << endl;
		exit(1);
	}

	// output to console
	cout << "		-- In hopper initialization, initializing cells and relaxing initial overlaps as SP particles" << endl;

	// initialize cell information
	cout << "		-- Ininitializing cell objects" << endl;
	for (ci=0; ci<NCELLS; ci++){
		// boundary information ( SET PBCS TO 0 )
		cell(ci).setL(L);
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

		// calculate a0 based on fact that they are regular polygons
		calA = nvtmp*tan(PI/nvtmp)/PI;
		cell(ci).setAsphericity(calA);
	}

	// initialize particle positions
	cout << "		-- Ininitializing cell positions" << endl;
	for (ci=0; ci<NCELLS; ci++){
		// set min and max values of positions
		xmin = -Lmin*L;
		xmax = -radii.at(ci);
		ymin = radii.at(ci);
		ymax = w0 - radii.at(ci);

		// get random location in hopper RESERVOIR
		xpos = (xmin-xmax)*drand48() + xmax;
		ypos = (ymax-ymin)*drand48() + ymin;

		// set as initial position of com
		cell(ci).setCPos(0,xpos);
		cell(ci).setCPos(1,ypos);

		// initialize vertices as a regular polygon
		cell(ci).regularPolygon();
	}

	// initialize phi
	cout << "		-- Ininitializing packing fraction...";
	phi = hopperPackingFraction(radii,w0,w,th);
	cout << "which is phi = " << phi << endl;

	// initial time scales (t_0^2 = mass*sigma/f_0, mass = 0.25*PI*sigma^2)
	cout << "		-- Ininitializing time scale" << endl;
	dt = 0.1*sqrt(0.25*PI);
	dt0 = dt;

	// use FIRE in hopper geometry to relax overlaps
	cout << "		-- Using FIRE to relax overlaps..." << endl;
	fireMinimizeHopperSP(radii,w0,w,th);
}

void cellPacking2D::fireMinimizeHopperSP(vector<double>& radii, double w0, double w, double th){
	// HARD CODE IN FIRE PARAMETERS
	const double alpha0 	= 0.25;
	const double finc 		= 1.1;
	const double fdec 		= 0.5;
	const double falpha 	= 0.99;
	const double dtmax 		= 10*dt0;
	const double dtmin 		= 0.02*dt0;
	const int NMIN 			= 20;
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
	double Ktol 			= 1e-24;
	bool converged 			= false;

	// local variables
	int ci,vi,d,itr,itrMax;
	double P,vstarnrm,fstarnrm,vtmp,ftmp,ptmp;
	double Knew, Pvirial;

	// reset time step
	dt = dt0;

	// initialize forces (neglect damping forces, only interactions)
	resetContacts();
	hopperForcesSP(radii,w0,w,th,0.0,closed);

	// initialize virial pressure from pressure from last time
	Pvirial = 0.5*(sigmaXX + sigmaYY);

	// update kinetic energy based on com velocity
	Knew = 0.0;
	for (ci=0; ci<NCELLS; ci++)
		Knew += 0.5*(PI*pow(radii.at(ci),2))*sqrt(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2));

	// reset velocities to 0
	for (ci=0; ci<NCELLS; ci++){
		for (d=0; d<NDIM; d++)
			cell(ci).setCVel(d,0.0);
	}

	// iterate through MD time until system converged
	itrMax = 5e5;
	for (itr=0; itr<itrMax; itr++){

		// output some information to console
		if (itr % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 	FIRE MINIMIZATION, itr = " << itr << endl << endl;
			cout << "===================================================" << endl;

			// print if object has been opened already
			if (packingPrintObject.is_open()){
				cout << "	* Printing SP center positions to file" << endl;
				printHopperSP(radii,w0,w,th,0.0);
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing SP energy to file" << endl;
				printSystemEnergy(itr,Pvirial,Knew);
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << Knew/Ktol << endl;
			cout << "	* Pvirial 	= " << Pvirial/Ptol << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << "	* alpha 	= " << alpha << endl;
			cout << "	* P 		= " << P << endl;
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
			if (npNeg > NNEGMAX)
				break;

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
		hopperPosVerletSP();

		// reset contacts before force calculation
		resetContacts();

		// calculate forces between disks (with door closed)
		hopperForcesSP(radii,w0,w,th,0.0,closed);

		// verlet velocity update
		hopperVelVerletSP(radii);

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

		// check for convergence
		converged = (abs(Pvirial) < Ptol && npPMIN > NMIN);
		converged = (converged || (abs(Pvirial) > 2*Ptol && Knew < Ktol));

		if (converged){
			cout << "	** FIRE has converged!" << endl;
			cout << "	** Knew = " << Knew << endl;
			cout << "	** itr = " << itr << ", t = " << t << endl;
			cout << "	** virial P = " << Pvirial << endl;
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


// calculate all forces on SP model particles
// 	** if closed = 1, orifice is closed off by wall
// 	** b is damping coefficient
// 	** g is force strength (i.e. generalized gravity)
void cellPacking2D::hopperForcesSP(vector<double>& radii, double w0, double w, double th, double g, int closed){
	// local variables
	int ci, cj, vi, d;
	double contactDistance = 0.0; 
	double centerDistance = 0.0; 
	double overlap = 0.0;
	double uv = 0.0;
	double ftmp, utmp;
	vector<double> distanceVec(NDIM,0.0);

	// reset virial stresses to 0
	sigmaXX = 0.0;
	sigmaXY = 0.0;
	sigmaYX = 0.0;
	sigmaYY = 0.0;

	// get disk-disk forces
	for (ci=0; ci<NCELLS; ci++){
		for (cj=ci+1; cj<NCELLS; cj++){
			// contact distance
			contactDistance = radii.at(ci) + radii.at(cj);

			// center-to-center distance
			centerDistance = 0.0;
			for (d=0; d<NDIM; d++){
				// vectorial quantity
				distanceVec.at(d) = cell(ci).cellDistance(cell(cj),d);

				// add to distance
				centerDistance += pow(distanceVec.at(d),2);
			}

			// check for contact
			if (contactDistance*contactDistance > centerDistance){
				// add to contact checking
				addContact(ci,cj);

				// get true distance
				centerDistance = sqrt(centerDistance);

				// overlap scale
				overlap = centerDistance/contactDistance;

				// force scale
				ftmp = (1 - overlap);

				// add to potential energy (energy should increase because particles are growing)
				utmp = 0.5*contactDistance*pow(1 - overlap,2);
				for (vi=0; vi<cell(ci).getNV(); vi++)
					cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());
				for (vi=0; vi<cell(cj).getNV(); vi++)
					cell(cj).setUInt(vi,cell(cj).uInt(vi) + utmp/cell(cj).getNV());

				// add to forces
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

	// wall forces
	hopperWallForcesSP(radii,w0,w,th,closed);

	// body force (in x direction)
	if (g > 1e-16){
		for (ci=0; ci<NCELLS; ci++)
			cell(ci).setCForce(0,cell(ci).cforce(0) + g*pow(radii.at(ci),2));
	}
}


// wall forces between cells as soft disks (SP model)
// 	** if closed = 1, orifice is closed off by wall
void cellPacking2D::hopperWallForcesSP(vector<double>& radii, double w0, double w, double th, int closed){
	// local variables
	int ci, vi; 							// indices
	double x, y;							// particle positions (IN UNITS OF SIGMA)
	double sigma; 							// particle diameter
	double t, c, s;							// tangent, cosine, sine
	double hPlus, hMinus;					// height of angled wall
	double dyPlus, dyMinus; 				// distance from top/bottom wall
	double yPlusMin, yMinusMax; 			// cutoff positions for wall forces
	double lw, lwx, lwy;					// elements of vector pointing from wall to particle
	double overlap;							// overlap of particle with wall
	double ftmp, utmp;						// force/energy of particle overlap with walls
	bool nearEdge = false;					// boolean to check if near the edge of the hopper

	// preliminary calculations
	t = tan(th);
	c = cos(th);
	s = sin(th);

	// loop over cells
	for (ci=0; ci<NCELLS; ci++){
		// get sigma (2*radius)
		sigma = 2*radii.at(ci);

		// get particle positions
		x = cell(ci).cpos(0);
		y = cell(ci).cpos(1);

		// check hopper walls
		if (x > -sigma*s){

			// in hopper, but not near cusp at edge
			nearEdge = (x > L - 0.5*sigma*s && x < L);
			nearEdge = (nearEdge && (y > 0.5*(w0 + w) - 0.5*sigma || y < 0.5*(w0 - w) + 0.5*sigma));

			// if particle in hopper, either in bulk or near outflow point
			if (x < L - 0.5*sigma*s || nearEdge){
				// check ymin for walls
				yPlusMin 	= w0 - x*t - 0.5*sigma/c;
				yMinusMax 	= x*t + 0.5*sigma/c;

				// if true, interacting with bottom wall
				if (y < yMinusMax){
					// vector to wall
					lwx = s*(x*s - y*c);
					lwy = c*(y*c - x*s);

					// distance
					lw = sqrt(lwx*lwx + lwy*lwy);

					if (lw < 0.5*sigma){
						// overlap with wall
						overlap = 2.0*lw/sigma;

						// force
						ftmp = 1 - overlap;
						cell(ci).setCForce(0,cell(ci).cforce(0) + ftmp*(lwx/lw));
						cell(ci).setCForce(1,cell(ci).cforce(1) + ftmp*(lwy/lw));

						// add to energies
						utmp = 0.25*sigma*pow(1 - overlap,2);
						for (vi=0; vi<cell(ci).getNV(); vi++)
							cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());

						// virial stresses
						sigmaXX += ftmp*(lwx/lw)*lwx;
						sigmaXY += ftmp*(lwx/lw)*lwy;
						sigmaYX += ftmp*(lwy/lw)*lwx;
						sigmaYY += ftmp*(lwy/lw)*lwy;
					}
				}


				// if true, interacting with top wall
				if (y > yPlusMin){
					// vector to wall
					lwx = s*(x*s + (y - w0)*c);
					lwy = c*((y - w0)*c + x*s);

					// distance
					lw = sqrt(lwx*lwx + lwy*lwy);

					if (lw < 0.5*sigma){
						// overlap with wall
						overlap = 2.0*lw/sigma;

						// force
						ftmp = 1 - overlap;
						cell(ci).setCForce(0,cell(ci).cforce(0) + ftmp*(lwx/lw));
						cell(ci).setCForce(1,cell(ci).cforce(1) + ftmp*(lwy/lw));

						// add to energies
						utmp = 0.25*sigma*pow(1 - overlap,2);
						for (vi=0; vi<cell(ci).getNV(); vi++)
							cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());

						// virial stresses
						sigmaXX += ftmp*(lwx/lw)*lwx;
						sigmaXY += ftmp*(lwx/lw)*lwy;
						sigmaYX += ftmp*(lwy/lw)*lwx;
						sigmaYY += ftmp*(lwy/lw)*lwy;
					}
				}
			}
		}

		// check reservoir walls
		if (x < 0){
			// check ymin for walls
			yPlusMin 	= w0 - 0.5*sigma;
			yMinusMax 	= 0.5*sigma;

			// if true, interacting with bottom wall
			if (y < yMinusMax){
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
			if (y > yPlusMin){
				// vector from particle to wall
				lwy = w0 - y;

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

		// if orifice closed, check orifice wall
		if (closed == 1){
			if (x > L - radii.at(ci)){

				// vector from particle to wall
				lwx = L - x;

				// overlap with wall
				overlap = 2.0*lwx/sigma;

				// add to x force ONLY (points in negative x direction)
				ftmp = 1 - overlap;
				cell(ci).setCForce(0,cell(ci).cforce(0) - ftmp);

				// add to energies
				utmp = 0.25*sigma*pow(1 - overlap,2);
				for (vi=0; vi<cell(ci).getNV(); vi++)
					cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());

				// virial stress in XX direction 
				sigmaXX -= ftmp*lwx;
			}
		}
		
	}
}


// wall forces between cells as droplets (DP model)
// 	** if closed = 1, orifice is closed off by wall
;void cellPacking2D::hopperWallForcesDP(double w0, double w, double th, int closed){
	
}







// function to run NVE dynamics using velocity-verlet to check energy conservation
void cellPacking2D::hopperSPNVE(vector<double>& radii, double w0, double w, double th, double T0){
	// local variables
	int t, ci, d;
	double Pvirial, K;
	int closed = 1;

	// check that NT has been set 
	if (NT <= 0){
		cout << "	** ERROR: in hopper SP NVE, sim length NT = " << NT << ", which is <= 0. ending." << endl;
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
			cout << " 		Soft Particle NVE, t = " << t << endl << endl;
			cout << "===================================================" << endl;

			// print if object has been opened already
			if (packingPrintObject.is_open()){
				cout << "	* Printing SP center positions to file" << endl;
				printHopperSP(radii,w0,w,th,0.0);
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing SP energy to file" << endl;
				printSystemEnergy(t,Pvirial,K);
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << K << endl;
			cout << "	* virial P 	= " << Pvirial << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << endl << endl;
		}

		// verlet position update
		hopperPosVerletSP();

		// reset contacts before force calculation
		resetContacts();

		// calculate forces between disks (with door closed)
		hopperForcesSP(radii,w0,w,th,0.0,closed);

		// verlet velocity update
		hopperVelVerletSP(radii);
	}
}



// function to flow cells through hopper as soft particles (SP) using body force
void cellPacking2D::flowHopperSP(vector<double>& radii, double w0, double w, double th, double g){
	// local variables
	int t, ci, vi, d, itr, itrMax;
	double Pvirial, K, veltmp, postmp, xmax, xmin;
	int closed = 0;

	// find max x position
	xmax = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		if (cell(ci).cpos(0) > xmax)
			xmax = cell(ci).cpos(0);
	}

	// flow until max x is near orifice
	itr = 0;
	itrMax = 1e6;
	while (xmax < L - 1.0 && itr < itrMax){
		// update virial pressure
		Pvirial = 0.5*(sigmaXX + sigmaYY);

		// update kinetic energy based on com velocity
		K = 0.0;
		for (ci=0; ci<NCELLS; ci++)
			K += 0.5*(PI*pow(radii.at(ci),2))*sqrt(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2));

		// new phi
		phi = hopperPackingFraction(radii,w0,w,th);

		// output some information to console
		if (itr % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 		Soft Particle FLOW TOWARD OPENING, g = " << g << ", itr = " << itr << endl << endl;
			cout << "===================================================" << endl;

			// print if object has been opened already
			if (packingPrintObject.is_open()){
				cout << "	* Printing SP center positions to file" << endl;
				printHopperSP(radii,w0,w,th,0.0);
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing SP energy to file" << endl;
				printSystemEnergy(itr,Pvirial,K);
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << K << endl;
			cout << "	* virial P 	= " << Pvirial << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << "	* g 		= " << g << endl;
			cout << "	* xmax/L 	= " << xmax/L << endl;
			cout << endl << endl;
		}

		// update positions based on forces (EULER)
		for (ci=0; ci<NCELLS; ci++){
			// loop over dimensions, update positions and reset forces for next time
			for (d=0; d<NDIM; d++){
				// get velocities (= forces in overdamped regime)
				veltmp = cell(ci).cforce(d);

				// update positions (EULER STEP)
				cell(ci).setCPos(d,cell(ci).cpos(d) + dt*veltmp);

				// update velocities
				cell(ci).setCVel(d,veltmp);

				// reset forces
				cell(ci).setCForce(d,0.0);
			}

			// set interaction energy to 0
			for (vi=0; vi<cell(ci).getNV(); vi++)
				cell(ci).setUInt(vi,0.0);
		}

		// reset contacts before force calculation
		resetContacts();

		// calculate forces
		hopperForcesSP(radii,w0,w,th,g,closed);

		// find max x position
		xmax = 0.0;
		for (ci=0; ci<NCELLS; ci++){
			if (cell(ci).cpos(0) > xmax)
				xmax = cell(ci).cpos(0);
		}

		// update iterator
		itr++;
	}

	// check for iterator error
	if (itr == itrMax){
		cout << "	** itr = itrMax = " << itrMax << ", particles did not flow to boundary in enough time, ending..." << endl;
		exit(1);
	}

	// loop over time, replace outflow to back of the hopper
	cout << "	** SYSTEM HAS FLOWED TO OPENING, STARTING FLOW SIMULATION " << endl;
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
			cout << " 		Soft Particle FLOW, g = " << g << ", t = " << t << endl << endl;
			cout << "===================================================" << endl;

			// print if object has been opened already
			if (packingPrintObject.is_open()){
				cout << "	* Printing SP center positions to file" << endl;
				printHopperSP(radii,w0,w,th,0.0);
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing SP energy to file" << endl;
				printSystemEnergy(t,Pvirial,K);
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << K << endl;
			cout << "	* virial P 	= " << Pvirial << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << "	* g 		= " << g << endl;
			cout << endl << endl;
		}

		// find left-most grain
		xmin = 2*L;
		for (ci=0; ci<NCELLS; ci++){
			if (cell(ci).cpos(0) < xmin)
				xmin = cell(ci).cpos(0);
		}

		// update positions based on forces (EULER)
		for (ci=0; ci<NCELLS; ci++){
			// loop over dimensions, update positions and reset forces for next time
			for (d=0; d<NDIM; d++){
				// get velocities (= forces in overdamped regime)
				veltmp = cell(ci).cvel(d);

				// if new position in outflow region, place back in hopper
				postmp = cell(ci).cpos(d) + dt*veltmp;

				if (d == 0 && postmp > L + 2.0*radii.at(ci)){
					// put particle somewhere near back of reservoir
					postmp = xmin + radii.at(ci);

					// pick random y value
					cell(ci).setCPos(1,drand48()*w0);

					// set x value
					veltmp = 2.0*g*pow(radii.at(ci),2);

					// set new y velocity to 0
					cell(ci).setCVel(1,0.0);
				}

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
		}

		// reset contacts before force calculation
		resetContacts();

		// calculate forces
		hopperForcesSP(radii,w0,w,th,g,closed);

		// update velocities based on forces
		for (ci=0; ci<NCELLS; ci++){
			for (d=0; d<NDIM; d++)
				cell(ci).setCVel(d,cell(ci).cforce(d));
		}
	}
}


// function to flow cells through hopper as deformable particles (DP) using body force of scale g
void cellPacking2D::flowHopperDP(double w0, double w, double th, double g){

}


// calculate packing fraction in hopper geometry (GIVEN RADII)
double cellPacking2D::hopperPackingFraction(vector<double>& radii, double w0, double w, double th){
	// local variables
	int ci;
	double aH, aR, aT;

	// reservoir area
	aR = w0*L;

	// hopper area
	aH = w*w0 + pow(w0-w,2)/(4.0*tan(th));

	// total area
	aT = aR + aH;

	// initialize phi to 0
	phi = 0.0;

	// loop over radii of cells in flow part of hopper
	for (ci=0; ci<NCELLS; ci++){
		if (cell(ci).cpos(0) > -L)
			phi += PI*pow(radii.at(ci),2);
	}

	// scale areas by total area
	phi /= aT;

	// return packing fraction
	return phi;
}


// set packing fraction in hopper geometry (GIVEN RADII)
void cellPacking2D::setHopperPackingFraction(vector<double>& radii, double phiNew, double w0, double w, double th){
	// local variables
	int ci;
	double phiOld = hopperPackingFraction(radii,w0,w,th);

	// get scale factor
	double scaleFactor = pow(phiNew/phiOld,1.0/NDIM);

	// scale particle radii
	for (ci=0; ci<NCELLS; ci++)
		radii.at(ci) *= scaleFactor;
}


// update positions of SP using velocity-verlet
void cellPacking2D::hopperPosVerletSP(){
	// local variables
	int ci, vi, d;
	double postmp, acctmp;

	// update com position
	for (ci=0; ci<NCELLS; ci++){
		// loop over positions
		for (d=0; d<NDIM; d++){
			// calculate com acceleration
			acctmp = 0.0;
			for (vi=0; vi<cell(ci).getNV(); vi++)
				acctmp += cell(ci).vacc(vi,d);

			// update new position based on acceleration
			postmp = cell(ci).cpos(d) + dt*cell(ci).cvel(d) + 0.5*dt*dt*acctmp;

			// update new positions
			cell(ci).setCPos(d,postmp);

			// set forces to 0
			cell(ci).setCForce(d,0.0);
		}

		// set interaction energy to 0
		for (vi=0; vi<cell(ci).getNV(); vi++)
			cell(ci).setUInt(vi,0.0);
	}
}

void cellPacking2D::hopperVelVerletSP(vector<double>& radii){
	// local variables
	int ci, vi, d;
	double veltmp, aold, anew, diskMass;

	// update com velocity
	for (ci=0; ci<NCELLS; ci++){
		// get disk mass
		diskMass = PI*pow(radii.at(ci),2);

		// loop over velocities
		for (d=0; d<NDIM; d++){
			// get current velocity
			veltmp = cell(ci).cvel(d);

			// calculate old com acceleration
			aold = 0.0;
			for (vi=0; vi<cell(ci).getNV(); vi++)
				aold += cell(ci).vacc(vi,d);

			// get new accelation
			anew = cell(ci).cforce(d)/diskMass;

			// update velocity
			veltmp += 0.5*dt*(anew + aold);

			// set new velocity and acceleration
			cell(ci).setCVel(d,veltmp);
			for (vi=0; vi<cell(ci).getNV(); vi++)
				cell(ci).setVAcc(vi,d,anew/cell(ci).getNV());
		}
	}
}


// function to print SP information to file
void cellPacking2D::printHopperSP(vector<double>& radii, double w0, double w, double th, double g){
	// local variables
	int ci,d;
	int w1 = 25;
	int w2 = 15;

	// check to see if file is open
	if (!packingPrintObject.is_open()) {
		cout << "	ERROR: packingPrintObject is not open in printHopperSP(), ending." << endl;
		exit(1);
	}

	// print information starting information
	packingPrintObject << setw(w1) << left << "NEWFR" << " " << endl;
	packingPrintObject << setw(w1) << left << "NUMCL" << setw(w2) << right << NCELLS << endl;

	// print hopper information
	packingPrintObject << setw(w1) << left << "HOPPR";
	packingPrintObject << setw(w2) << right << w0;
	packingPrintObject << setw(w2) << right << w;
	packingPrintObject << setw(w2) << setprecision(6) << right << th;
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

		// print sp positions
		for (d=0; d<NDIM; d++)
			packingPrintObject << setw(w2) << right << cell(ci).cpos(d);

		// print sp velocities
		for (d=0; d<NDIM; d++)
			packingPrintObject << setw(w2) << right << cell(ci).cvel(d);

		// print sp forces
		for (d=0; d<NDIM; d++)
			packingPrintObject << setw(w2) << right << cell(ci).cforce(d);

		// print new line
		packingPrintObject << endl;
	}

	// print end frame
	packingPrintObject << setw(w1) << left << "ENDFR" << " " << endl;
}

// function to print DP information to file
void cellPacking2D::printHopperDP(double w0, double w, double th, double g){}










