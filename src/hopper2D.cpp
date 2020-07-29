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
// 	** NCELLS, NPRINT, L.at(0) (using w,w0,th)
// 	** SET PBCS TO 0


void cellPacking2D::initializeHopperSP(vector<double>& radii, double w0, double w, double th, double Lmin){
	// local variables
	int ci, vi, d, nvtmp;
	double a0tmp, l0tmp, calA0tmp;
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
		// boundary information ( SET PBCS TO 1 )
		for (d=0; d<NDIM; d++){
			cell(ci).setL(d,L.at(0));
			cell(ci).setpbc(d,0);
		}

		// number of vertices ( SIGMA SETS # OF VERTS )
		nvtmp = round(2.0*radii.at(ci)*nvmin);
		if (nvtmp > nvmin)
 			cell(ci).setNV(nvtmp);
		else
			cell(ci).setNV(nvmin);

		// array information
		cell(ci).initializeVertices();
		cell(ci).initializeCell();

		// initialize cells as regular polygons
		calA0tmp = nvtmp*tan(PI/nvtmp)/PI;

		// preferred area is regular polygon with slightly smaller radius
		a0tmp = 0.5*nvtmp*pow(radii.at(ci),2.0)*sin(2.0*PI/nvtmp);

		// initial length of polygon side
		l0tmp = sqrt(4.0*PI*a0tmp*calA0tmp)/nvtmp;

		// set preferred area and length 
		cell(ci).seta0(a0tmp);
		cell(ci).setl0(l0tmp);
		cell(ci).setdel(1.0);
	}

	// initialize particle positions
	cout << "		-- Ininitializing cell positions" << endl;
	for (ci=0; ci<NCELLS; ci++){
		// set min and max values of positions
		xmin = -Lmin*L.at(0);
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



// FUNCTION TO INITALIZE PARTICLE POSITIONS AS IF THEY WERE SOFT PARTICLES
// WITH DIAMETER SIGMA
// 
// THEN turn particles to deformable particles
// 
// assume following will already be initialized:
// 	** NCELLS, NPRINT, L.at(0) (using w,w0,th)
// 	** SET PBCS TO 0


void cellPacking2D::initializeHopperDP(vector<double>& radii, double w0, double w, double th, double Lmin, int NV){
	// local variables
	int ci, vi, d, nvtmp;
	double a0tmp, l0tmp, calA0tmp;
	double xpos, ypos;
	double xmin, xmax, ymin, ymax;
	double vrmin, Ltmp;

	// minimum number of vertices
	const int nvmin = 12;

	// random number generator
	srand48(10*seed);

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
		for (d=0; d<NDIM; d++)
			cell(ci).setpbc(d,0);

		// number of vertices ( SIGMA SETS # OF VERTS )
		nvtmp = round(2.0*radii.at(ci)*NV);
		if (nvtmp > nvmin)
 			cell(ci).setNV(nvtmp);
		else
			cell(ci).setNV(nvmin);

		// array information
		cell(ci).initializeVertices();
		cell(ci).initializeCell();

		// initialize cells as regular polygons
		calA0tmp = nvtmp*tan(PI/nvtmp)/PI;

		// preferred area is regular polygon with slightly smaller radius
		a0tmp = 0.5*nvtmp*pow(radii.at(ci),2.0)*sin(2.0*PI/nvtmp);

		// initial length of polygon side
		l0tmp = sqrt(4.0*PI*a0tmp*calA0tmp)/nvtmp;

		// set preferred area and length 
		cell(ci).seta0(a0tmp);
		cell(ci).setl0(l0tmp);
		cell(ci).setdel(1.0);
	}

	// initialize L based on smallest vertex radius
	vrmin = 0.5*cell(0).getl0();
	Ltmp = 0.5*(w0 - w)*tan(th) + vrmin*((1.0/cos(th)) + 1.0 - tan(th));
	L.at(0) = Ltmp;
	for (ci=0; ci<NCELLS; ci++){
		for (d=0; d<NDIM; d++)
			cell(ci).setL(d,Ltmp);
	}

	// initialize particle positions
	cout << "		-- Ininitializing cell positions" << endl;
	for (ci=0; ci<NCELLS; ci++){
		// set min and max values of positions
		xmin = -Lmin*L.at(0);
		xmax = -radii.at(ci);
		ymin = radii.at(ci);
		ymax = w0 - radii.at(ci);

		// get random location in hopper reservoir
		xpos = (xmin-xmax)*drand48() + xmax;
		ypos = (ymax-ymin)*drand48() + ymin;

		// set as initial position of com
		cell(ci).setCPos(0,xpos);
		cell(ci).setCPos(1,ypos);
	}

	// initialize phi
	cout << "		-- Ininitializing packing fraction...";
	phi = hopperPackingFraction(radii,w0,w,th);
	cout << "which is phi = " << phi << endl;

	// initial time scales
	cout << "		-- Ininitializing time scale" << endl;
	dt = 0.05;
	dt0 = dt;

	// slightly increase radii to give more space for dp
	double rscaleForDP = 1.1;
	for (ci=0; ci<NCELLS; ci++)
		radii.at(ci) *= rscaleForDP;

	// use FIRE in hopper geometry to relax overlaps
	cout << "		-- Using FIRE to relax initial overlaps..." << endl;
	fireMinimizeHopperSP(radii,w0,w,th);

	// shrink radii back down for SP runs
	for (ci=0; ci<NCELLS; ci++)
		radii.at(ci) /= rscaleForDP;

	// update vertex positions based on cell positions
	for (ci=0; ci<NCELLS; ci++){		

		// initialize vertices as a regular polygon
		cell(ci).regularPolygon();

		// update real-space positions
		for (vi=0; vi<cell(ci).getNV(); vi++){
			for (d=0; d<NDIM; d++)
				cell(ci).setVPos(vi,d,cell(ci).cpos(d) + cell(ci).vrel(vi,d));
		}
	}
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
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << Knew << endl;
			cout << "	* Pvirial 	= " << Pvirial << endl;
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
		converged = (converged || (abs(Pvirial) > Ptol && Knew < Ktol));

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

// calculate all forces on DP model particles
// 	** if closed = 1, orifice is closed off by wall
// 	** g is force strength (i.e. generalized gravity)
void cellPacking2D::hopperForcesDP(double w0, double w, double th, double g, int closed){
	// local variables
	int ci,cj,vi,d,dd,inContact;

	// reset virial stresses to 0
	sigmaXX = 0.0;
	sigmaXY = 0.0;
	sigmaYX = 0.0;
	sigmaYY = 0.0;

	// reset contacts before force calculation
	resetContacts();
	Ncc = 0;
	Nvv = 0;

	// reset forces
	for (ci=0; ci<NCELLS; ci++){
		// reset center of mass forces
		for (d=0; d<NDIM; d++)
			cell(ci).setCForce(d,0.0);

		// reset vertex forces and interaction energy
		for (vi=0; vi<cell(ci).getNV(); vi++){
			// forces
			for (d=0; d<NDIM; d++)
				cell(ci).setVForce(vi,d,0.0);

			// energies
			cell(ci).setUInt(vi,0.0);
		}
	}


	// loop over cells and cell pairs, calculate shape and interaction forces
	for (ci=0; ci<NCELLS; ci++){
		// loop over pairs, add info to contact matrix
		for (cj=ci+1; cj<NCELLS; cj++){
			// calculate forces, add to number of vertex-vertex contacts
			inContact = cell(ci).vertexForce(cell(cj),sigmaXX,sigmaXY,sigmaYX,sigmaYY);
			if (inContact > 0){
				// add to cell-cell contacts
				addContact(ci,cj);
				Ncc++;

				// increment vertex-vertex contacts
				Nvv += inContact;
			}
		}

		// forces on vertices due to shape
		cell(ci).shapeForces();

		// gravity force (in +x direction)
		cell(ci).gravityForces(g,0);
	}

	// reset vstress to 0, for hopper sims used to compute net force on top (x) and bottom (y) walls
	sigmaXX = 0.0;
	sigmaXY = 0.0;
	sigmaYX = 0.0;
	sigmaYY = 0.0;

	// wall forces
	hopperWallForcesDP(w0,w,th,closed);
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
	double yline;							// line separating edge force from wall force

	// preliminary calculations
	t = tan(th);
	c = cos(th);
	s = sin(th);

	// reset wall forces 
	// 	** forces on top wall in x/y direction stored in sigmaXX/sigmaXY
	// 	** forces on bottom wall in x/y direction stored in sigmaYX/sigmaYY
	sigmaXX = 0.0;
	sigmaXY = 0.0;
	sigmaYX = 0.0;
	sigmaYY = 0.0;

	// loop over cells
	for (ci=0; ci<NCELLS; ci++){
		// get sigma (2*radius)
		sigma = 2*radii.at(ci);

		// get particle positions
		x = cell(ci).cpos(0);
		y = cell(ci).cpos(1);

		// check hopper walls
		if (x > -sigma*s){
			// if particle in hopper bulk
			if (x < L.at(0) - 0.5*sigma*s){
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

						// add to net force on wall
						sigmaYX -= ftmp*(lwx/lw);
						sigmaYY -= ftmp*(lwy/lw);
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

						// add to net force on wall
						sigmaXX -= ftmp*(lwx/lw);
						sigmaXY -= ftmp*(lwy/lw);
					}
				}
			}

			// if particle is near edge; either do wall force or force due to edge
			else if (x > L.at(0) - 0.5*sigma*s && x < L.at(0)){
				// check on top wall
				if (y > 0.5*w0){
					// define line separating wall force and edge force regime
					yline = 0.5*(w0 + w) - ((L.at(0) - x)/t);

					// if above yline, use wall force
					if (y > yline){
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

							// add to net force on wall
							sigmaXX -= ftmp*(lwx/lw)*lwx;
							sigmaXY -= ftmp*(lwy/lw)*lwy;
						}
					}
					// else if below yline but above radius cutoff, use edge force
					else if (y < yline && y > 0.5*(w0 + w - sigma)){
						// vector to edge point on top lip
						lwx = x - L.at(0);
						lwy = y - 0.5*(w0 + w);

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

							// add to net force on wall
							sigmaXX -= ftmp*(lwx/lw)*lwx;
							sigmaXY -= ftmp*(lwy/lw)*lwy;
						}
					}
				}
				// check on bottom wall
				else{
					// define line separating wall force and edge force regime
					yline = 0.5*(w0 - w) + ((L.at(0) - x)/t);

					// if below yline, use wall force
					if (y < yline){
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

							// add to net force on wall
							sigmaYX -= ftmp*(lwx/lw);
							sigmaYY -= ftmp*(lwy/lw);
						}
					}
					// else if above yline but below radius cutoff, use edge force
					else if (y > yline && y < 0.5*(w0 - w + sigma)){
						// vector to edge point on bottom lip
						lwx = x - L.at(0);
						lwy = y - 0.5*(w0 - w);

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

							// add to net force on wall
							sigmaYX -= ftmp*(lwx/lw);
							sigmaYY -= ftmp*(lwy/lw);
						}
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

				// add to net force on wall in y direction
				sigmaYY -= ftmp;
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

				// add to net force on wall in y direction
				sigmaXY += ftmp;
			}
		}

		// if orifice closed, check orifice wall
		if (closed == 1){
			if (x > L.at(0) - radii.at(ci)){

				// vector from particle to wall
				lwx = L.at(0) - x;

				// overlap with wall
				overlap = 2.0*lwx/sigma;

				// add to x force ONLY (points in negative x direction)
				ftmp = 1 - overlap;
				cell(ci).setCForce(0,cell(ci).cforce(0) - ftmp);

				// add to energies
				utmp = 0.25*sigma*pow(1 - overlap,2);
				for (vi=0; vi<cell(ci).getNV(); vi++)
					cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());
			}
		}
		
	}
}


// wall forces between cells as droplets (DP model)
// 	** if closed = 1, orifice is closed off by wall
void cellPacking2D::hopperWallForcesDP(double w0, double w, double th, int closed){
	// local variables
	int ci, vi; 							// indices
	double x, y;							// vertex positions
	double sigma; 							// vertex diameter
	double sb, sib, xtb, ytb, xbb, ybb;		// bead information
	double Lx, xedge;						// hopper nozzle length variables
	double t, c, s;							// tangent, cosine, sine
	double hPlus, hMinus;					// height of angled wall
	double dyPlus, dyMinus; 				// distance from top/bottom wall
	double yPlusMin, yMinusMax; 			// cutoff positions for wall forces
	double lw, lwx, lwy;					// elements of vector pointing from wall to vertex
	double overlap;							// overlap of vertex with wall
	double ftmp, utmp;						// force/energy of particle overlap with walls
	double yline;							// line separating edge force from wall force

	// hopper nozzle length
	Lx = L.at(0);

	// trig factors
	t = tan(th);
	c = cos(th);
	s = sin(th);

	// edge bead information
	sb 		= cell(0).getl0()*cell(0).getdel();
	xtb 	= Lx - 0.5*sb;
	ytb 	= 0.5*(w0 + w + sb);
	xbb 	= xtb;
	ybb 	= 0.5*(w0 - w - sb);


	// loop over cells and vertices
	for (ci=0; ci<NCELLS; ci++){
		for (vi=0; vi<cell(ci).getNV(); vi++){
			// get vertex diameter
			sigma = cell(ci).getl0()*cell(ci).getdel();

			// determine sigma_ij with edge bead
			sib = 0.5*(sigma + sb);

			// x cutoff for bead interaction
			xedge = xtb - sib*c;

			// get particle positions
			x = cell(ci).vpos(vi,0);
			y = cell(ci).vpos(vi,1);

			// check hopper walls
			if (x > -sigma*s){
				// if vertex in hopper bulk
				if (x < xedge){
					// check ymin for walls
					yPlusMin 	= w0 - (x/t) - (0.5*sigma/s);
					yMinusMax 	= (x/t) + (0.5*sigma/s);

					// if true, interacting with bottom wall
					if (y < yMinusMax){
						// distance to wall
						lw = (y - x/t)*s;

						if (lw < 0.5*sigma){
							// overlap with wall
							overlap = 2.0*lw/sigma;

							// force
							ftmp = (2.0/sigma)*(1 - overlap);
							cell(ci).setVForce(vi,0,cell(ci).vforce(vi,0) - ftmp*c);
							cell(ci).setVForce(vi,1,cell(ci).vforce(vi,1) + ftmp*s);

							// add to energies
							utmp = 0.5*pow(1 - overlap,2);
							cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp);

							// add to net force on bottom wall
							sigmaYX += ftmp*c;
							sigmaYY -= ftmp*s;
						}
					}


					// if true, interacting with top wall
					if (y > yPlusMin){
						// distance to wall
						lw = (w0 - x/t - y)*s;

						if (lw < 0.5*sigma){
							// overlap with wall
							overlap = 2.0*lw/sigma;

							// force
							ftmp = (2.0/sigma)*(1 - overlap);
							cell(ci).setVForce(vi,0,cell(ci).vforce(vi,0) - ftmp*c);
							cell(ci).setVForce(vi,1,cell(ci).vforce(vi,1) - ftmp*s);

							// add to energies
							utmp = 0.5*pow(1 - overlap,2);
							cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp);

							// add to net force on top wall
							sigmaXX += ftmp*c;
							sigmaXY += ftmp*s;
						}
					}
				}

				// if particle is near edge; either do wall force or force due to edge
				else if (x > xedge && x < Lx){
					// check on top wall
					if (y > 0.5*w0){
						// define line separating wall force and edge force regime
						yline = (x - xtb)/t + ytb;

						// if above yline, use wall force from top wall
						if (y > yline){
							// distance
							lw = (w0 - (x/t) - y)*s;

							if (lw < 0.5*sigma){
								// overlap with wall
								overlap = 2.0*lw/sigma;

								// force
								ftmp = (2.0/sigma)*(1 - overlap);
								cell(ci).setVForce(vi,0,cell(ci).vforce(vi,0) - ftmp*c);
								cell(ci).setVForce(vi,1,cell(ci).vforce(vi,1) - ftmp*s);

								// add to energies
								utmp = 0.5*pow(1 - overlap,2);
								cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp);

								// add to net force on top wall
								sigmaXX += ftmp*c;
								sigmaXY += ftmp*s;
							}
						}
						// else, check overlap with bead
						else{
							// vector to edge bead
							lwx = x - xtb;
							lwy = y - ytb;

							// distance to bead center edge
							lw = sqrt(lwx*lwx + lwy*lwy) - 0.5*sb;

							if (lw < 0.5*sigma){
								// overlap with wall (use bead edge, not center)
								overlap = 2.0*lw/sigma;

								// force
								ftmp = (2.0/sigma)*(1 - overlap);
								cell(ci).setVForce(vi,0,cell(ci).vforce(vi,0) + ftmp*(lwx/lw));
								cell(ci).setVForce(vi,1,cell(ci).vforce(vi,1) + ftmp*(lwy/lw));

								// add to energies
								utmp = 0.5*pow(1 - overlap,2);
								cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp);

								// add to net force on top wall
								sigmaXX -= ftmp*(lwx/lw);
								sigmaXY -= ftmp*(lwy/lw);
							}
						}
					}
					// check on bottom wall
					else{
						// define line separating wall force and edge force regime
						yline = (xbb - x)/t + ybb;

						// if below yline, use wall force
						if (y < yline){
							// distance to wall
							lw = (y - (x/t))*s;

							if (lw < 0.5*sigma){
								// overlap with wall
								overlap = 2.0*lw/sigma;

								// force
								ftmp = (2.0/sigma)*(1 - overlap);
								cell(ci).setVForce(vi,0,cell(ci).vforce(vi,0) - ftmp*c);
								cell(ci).setVForce(vi,1,cell(ci).vforce(vi,1) + ftmp*s);

								// add to energies
								utmp = 0.5*pow(1 - overlap,2);
								cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp);

								// add to net force on bottom wall
								sigmaYX += ftmp*c;
								sigmaYY -= ftmp*s;
							}
						}
						// else, check overlap with bottom bead						
						else{
							// vector to bottom bead
							lwx = x - xbb;
							lwy = y - ybb;

							// distance
							lw = sqrt(lwx*lwx + lwy*lwy) - 0.5*sb;

							if (lw < 0.5*sigma){
								// overlap with wall (use bead edge, not center)
								overlap = 2.0*lw/sigma;

								// force
								ftmp = (2.0/sigma)*(1 - overlap);
								cell(ci).setVForce(vi,0,cell(ci).vforce(vi,0) + ftmp*(lwx/lw));
								cell(ci).setVForce(vi,1,cell(ci).vforce(vi,1) + ftmp*(lwy/lw));

								// add to energies
								utmp = 0.5*pow(1 - overlap,2);
								cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp);

								// add to net force on bottom wall
								sigmaYX -= ftmp*(lwx/lw);
								sigmaYY -= ftmp*(lwy/lw);
							}
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
					ftmp = (2.0/sigma)*(1 - overlap);
					cell(ci).setVForce(vi,1,cell(ci).vforce(vi,1) + ftmp);

					// add to energies
					utmp = 0.5*pow(1 - overlap,2);
					cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp);	

					// add to net force on bottom wall
					sigmaYY -= ftmp;
				}

				// if true, interacting with top wall
				if (y > yPlusMin){
					// vector from particle to wall
					lwy = w0 - y;

					// overlap with wall
					overlap = 2.0*lwy/sigma;

					// add to y force ONLY (points in negative y direction)
					ftmp = (2.0/sigma)*(1 - overlap);
					cell(ci).setVForce(vi,1,cell(ci).vforce(vi,1) - ftmp);

					// add to energies
					utmp = 0.5*pow(1 - overlap,2);
					cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp);	

					// add to net force on top wall
					sigmaXY += ftmp;
				}
			}

			// if orifice closed, check orifice wall
			if (closed == 1){
				if (x > Lx - 0.5*sigma){

					// vector from particle to wall
					lwx = Lx - x;

					// overlap with wall
					overlap = 2.0*lwx/sigma;

					// add to x force ONLY (points in negative x direction)
					ftmp = (2.0/sigma)*(1 - overlap);
					cell(ci).setVForce(vi,0,cell(ci).vforce(vi,0) - ftmp);

					// add to energies
					utmp = 0.5*pow(1 - overlap,2);
					cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp);	
				}
			}

			// add vertical edge wall
			if (x > Lx && x < Lx + 0.5*sigma){
				// check if overlap with edge bead or wall
				if (y > 0.5*w0){
					// check top bead or vert wall
					if (y < ytb){
						// interacting with top edge bead but outside of hopper

						// vector to edge bead
						lwx = x - xtb;
						lwy = y - ytb;

						// distance
						lw = sqrt(lwx*lwx + lwy*lwy) - 0.5*sb;

						if (lw < 0.5*sigma){
							// overlap with wall
							overlap = 2.0*lw/sigma;

							// force
							ftmp = (2.0/sigma)*(1 - overlap);
							cell(ci).setVForce(vi,0,cell(ci).vforce(vi,0) + ftmp*(lwx/lw));
							cell(ci).setVForce(vi,1,cell(ci).vforce(vi,1) + ftmp*(lwy/lw));

							// add to energies
							utmp = 0.5*pow(1 - overlap,2);
							cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp);

							// add to net force on top wall
							sigmaXX -= ftmp*(lwx/lw);
							sigmaXY -= ftmp*(lwy/lw);
						}
					}
					else{
						// interacting with vertical wall, force only in +x direction

						// vector from wall to particle
						lwx = x - Lx;

						// overlap with wall
						overlap = 2.0*lwx/sigma;

						// add to y force ONLY (points in positive y direction)
						ftmp = (2.0/sigma)*(1 - overlap);
						cell(ci).setVForce(vi,0,cell(ci).vforce(vi,0) + ftmp);

						// add to energies
						utmp = 0.5*pow(1 - overlap,2);
						cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp);	

						// add to net force on top wall
						sigmaXX -= ftmp;
					}
				}
				else if (y < 0.5*w0){
					// check if bottom bead or vert wall
					if (y > ybb){
						// vector to bottom bead
						lwx = x - xbb;
						lwy = y - ybb;

						// distance
						lw = sqrt(lwx*lwx + lwy*lwy) - 0.5*sb;

						if (lw < 0.5*sigma){
							// overlap with wall
							overlap = 2.0*lw/sigma;

							// force
							ftmp = (2.0/sigma)*(1 - overlap);
							cell(ci).setVForce(vi,0,cell(ci).vforce(vi,0) + ftmp*(lwx/lw));
							cell(ci).setVForce(vi,1,cell(ci).vforce(vi,1) + ftmp*(lwy/lw));

							// add to energies
							utmp = 0.5*pow(1 - overlap,2);
							cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp);

							// add to net force on bottom wall
							sigmaYX -= ftmp*(lwx/lw);
							sigmaYY -= ftmp*(lwy/lw);
						}
					}
					else{
						// interacting with vertical wall, force only in +x direction

						// vector from wall to particle
						lwx = x - Lx;

						// overlap with wall
						overlap = 2.0*lwx/sigma;

						// add to y force ONLY (points in positive y direction)
						ftmp = (2.0/sigma)*(1 - overlap);
						cell(ci).setVForce(vi,0,cell(ci).vforce(vi,0) + ftmp);

						// add to energies
						utmp = 0.5*pow(1 - overlap,2);
						cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp);	

						// add to net force on bottom wall
						sigmaYX -= ftmp;
					}
				}
			}


		}
	}
}







// function to run NVE dynamics on DP particles using velocity-verlet to check energy conservation
void cellPacking2D::hopperDPNVE(double w0, double w, double th, double g, double T0){
	// local variables
	int closed = 0;
	int t, ci, vi, vip1, nvtmp, d;
	double Pvirial, K;
	double aH, aR, aT;
	double cxtmp, xi, xip1, yi, yip1, utmp;

	// reservoir area
	aR = w0*L.at(0);

	// hopper area
	aH = w*w0 + pow(w0-w,2)/(4.0*tan(th));

	// total area
	aT = aR + aH;

	// check that NT has been set 
	if (NT <= 0){
		cout << "	** ERROR: in hopper DP NVE, sim length NT = " << NT << ", which is <= 0. ending." << endl;
		exit(1);
	}

	// initialize velocities
	initializeVelocities(T0);

	// loop over time, run NVE dynamics
	for (t=0; t<NT; t++){
		// update virial pressure
		Pvirial = 0.5*(sigmaXX + sigmaYY)/(NCELLS*aT);

		// update kinetic energy based on com velocity
		K = totalKineticEnergy();

		// output some information to console
		if (t % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 		Deformable Particle NVE in hopper geometry, t = " << t << endl << endl;
			cout << "===================================================" << endl;

			// add gravitiational potential energy to uint
			for (ci=0; ci<NCELLS; ci++){
				// get com position
				cxtmp = 0.0;
				nvtmp = cell(ci).getNV();
				for (vi=0; vi<nvtmp; vi++){
					// next index
					vip1 = (vi+1) % nvtmp;

					// get vertex coordinates
					xi = cell(ci).vpos(vi,0);
					xip1 = cell(ci).vpos(vip1,0);

					yi = cell(ci).vpos(vi,1);
					yip1 = cell(ci).vpos(vip1,1);

					// add to com
					cxtmp += (1.0/6.0)*((xi + xip1)*(xi*yip1 - xip1*yi));
				}

				// add to potential energy
				utmp = -g*cxtmp;
				for (vi=0; vi<nvtmp; vi++)
					cell(ci).setUInt(vi,cell(ci).uInt(vi) + (utmp/nvtmp));
			}

			// print if object has been opened already
			if (packingPrintObject.is_open()){
				cout << "	* Printing DP center positions to file" << endl;
				printHopperDP(w0, w, th);
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing DP energy to file" << endl;
				printSystemEnergy();
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* U 		= " << totalPotentialEnergy() << endl;
			cout << "	* K 		= " << totalKineticEnergy() << endl;
			cout << "	* E 		= " << totalPotentialEnergy() + totalKineticEnergy() << endl;
			cout << "	* virial P 	= " << Pvirial << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << endl << endl;
		}

		// VV update in FIRE 2.0: position update
		for (ci=0; ci<NCELLS; ci++){
			cell(ci).verletPositionUpdate(dt);
			cell(ci).updateCPos();
		}

		// calculate forces
		hopperForcesDP(w0, w, th, g, closed);

		// VV update in FIRE 2.0: Velocity update 2
		for (ci=0; ci<NCELLS; ci++)
			cell(ci).verletVelocityUpdate(dt);
	}
}


// function to run NVE dynamics on SP particles using velocity-verlet to check energy conservation
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
	while (xmax < L.at(0) - 2.5*radii.at(0) && itr < itrMax){
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
				printHopperSP(radii,w0,w,th,g);
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing SP energy to file" << endl;
				printSystemEnergy(itr,Pvirial,K);
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << K << endl;
			cout << "	* Ftop  	= " << sqrt(sigmaXX*sigmaXX + sigmaXY*sigmaXY) << endl;
			cout << "	* Fbottom  	= " << sqrt(sigmaYX*sigmaYX + sigmaYY*sigmaYY) << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << "	* g 		= " << g << endl;
			cout << "	* xmax/L 	= " << xmax/L.at(0) << endl;
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
		xmin = 2*L.at(0);
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

				if (d == 0 && postmp > L.at(0) + 2.0*radii.at(ci)){
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
void cellPacking2D::flowHopperDP(double w0, double w, double th, double g, double b){
	// local variables
	int closed = 0;
	int t, ci, cj, vi, vip1, nvtmp, d;
	double Pvirial, K;
	double aH, aR, aT;
	double cxtmp, cytmp, xi, xip1, yi, yip1, utmp;

	// replacement check
	double minx, dx, dy, wmax, wmin;

	// reservoir area
	aR = w0*L.at(0);

	// hopper area
	aH = w*w0 + pow(w0-w,2)/(4.0*tan(th));

	// total area
	aT = aR + aH;

	// check that NT has been set 
	if (NT <= 0){
		cout << "	** ERROR: in hopper DP NVE, sim length NT = " << NT << ", which is <= 0. ending." << endl;
		exit(1);
	}

	// loop over time, run NVE dynamics
	for (t=0; t<NT; t++){
		// update virial pressure
		Pvirial = 0.5*(sigmaXX + sigmaYY)/(NCELLS*aT);

		// update kinetic energy based on com velocity
		K = totalKineticEnergy();

		// output some information to console
		if (t % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 		Deformable Particle flow in hopper geometry, t = " << t << endl << endl;
			cout << "===================================================" << endl;

			// add gravitiational potential energy to uint
			for (ci=0; ci<NCELLS; ci++){
				// get com position
				cxtmp = 0.0;
				nvtmp = cell(ci).getNV();
				for (vi=0; vi<nvtmp; vi++){
					// next index
					vip1 = (vi+1) % nvtmp;

					// get vertex coordinates
					xi = cell(ci).vpos(vi,0);
					xip1 = cell(ci).vpos(vip1,0);

					yi = cell(ci).vpos(vi,1);
					yip1 = cell(ci).vpos(vip1,1);

					// add to com
					cxtmp += (1.0/6.0)*((xi + xip1)*(xi*yip1 - xip1*yi));
				}

				// add to potential energy
				utmp = -g*cxtmp;
				for (vi=0; vi<nvtmp; vi++)
					cell(ci).setUInt(vi,cell(ci).uInt(vi) + (utmp/nvtmp));
			}

			// print if object has been opened already
			if (packingPrintObject.is_open()){
				cout << "	* Printing DP center positions to file" << endl;
				printHopperDP(w0, w, th);
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing DP energy to file" << endl;
				printSystemEnergy();
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* U 		= " << totalPotentialEnergy() << endl;
			cout << "	* K 		= " << totalKineticEnergy() << endl;
			cout << "	* E 		= " << totalPotentialEnergy() + totalKineticEnergy() << endl;
			cout << "	* virial P 	= " << Pvirial << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << endl << endl;
		}

		// VV position update
		for (ci=0; ci<NCELLS; ci++){
			cell(ci).verletPositionUpdate(dt);
			cell(ci).updateCPos();

			// replace in back of hopper
			if (cell(ci).cpos(0) > L.at(0) + 2.0*sqrt(cell(ci).geta0())){
				// find minimum xcoordinate
				minx = L.at(0);
				for (cj=0; cj<NCELLS; cj++){
					if (cell(cj).cpos(0) < minx)
						minx = cell(cj).cpos(0);
				}

				// place particle behind minimum
				cxtmp = minx - 1.5*sqrt(cell(ci).geta0());

				// give random y (if cx tmp is in nozzle, use ybounds from nozzle heights)
				if (cxtmp < 0)
					cytmp = (w0 - 2.0*sqrt(cell(ci).geta0()))*drand48() + sqrt(cell(ci).geta0());
				else{
					wmin = cxtmp*tan(th) + sqrt(cell(ci).geta0());
					wmax = w0 - wmin - sqrt(cell(ci).geta0());
					cytmp = (wmax - wmin)*drand48() + wmin;
				}


				// compute change in position
				dx = cxtmp - cell(ci).cpos(0);
				dy = cytmp - cell(ci).cpos(1);

				// update real-space positions
				for (vi=0; vi<cell(ci).getNV(); vi++){
					cell(ci).setVPos(vi,0,cell(ci).vpos(vi,0) + dx);
					cell(ci).setVPos(vi,1,cell(ci).vpos(vi,1) + dy);
				}
				cell(ci).updateCPos();
			}
		}
		

		// calculate forces
		hopperForcesDP(w0, w, th, g, closed);

		// VV update 2 with damping
		for (ci=0; ci<NCELLS; ci++)
			cell(ci).verletVelocityUpdate(dt,b);
	}
}


// calculate packing fraction in hopper geometry (GIVEN RADII)
double cellPacking2D::hopperPackingFraction(vector<double>& radii, double w0, double w, double th){
	// local variables
	int ci;
	double aH, aR, aT;

	// reservoir area
	aR = w0*L.at(0);

	// hopper area
	aH = w*w0 + pow(w0-w,2)/(4.0*tan(th));

	// total area
	aT = aR + aH;

	// initialize phi to 0
	phi = 0.0;

	// loop over radii of cells in flow part of hopper
	for (ci=0; ci<NCELLS; ci++){
		if (cell(ci).cpos(0) > -L.at(0))
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
void cellPacking2D::printHopperDP(double w0, double w, double th){
	// local variables
	int ci,d;
	int w1 = 12;
	int w2 = 6;
	int w3 = 30;

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
	packingPrintObject << setw(w1) << right << w0;
	packingPrintObject << setw(w1) << right << w;
	packingPrintObject << setw(w1) << setprecision(6) << right << th;
	packingPrintObject << endl;

	// print stress information
	packingPrintObject << setw(w1) << left << "WLFRC";
	packingPrintObject << setw(w3) << right << sigmaXX;
	packingPrintObject << setw(w3) << right << sigmaXY;
	packingPrintObject << setw(w3) << right << sigmaYX;
	packingPrintObject << setw(w3) << right << sigmaYY;
	packingPrintObject << endl;

	// print contact information
	packingPrintObject << setw(w1) << left << "NCNTS";
	packingPrintObject << setw(w1) << right << Ncc;
	packingPrintObject << setw(w1) << right << Nvv;
	packingPrintObject << endl;

	// print contact matrix
	packingPrintObject << setw(w1) << left << "CMMAT";

	// print contact matrix
	for (int ci=0; ci<NCELLS; ci++){
		for (int cj=ci+1; cj<NCELLS; cj++)
			packingPrintObject << setw(w2) << contacts(ci,cj);
	}
	packingPrintObject << endl;

	// print info for rest of the cells
	for (int ci=0; ci<NCELLS; ci++)
		cell(ci).printVertexPositions(packingPrintObject,ci);

	// print end frame
	packingPrintObject << setw(w1) << left << "ENDFR" << " " << endl;
}










