/*

	Methods file for 2D sticky particles

*/


// include file
#include "deformableParticles2D.h"
#include "cellPacking2D.h"

// namespace
using namespace std;

// constants
const double PI = 4*atan(1);

// 2D STICKY PARTICLE FUNCTIONS

// FUNCTION TO INITALIZE PARTICLE POSITIONS AS IF THEY WERE SOFT PARTICLES
// WITH DIAMETER SIGMA
// 
// 	NOTE: radii are input in units of sigma, the mean particle diameter
void cellPacking2D::initializeStickySP(vector<double>& radii, double phiDisk, double sizeDispersion){
	// local variables
	int ci, vi, d, nvtmp;
	double r1, r2, g1, radsum;
	double xpos, ypos;
	double xmin, xmax, ymin, ymax;
	double calA;

	// minimum number of vertices
	const int nvmin = 12;
	const int NV = 24;

	// check inputs to stick SP initialization
	if (radii.size() < NCELLS){
		cout << "	** ERROR: in initializing sticky SP, input radii vector size = " << radii.size() << ", which is != NCELLS (= " << NCELLS << "). ending." << endl;
		exit(1);
	}

	// output to console
	cout << "		-- In stickySP initialization, initializing SP particles" << endl;

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
	for (d=0; d<NDIM; d++)
		L.at(d) = sqrt(PI*radsum/phiDisk);

	// set phi to input
	phi = phiDisk;

	// reseed rng
	srand48(56835698*seed);

	// initialize cell information
	cout << "		-- Ininitializing cell objects" << endl;
	for (ci=0; ci<NCELLS; ci++){
		// boundary information ( SET PBCS TO 1 )
		for (d=0; d<NDIM; d++){
			cell(ci).setL(d,L.at(d));
			cell(ci).setpbc(d,1);
		}

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
		xmin = radii.at(ci);
		xmax = L.at(0) - radii.at(ci);
		ymin = radii.at(ci);
		ymax = L.at(1) - radii.at(ci);

		// get random location in hopper RESERVOIR
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
	dt = 0.1*sqrt(0.25*PI);
	dt0 = dt;

	// use FIRE in hopper geometry to relax overlaps
	cout << "		-- Using FIRE to relax overlaps..." << endl;
	fireMinimizeSP(radii,0.0);
}

void cellPacking2D::stickySPTriangularLattice(vector<double>& radii, double phiDisk){
	// local variables
	int ci, d, nvtmp;
	double xmin, xmax, ymin, ymax;
	double xpos, ypos, calA, rad, radsum;

	// output to console
	cout << "		-- In stickySP initialization, initializing on triangular lattice" << endl;

	// minimum number of vertices
	const int nvmin = 12;
	const int NV = 24;

	// ratio of box length in y direction
	const double ly = 0.5*sqrt(3);

	// loop over radii (all have unit diameters)
	rad = 0.5;
	for (ci=0; ci<NCELLS; ci++)
		radii.at(ci) = 0.5;

	// get radii^2 sum
	radsum = NCELLS*pow(rad,2.0);

	// determine box length from particle sizes and input packing fraction
	L.at(0) = sqrt(PI*radsum/(ly*phiDisk));
	L.at(1) = sqrt(ly*PI*radsum/phiDisk);

	// initialize cell information
	cout << "		-- Ininitializing cell objects" << endl;
	for (ci=0; ci<NCELLS; ci++){
		// boundary information ( SET PBCS TO 1 )
		for (d=0; d<NDIM; d++){
			cell(ci).setL(d,L.at(d));
			cell(ci).setpbc(d,1);
		}

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
		xmin = radii.at(ci);
		xmax = L.at(0) - radii.at(ci);
		ymin = radii.at(ci);
		ymax = L.at(1) - radii.at(ci);

		// place on triangular lattice
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
	dt = 0.1*sqrt(0.25*PI);
	dt0 = dt;

	// use FIRE in hopper geometry to relax overlaps
	cout << "		-- Using FIRE to relax overlaps..." << endl;
	fireMinimizeSP(radii,0.0);

	// Use NVE to jumble positions
	initializeVelocities(10.0);
	spNVE(radii,1000);

	// use FIRE in hopper geometry to relax overlaps
	cout << "		-- Using FIRE to relax to new state..." << endl;
	fireMinimizeSP(radii,0.0);
}

void cellPacking2D::fireMinimizeSP(vector<double>& lenscales, double attractiveParam){
	// HARD CODE IN FIRE PARAMETERS
	const double alpha0 	= 0.25;
	const double finc 		= 1.01;
	const double fdec 		= 0.5;
	const double falpha 	= 0.99;
	const double dtmax 		= 10*dt0;
	const double dtmin 		= 0.02*dt0;
	const int NMIN 			= 2000;
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

	// initialize forces (neglect damping forces, only interactions)
	resetContacts();
	spForces(lenscales);

	// initialize virial pressure from pressure from last time
	Pvirial = 0.5*(sigmaXX + sigmaYY);

	// update kinetic energy based on com velocity
	Knew = 0.0;
	for (ci=0; ci<NCELLS; ci++)
		Knew += 0.5*(PI*pow(lenscales.at(ci),2))*sqrt(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2));

	// calc check variables
	Kcheck = Knew/(NDIM*NCELLS);
	Pcheck = Pvirial/(NDIM*NCELLS);

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
				printPositionsStickySP(lenscales);
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
			cout << "	* K 		= " << Kcheck << endl;
			cout << "	* Pvirial 	= " << Pcheck << endl;
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
		spPosVerlet();

		// reset contacts before force calculation
		resetContacts();

		// calculate forces between disks (with door closed)
		spAttractiveForces(lenscales,attractiveParam);

		// verlet velocity update
		spVelVerlet(lenscales);

		// update t
		t += dt;

		// update virial pressure
		Pvirial = 0.5*(sigmaXX + sigmaYY);

		// update kinetic energy based on com velocity
		Knew = 0.0;
		for (ci=0; ci<NCELLS; ci++)
			Knew += 0.5*(PI*pow(lenscales.at(ci),2))*sqrt(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2));

		// update if Pvirial under tol
		if (abs(Pvirial) < Ptol)
			npPMIN++;
		else
			npPMIN = 0;

		// calc check variables
		Kcheck = Knew/(NDIM*NCELLS);
		Pcheck = Pvirial/(NDIM*NCELLS);

		// check for convergence
		converged = (abs(Pcheck) < Pvirial && npPMIN > NMIN);
		converged = (converged || (abs(Pcheck) > Ptol && Kcheck < Ktol));

		if (converged){
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

void cellPacking2D::spAttractiveForces(vector<double>& lenscales, double attractiveParam){
	// local variables
	int ci, cj, vi, d;
	double contactDistance = 0.0; 
	double centerDistance = 0.0;
	double energyScale;
	double overlap = 0.0;
	double uv = 0.0;
	double ftmp, utmp;
	double p1 = 1.0 + attractiveParam;
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
			contactDistance = lenscales.at(ci) + lenscales.at(cj);

			// center-to-center distance
			centerDistance = 0.0;
			for (d=0; d<NDIM; d++){
				// vectorial quantity
				distanceVec.at(d) = cell(ci).cellDistance(cell(cj),d);

				// add to distance
				centerDistance += pow(distanceVec.at(d),2);
			}

			centerDistance = sqrt(centerDistance);

			// check if within interaction zone
			if (centerDistance < contactDistance*p1){
				// add to contact checking
				addContact(ci,cj);

				// overlap scale
				overlap = centerDistance/contactDistance;

				// energy scale
				energyScale = contactDistance;

				// IF in zone to use repulsive force (and, if attractiveParam > 0, bottom of attractive well)
				if (centerDistance < contactDistance){

					// add to potential energy (energy should increase because particles are growing)
					utmp = 0.5 * energyScale * pow(1 - overlap,2) - (energyScale*attractiveParam*attractiveParam)/6.0;
					for (vi=0; vi<cell(ci).getNV(); vi++)
						cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());
					for (vi=0; vi<cell(cj).getNV(); vi++)
						cell(cj).setUInt(vi,cell(cj).uInt(vi) + utmp/cell(cj).getNV());

					// add to forces
					ftmp = 1.0 - overlap;
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

				// IF attractiveParam > 0, in attractive well
				else if (centerDistance >= contactDistance && centerDistance < contactDistance*p1 && attractiveParam > 0.0){

					// add to potential energy (energy should increase because particles are growing)
					utmp = (-energyScale/(6.0*attractiveParam)) * pow(overlap - 1.0,2) * (2.0*(overlap - 1.0) - 3.0*attractiveParam) - (energyScale*attractiveParam*attractiveParam)/6.0;
					for (vi=0; vi<cell(ci).getNV(); vi++)
						cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());
					for (vi=0; vi<cell(cj).getNV(); vi++)
						cell(cj).setUInt(vi,cell(cj).uInt(vi) + utmp/cell(cj).getNV());

					// add to forces
					ftmp = (1.0/attractiveParam) * (overlap - 1.0) * (overlap - p1);
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
	}
}


void cellPacking2D::stickySPGelationQS(vector<double>& radii, double phiGel, double dphiGel, double attractiveParam){
	// local variables
	int ci;
	double phitmp, phiSc, t;

	// initial relaxation
	fireMinimizeSP(radii,attractiveParam);

	// initial phi
	phitmp = 0.0;
	for (ci=0; ci<NCELLS; ci++)
		phitmp += PI*radii.at(ci)*radii.at(ci)/(L.at(0)*L.at(1));
	phi = phitmp;

	// loop time until packing fraction below threshold
	t = 0.0;
	while (phitmp > phiGel){

		// decrease phi
		phiSc = pow((phi - dphiGel)/phi,1.0/NDIM);
		for (ci=0; ci<NCELLS; ci++)
			radii.at(ci) *= phiSc;

		// updated phi
		phitmp = 0.0;
		for (ci=0; ci<NCELLS; ci++)
			phitmp += PI*radii.at(ci)*radii.at(ci)/(L.at(0)*L.at(1));
		phi = phitmp;

		// print statement
		cout << "===================================================" << endl << endl;
		cout << " 	CHANGING PACKING FRACTION FROM " << phi + dphiGel << " to " << phi << ": t = " << t << endl << endl;
		cout << "===================================================" << endl;
		
		cout << "	* Run data:" << endl;
		cout << "	* old phi 	= " << phi + dphiGel << endl;
		cout << "	* new phi 	= " << phi << endl;
		cout << "	* dt 		= " << dt << endl;
		cout << "	* nc 		= " << totalNumberOfContacts() << endl;
		cout << endl << endl;

		// relax overlaps
		fireMinimizeSP(radii,attractiveParam);

		// increment
		t += dt;
	}
}


void cellPacking2D::stickySPGelationRate(vector<double>& radii, double phiGel, double gelRate, double attractiveParam, double timeStepMag){
	// local variables
	int ci, vi, d, k;
	double phi0, phitmp, phiSc, t;
	double veltmp;

	// initial relaxation
	fireMinimizeSP(radii,attractiveParam);

	// initial phi
	phitmp = 0.0;
	for (ci=0; ci<NCELLS; ci++)
		phitmp += PI*radii.at(ci)*radii.at(ci)/(L.at(0)*L.at(1));
	phi = phitmp;
	phi0 = phi;

	// set time step
	dt = timeStepMag*sqrt(0.25*PI);

	// loop time until packing fraction below threshold
	t = 0.0;
	k = 0;
	while (phitmp > phiGel){
		// decrease phi according to area expansion
		phitmp = phi0/(gelRate*t + 1.0);

		// scale lengths
		phiSc = pow(phitmp/phi,1.0/NDIM);
		for (ci=0; ci<NCELLS; ci++)
			radii.at(ci) *= phiSc;

		// update phitmp
		phitmp = 0.0;
		for (ci=0; ci<NCELLS; ci++)
			phitmp += PI*radii.at(ci)*radii.at(ci)/(L.at(0)*L.at(1));
		
		// print statement
		if (k % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 	CHANGING PACKING FRACTION FROM " << phi << " to " << phitmp << ": t = " << t << endl << endl;
			cout << "===================================================" << endl;
			cout << "	* Run data:" << endl;
			cout << "	* old phi 	= " << phi << endl;
			cout << "	* new phi 	= " << phitmp << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << "	* nc 		= " << totalNumberOfContacts() << endl;
			cout << endl << endl;

			// print config and energy
			if (packingPrintObject.is_open()){
				cout << "	* Printing sticky SP positions to file" << endl;
				printPositionsStickySP(radii);
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing sticky SP energy to file" << endl;
				printEnergyStickySP();
			}

			if (statPrintObject.is_open()){
				cout << "	* Printing sticky SP contact network to file" << endl;
				printSystemContacts();
			}
		}

		// update phi
		phi = phitmp;

		// update iterator
		k++;

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
		spAttractiveForces(radii,attractiveParam);

		// increment
		t += dt;
	}
}



/*


	PRINT FUNCTIONS


*/

// function to print sticky SP position information to file
void cellPacking2D::printPositionsStickySP(vector<double>& radii){
	// local variables
	int ci,d;
	int w1 = 25;
	int w2 = 15;
	double cpostmp;

	// check to see if file is open
	if (!packingPrintObject.is_open()) {
		cout << "	ERROR: packingPrintObject is not open in printPositionsStickySP(), ending." << endl;
		exit(1);
	}

	// print information starting information
	packingPrintObject << setw(w1) << left << "NEWFR" << " " << endl;
	packingPrintObject << setw(w1) << left << "NUMCL" << setw(w2) << right << NCELLS << endl;

	// print hopper information
	packingPrintObject << setw(w1) << left << "BOXSZ";
	packingPrintObject << setw(w2) << right << L.at(0);
	packingPrintObject << setw(w2) << right << L.at(1);
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
		for (d=0; d<NDIM; d++){
			cpostmp = cell(ci).cpos(d);
			if (cpostmp > L.at(d)){
				while(cpostmp > L.at(d))
					cpostmp -= L.at(d);
			}
			else if (cpostmp < 0){
				while(cpostmp < 0)
					cpostmp += L.at(d);
			}

			packingPrintObject << setw(w2) << right << cpostmp;
		}

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

// function to print sticky SP position information to file
void cellPacking2D::printEnergyStickySP(){
	// check to see if file is open
	if (!energyPrintObject.is_open()) {
		cout << "	ERROR: energyPrintObject is not open in printSystemEnergy(), ending." << endl;
		exit(1);
	}

	// loop over particles, print cell energy
	energyPrintObject << setw(30) << setprecision(16) << right << interactionPotentialEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << totalKineticEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaXX;
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaXY;
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaYX;
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaYY;
	energyPrintObject << setw(30) << setprecision(16) << right << packingFraction();
	energyPrintObject << endl;
}







