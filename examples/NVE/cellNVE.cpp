// EXAMPLE File to check DPM energy conservation in NVE

// include files
#include "cellPacking2D.h"
#include "deformableParticles2D.h"
#include <sstream>

// use std name space
using namespace std;

// define PI
const double PI = 4.0*atan(1);

// length paramaters
const int NT 					= 1e5;				// number of time steps
const int NPRINT				= 100;				// number of steps between prints

// simulation constants
const double timeStepMag 		= 0.01;				// time step in MD units (zeta * lenscale / forcescale)

// force parameters
const double kl 				= 1.0;				// perimeter force constant
const double ka 				= 1.0;				// area force constant
const double kb 				= 0.01;				// bending energy constant
const double kint 				= 0.1;				// interaction energy constant
const double gam 				= 0.0;				// surface tension force constant
const double del 				= 1.0;				// vertex diameter in units of l0
const double a 					= 0.0;				// inter-vertex attraction parameter

// main function
int main()
{
	// local variables

	// output files
	string posFile = "pos.test";
	string enFile = "en.test";

	// initialize system variables
	int NCELLS, NV, seed;
	double T0, phiDisk, sizeRatio, sizeFraction, calA0;

	// set system variables
	NCELLS 			= 32;		// number of particles
	NV 				= 24;		// number of vertices on smaller particle
	seed 			= 1;		// initial seed for random number generator

	T0 				= 1e-4;		// initial temperature
	phiDisk 		= 0.7;		// initial packing fraction
	sizeRatio 		= 1.4;		// ratio of large radii to small radii
	sizeFraction	= 0.5;		// fraction of small particles
	calA0 			= 1.1;		// preferred shape parameter

	// instantiate object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,1.0,seed); 	

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initializing bidisperse system at phiDisk = " << phiDisk << " using SP model" << endl;
	packingObject.initializeBidisperse(NV, phiDisk, sizeRatio, sizeFraction, del);

	// set preferred shape parameter, force values
	packingObject.forceVals(calA0,ka,kl,gam,kb,kint,del,a);

	// update time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// open position output file
	packingObject.openPackingObject(posFile);
	packingObject.openEnergyObject(enFile);

	// run test NVE
	packingObject.initializeVelocities(T0);
	packingObject.cellNVE();

	// end
	return 0;
}





