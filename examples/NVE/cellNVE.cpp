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
const int NT 					= 2e4;			// number of time steps
const int NPRINT				= 100;			// number of steps between prints

// simulation constants
const double timeStepMag 		= 0.01;			// time step in MD units (zeta * lenscale / forcescale)

// disk constants
const double phiDisk	 		= 0.3;			// initial packing fraction of disks (sets boundary)

// force parameters
const double kl 				= 1.0;				// perimeter force constant
const double ka 				= 1.0;				// area force constant
const double gam 				= 0.0;				// surface tension force constant
const double kb 				= 1.0;				// bending energy constant
const double kint 				= 0.1;				// interaction energy constant
const double del 				= 1.0;				// width of vertices in units of l0, vertex sep on regular polygon
const double a 					= 0.0;				// attraction parameter

// deformability
const double calA0 				= 1.01;				// ratio of preferred perimeter^2 to preferred area

// main function
int main()
{
	// local variables

	// output files
	string posFile = "pos.test";
	string enFile = "en.test";

	// system details
	int NCELLS 		= 4;
	int NV			= 24;
	int seed 		= 1;
	double Ltmp 	= 1.0;
	double T0 		= 1e-3;

	double sizeRatio = 1.4;
	double sizeFraction = 0.5;

	// instantiate object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed); 	

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initializing gel at phiDisk = " << phiDisk << " using SP model" << endl;
	packingObject.initializeBidisperse(NV, phiDisk, sizeRatio, sizeFraction, del);

	// set deformability, force values
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





