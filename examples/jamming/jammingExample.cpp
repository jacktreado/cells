/*

	Example .cpp file to compress NCELLS to jamming
	and to compute shear modulus

*/

// include files
#include "cellPacking2D.h"
#include "deformableParticles2D.h"
#include <sstream>

// use std name space
using namespace std;

// define PI
const double PI = 4.0*atan(1);

// simulation constants
const int NT 					= 1e7; 			// number of time steps
const int NPRINT 				= 5e2;			// number of time steps between prints
const double timeStepMag 		= 0.03;			// time step in MD unit
const double phiDisk 			= 0.5;			// initial phi of SP disks
const double deltaPhi0 			= 2e-3;			// initial delta phi
const double sizeRatio 			= 1.4;			// ratio between small and large particles
const double sizeFraction		= 0.5;			// fraction of small particles

// force parameters
const double ka 			= 1.0;			// area force constant (should be = 1)
const double gam 			= 0.0;			// surface tension force constant
const double kint 			= 1.0;			// interaction energy constant
const double a 				= 0.0;			// attraction parameter 
const double del 			= 1.0;			// radius of vertices in units of l0

// tolerances
const double Ftol 			= 1e-12;		// force tolerance (for FIRE min)
const double Ptol 			= 1e-8;			// pressure tolerance

// main function
int main()
{
	// local variables

	// output files
	string posFile = "pos_sample.test";
	string enFile = "en_sample.test";
	string jamFile = "jam_sample.test";
	string vdosFile = "vdos_sample.test";

	// system details
	int NCELLS 		= 8;
	int NV			= 16;
	int seed 		= 1;
	double Ltmp 	= 1.0;

	// mechanical parameters
	double kl = 1.0;
	double kb = 0.0;
	double calA0 = 1.04;

	// instantiate main packing object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,1.0,seed);

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initializing gel at phiDisk = " << phiDisk << " using SP model" << endl;
	packingObject.initializeBidisperse(NV, phiDisk, sizeRatio, sizeFraction, del);

	// set deformability, force values
	packingObject.forceVals(calA0,ka,kl,gam,kb,kint,del,a);

	// update time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// open output files
	packingObject.openPackingObject(posFile);
	packingObject.openEnergyObject(enFile);
	packingObject.openJamObject(jamFile);
	packingObject.openStatObject(vdosFile);

	// compress to set packing fraction using FIRE, pressure relaxation
	cout << "	** jamming protocol with Ftol = " << Ftol << ", Ptol = " << Ptol << endl;
	packingObject.enthalpyMin(deltaPhi0, Ftol, Ptol);


	cout << "	** computing VDOS, printing to " << vdosFile << endl;
	packingObject.vdos();

	return 0;
}
























