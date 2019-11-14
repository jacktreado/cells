/*

	Example file to generate a gel of cells
	from an initial bidisperse sphere packing,
	and to decompress at a fixed rate
	rather than QS

*/

// include files
#include "cellPacking2D.h"
#include "deformableParticles2D.h"
#include <sstream>

// use std name space
using namespace std;

// define PI
const double PI = 4.0*atan(1);

// length paramaters
const int NT 					= 5e6;
const int NPRINT				= 1e3;

// simulation constants
const double sizeDispersion 	= 0.125;		// size dispersion (std dev of cell sizes)
const double timeStepMag		= 0.05;			// time step scale

// disk constants
const double phiDisk	 		= 0.05;		// initial packing fraction of disks (sets boundary)

// force constants
const double a 					= 0.0;			// attractive parameter

// main function
int main()
{
	// local variables

	// output files
	string posFile = "pos.test";
	string enFile = "en.test";
	string cmFile = "cm.test";

	// system details
	int NCELLS 		= 1; 
	int seed 		= 1;
	int NV 			= 50;
	double Ltmp 	= 10.0;
	double v0 		= 1e-2;
	double stddev 	= 0.1*(2.0*PI);
	double vari		= pow(stddev,2.0);
	double Dc 		= 0.5;
	double Dv 		= 0.1;

	// force constants
	double calA0 	= 1.25;
	double kl		= 1.0;
	double ka 		= 1.0;
	double gam 		= 0.0;
	double kb 		= 1e-4;
	double kint 	= 1.0;
	double del 		= 1.0;
	double a 		= 0.0;

	// instantiate object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed); 	// NOTE: NEED TO MAKE NEW CONSTRUCTOR, EVERYTHING ELSE DONE IN initializeGel AND regularPolygon FUNCTIONS

	// open print files
	packingObject.openPackingObject(posFile);
	packingObject.openEnergyObject(enFile);

	// set force values
	packingObject.gelForceVals(calA0,kl,ka,gam,kb,kint,del,a);

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initializing gel at phiDisk = " << phiDisk << " using SP model" << endl;
	packingObject.singleActiveCell(NV,phiDisk,calA0,Dc,vari,v0);
	return 0;
}
























