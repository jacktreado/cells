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
const int NT 					= 1e7;
const int NPRINT				= 50;

// simulation constants
const double sizeDispersion 	= 0.125;		// size dispersion (std dev of cell sizes)
const double timeStepMag		= 0.05;			// time step scale

// disk constants
const double phiDisk	 		= 0.6;			// initial packing fraction of disks (sets boundary)

// force constants
const double a 					= 0.05;			// attractive parameter

// main function
int main()
{
	// local variables

	// output files
	string posFile = "pos.test";
	string enFile = "en.test";
	string cmFile = "cm.test";

	// system details
	int NCELLS 		= 128; 
	int seed 		= 3;
	int NV 			= 16;
	double Ltmp 	= 10.0;
	double Lscale 	= 2.0;
	double T0 		= 1e-2;

	// vector of radii
	vector<double> radii(NCELLS,0.0);

	// instantiate object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed); 	// NOTE: NEED TO MAKE NEW CONSTRUCTOR, EVERYTHING ELSE DONE IN initializeGel AND regularPolygon FUNCTIONS

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initializing gel at phiDisk = " << phiDisk << " using SP model" << endl;
	packingObject.initializeActiveStickySP(radii,NV,phiDisk,sizeDispersion,Lscale);

	// set attraction of each particle
	for (int ci=0; ci<NCELLS; ci++)
		packingObject.cell(ci).seta(a);

	// open position output file
	packingObject.openPackingObject(posFile);
	packingObject.openEnergyObject(enFile);
	packingObject.openStatObject(cmFile);

	// run sticky SP NVE
	cout << "	** Running NVE using sticky SP model" << endl;
	packingObject.spActivePipeNVE(radii, T0);

	return 0;
}
























