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
const int NPRINT				= 200;

// simulation constants
const double sizeDispersion 	= 0.125;			// size dispersion (std dev of cell sizes)
const double timeStepMag		= 0.05;				// time step scale

// disk constants
const double phiDisk	 		= 0.91;			// initial packing fraction of disks (sets boundary)

// gelation constants
const double phiGel 			= 0.3;			// final packing fraction
const double dphiGel 			= 5e-4;			// change in packing fraction
const double aGelation			= 0.1;			// attraction parameter during gelation sim

// main function
int main()
{
	// local variables

	// output files
	string posFile = "pos.test";
	string enFile = "en.test";
	string cmFile = "cm.test";

	// system details
	int NCELLS 		= 48; 
	int seed 		= 3;
	double Ltmp 	= 10.0;

	// vector of radii
	vector<double> radii(NCELLS,0.0);

	// instantiate object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed); 	// NOTE: NEED TO MAKE NEW CONSTRUCTOR, EVERYTHING ELSE DONE IN initializeGel AND regularPolygon FUNCTIONS

	// open position output file
	packingObject.openPackingObject(posFile);
	packingObject.openEnergyObject(enFile);
	packingObject.openStatObject(cmFile);

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initializing gel at phiDisk = " << phiDisk << " using SP model" << endl;
	packingObject.stickySPTriangularLattice(radii, phiDisk);

	// -- decrease phi quasistatically
	cout << "	** Running gel extension simulation with gelRate = " << dphiGel << ", phiGel = " << phiGel << endl;
	packingObject.stickySPGelationRate(radii, phiGel, dphiGel, aGelation, timeStepMag);

	return 0;
}
























