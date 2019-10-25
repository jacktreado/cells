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
const int NPRINT				= 500;

// simulation constants
const double timeStepMag 		= 0.0001;			// time step in MD units (zeta * lenscale / forcescale)

// gelation constants
const double phiGel 			= 0.01;		// final packing fraction
const double gelRate 			= 5e-3;			// rate of size decrease (i.e. area loss relative to initial box area)
const double aGelation			= 0.1;			// attraction parameter during gelation sim

// force parameters
const double kl 			= 0.5;				// perimeter force constant
const double ka 			= 0.5;				// area force constant
const double gam 			= 0.0;				// surface tension force constant
const double kb 			= 0.01;				// bending energy constant
const double kint 			= 1.0;				// interaction energy constant
const double del 			= 1.0;				// width of vertices in units of l0, vertex sep on regular polygon
const double aInitial 		= 0.01;				// attraction parameter to start
const double da 			= 0.001;			// attraction increment

// deformability
const double calA0 			= 1.05;				// ratio of preferred perimeter^2 to preferred area

// main function
int main()
{
	// local variables

	// output files
	string posFile = "pos.test";
	string enFile = "en.test";

	// system details
	int NCELLS 		= 2;
	int NV			= 24;
	int seed 		= 1;
	double Ltmp 	= 1.0;

	// instantiate object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed);

	// open position output file
	packingObject.openPackingObject(posFile);
	packingObject.openEnergyObject(enFile);

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initially 2 same-sized cells at slight overlap" << endl;
	packingObject.twoParticleContact(NV);

	// set deformability, force values
	packingObject.gelForceVals(calA0,kl,ka,gam,kb,kint,del,aInitial);

	// -- ramp attraction
	cout << "	** Ramping attraction to a = " << aGelation << endl;
	// packingObject.attractionRamp(aGelation,da);

	// -- decrease phi as if boundary was growing: phi(t) = phi(0)/(1 + a*t)
	cout << "	** Running gel extension simulation with gelRate = " << gelRate << ", phiGel = " << phiGel << endl;
	packingObject.gelRateExtension(phiGel,gelRate,timeStepMag);

	return 0;
}
























