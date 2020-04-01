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
const int NPRINT				= 2000;

// simulation constants
const double sizeDispersion 	= 0.1;			// size dispersion (std dev of cell sizes)
const double timeStepMag 		= 0.005;		// time step in MD units (zeta * lenscale / forcescale)

// disk constants
const double phiDisk	 		= 0.7;			// initial packing fraction of disks

// compression constants
const double phiTarget			= 1.0;			// cell packing fraction (regardless of final pressure)
const double deltaPhi			= 0.001;		// compression step size

// force parameters
const double kl 			= 1.0;				// perimeter force constant
const double ka 			= 2.0;				// area force constant
const double gam 			= 0.0;				// surface tension force constant
const double kb 			= 0.0;				// bending energy constant
const double kint 			= 0.5;				// interaction energy constant
const double del 			= 1.0;				// width of vertices in units of l0, vertex sep on regular polygon
const double aInitial 		= 0.0;				// attraction parameter to start

// deformability
const double calA0 			= 1.04;				// ratio of preferred perimeter^2 to preferred area

// tolerances
const double Ftol 			= 1e-8;			// force tolerance (for FIRE min)
const double Ktol 			= 1e-8; 			// kinetic energy tolerance
const double Ptol 			= 1e-6;				// pressure tolerance

// main function
int main()
{
	// local variables

	// output files
	string posFile = "pos.test";
	string enFile = "en.test";
	string jamFile = "jam.test";

	// system details
	int NCELLS 		= 12;
	int NV			= 16;
	int seed 		= 5;
	double Ltmp 	= 1.0;

	// instantiate object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed); 	// NOTE: NEED TO MAKE NEW CONSTRUCTOR, EVERYTHING ELSE DONE IN initializeGel AND regularPolygon FUNCTIONS

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initializing gel at phiDisk = " << phiDisk << " using SP model" << endl;
	packingObject.initializeGel(NV, phiDisk, sizeDispersion, del, ka);

	// set deformability, force values
	packingObject.gelForceVals(calA0,kl,ka,gam,kb,kint,del,aInitial);

	// update time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// open position output file
	packingObject.openJamObject(jamFile);

	// open energy and position file
	packingObject.openPackingObject(posFile);
	packingObject.openEnergyObject(enFile);

	// compress to set packing fraction using FIRE, pressure relaxation
	cout << "	** BEGIN jamming protocol " << endl;
	packingObject.findJamming(deltaPhi, Ktol, Ftol, Ptol);

	return 0;

	// get packing fraction, test to see if we should keep compressing
	double phiJ = packingObject.packingFraction();
	double phiTmp, phiTargetTmp, deltaPhiTmp;

	// check vs phiTarget
	if (phiJ < phiTarget){
		// set first phi targets to be slightly greater than jamming
		phiTargetTmp = phiJ + 1e-6;
		deltaPhiTmp = 1e-8;
		cout << "	** QS compresison protocol to first phiTarget = " << phiTargetTmp << endl;
		packingObject.qsIsoCompression(phiTargetTmp,deltaPhiTmp, Ftol, Ktol);

		// set next phi targets to be consecutively larger
		phiTargetTmp = phiJ + 1e-4;
		deltaPhiTmp = 1e-6;
		cout << "	** QS compresison protocol to second phiTarget = " << phiTargetTmp << endl;
		packingObject.qsIsoCompression(phiTargetTmp,deltaPhiTmp, Ftol, Ktol);

		phiTargetTmp = phiJ + 1e-2;
		deltaPhiTmp = 1e-4;
		cout << "	** QS compresison protocol to third phiTarget = " << phiTargetTmp << endl;
		packingObject.qsIsoCompression(phiTargetTmp,deltaPhiTmp, Ftol, Ktol);

		// check if still under confluent phiTarget
		phiTmp = packingObject.packingFraction();
		if (phiTmp < phiTarget){
			phiTargetTmp = phiTarget;
			deltaPhiTmp = 1e-3;
			cout << "	** QS compresison protocol to final phiTarget = " << phiTargetTmp << endl;
			packingObject.qsIsoCompression(phiTargetTmp,deltaPhiTmp, Ftol, Ktol);
		}
	}

	// Print final confluent config to jammed file
	cout << "	** Printing final confluent state at phi = " << phiTmp << " to jam file" << endl;
	packingObject.printJammedConfig();

	cout << "	** FINISHED COMPRESSING ABOVE JAMMING, ENDING MAIN FILE" << endl;
	return 0;
}
























