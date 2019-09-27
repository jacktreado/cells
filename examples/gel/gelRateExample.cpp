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

// simulation constants
const double sizeRatio 			= 1.4;			// size ratio between large and small particles
const double timeStepMag 		= 0.005;			// time step in MD units
const double Ktol				= 1e-12;		// kinetic tolerance for compression protocol
const double Ptol 				= 1e-4;			// pressure tolerance for compression protocol

// force parameters
const double kl 			= 1.0;			// perimeter force constant
const double ka 			= 1.0;			// area force constant
const double gam 			= 0.0;			// surface tension force constant
const double kb 			= 1.0;			// bending energy constant
const double kint 			= 1.0;			// interaction energy constant
const double del 			= 1.0;			// width of vertices (WHEN = 1, ITS A VERTEX FORCE!)
const double aInitial 		= 0.0;			// attraction parameter to start


// main function
int main()
{
	// output files
	string posFile = "pos.test";
	string enFile = "en.test";

	// system details
	int NCELLS 		= 12;
	int NV			= 24;
	int seed 		= 10;
	int initNT 		= 5e7;
	int initNPRINT 	= 1e3;
	double L 		= 10.0*NCELLS;

	// instantiate object
	cout << "	** Instantiating object" << endl;
	cellPacking2D packingObject(NCELLS,initNT,initNPRINT,L,seed);

	// set values in setting up packing
	cout << "	** Initializing particle quantities, particles are initially regular polygons" << endl;
	packingObject.initializeBidisperse(NV,sizeRatio);
	packingObject.initializeForceConstants(kl,ka,gam,kb,kint);
	packingObject.initializeInteractionParams(del,aInitial);

	// set time step
	packingObject.setdt(timeStepMag);

	// initialize particle positions by BIDISPERSE DISK PACKING
	double initPhiTarget = 0.2;
	cout << "	** Initializing particle positions as disk packing to packing fraction = " << initPhiTarget << endl;
	packingObject.bidisperseDisks(sizeRatio,initPhiTarget);
	packingObject.cell(0).vertexPerturbation(0.1);

	// open output files
	cout << "	** Opening print objects" << endl;
	packingObject.openPackingObject(posFile);
	packingObject.openEnergyObject(enFile);

	// print positions
	// packingObject.printSystemPositions(0);


	/****************

		Simulation

	*****************/

	// REset time step (in case initial conditions changed value)
	packingObject.setdt(timeStepMag);

	// cout << packingObject.getdt() << endl;
	// cout << packingObject.timeScale() << endl;
	// cout << packingObject.meanAsphericity() << endl;
	// cout << packingObject.cell(0).getkint() << endl;
	// return 0;

	// run simulation 
	double deltaPhi = 0.001;
	double phiTarget = 0.79;
	double asphericityTarget = 1.05;
	int frameCount = 0;
	cout << "	** Compressing to a target phi = " << phiTarget << endl;
	packingObject.compressToTarget(deltaPhi,phiTarget,asphericityTarget,Ktol,Ptol,1,frameCount);

	// relax shape and attraction
	double aTarget = 0.01;
	double dAttraction = 0.001;
	asphericityTarget = 1.1;
	cout << "	** Relaxing attraction and shape" << endl;
	packingObject.attractionRamp(aTarget, dAttraction, asphericityTarget, 1, frameCount);

	// run gel simulation
	phiTarget = 0.2;
	deltaPhi = 0.001;
	packingObject.isoExtensionQS(1, frameCount, phiTarget, deltaPhi);

	return 0;
}