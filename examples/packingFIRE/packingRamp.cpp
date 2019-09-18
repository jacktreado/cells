/*

	Test .cpp file to test DP packing

*/

// include file
#include "cellPacking2D.h"
#include "deformableParticles2D.h"

// use std name space
using namespace std;

// define PI
const double PI = 4.0*atan(1);

int main(){
	// local main variables
	int i,t,d;
	double postmp,veltmp,anew;
	double dx,dy,dnorm,dscale;
	string positionsStr = "examples/packingFIRE/packingRampPositions.test";
	string energyStr = "examples/packingFIRE/packingRampEnergy.test";

	// input params
	int NCELLS 	= 12;
	int NT 		= 1e7;
	int NPRINT 	= 5e2;
	double L 	= 10.0*NCELLS;

	// instantiate object
	cout << "	** Instantiating object" << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,L,1);

	// number of vertices on small cells
	int NV = 12;

	// Basic cell variables (cells will be identical)
	double kl,ka,gam,kb,kint,del,C,l;
	double asphericity,T0,sizeRatio;

	// force and rest params
	kl 			= 2.0;
	ka 			= 2.0;
	gam 		= 0.0;
	kb 			= 0.0;
	kint 		= 1.0;
	del 		= 1.0;		// NOTE: NEEDS TO BE = 1.0 IF DOING VERTEX-VERTEX FORCE ONLY
	a 			= 0.0;

	// asphericity
	asphericity = 1.08;

	// temperature
	T0 = 1e-4;

	// size ratio
	sizeRatio = 1.4;

	// set values in setting up packing
	cout << "	** Initializing particle quantities" << endl;
	packingObject.initializeBidisperse(NV,sizeRatio);
	packingObject.initializeForceConstants(kl,ka,gam,kb,kint);
	packingObject.initializeInteractionParams(del,a);

	// initialize particle positions
	cout << "	** Initializing particle positions on a square lattice" << endl;
	packingObject.squareLattice();

	cout << "	** Initializing particle velocities with temperature T0 = " << T0 << endl;
	packingObject.initializeVelocities(T0);

	// open print objects
	packingObject.openPackingObject(positionsStr);
	packingObject.openEnergyObject(energyStr);
	// packingObject.printSystemPositions(0);
	// return 0;


	/****************

		Simulation

	*****************/

	// simulation params
	double timeStepMag 			= 0.01;
	double initialPhi 			= 0.6;
	double deltaPhi 			= 0.002;
	double kineticTol 			= 1e-24;
	double potentialTol 		= 1e-16;
	double fixedTemperature		= 1e-2;	
	double dCalA				= 0.001;

	// set time step
	packingObject.setdt(timeStepMag);

	// set initial packing fraction
	packingObject.setPackingFraction(initialPhi);

	// print simulation header
	cout << endl << endl;
	cout << "================================================" << endl << endl;
	cout << "	Beginning compression algorithm using velocity-Verlet and FIRE" << endl;
	cout << " 		* NCELLS = " << NCELLS << " with NV = " << NV << "vertices each" << endl;
	cout << "		* NT total = " << NT << endl;
	cout << "		* NPRINT = " << NPRINT << endl;
	cout << "		* asphericity = " << asphericity << endl;
	cout << " 		* timeScale = " << packingObject.timeScale() << endl;
	cout << "		* timeStepMag = " << timeStepMag << endl;
	cout << "================================================" << endl;
	cout << endl << endl;

	// perform initial overlap relief
	packingObject.overlapRelief();
	packingObject.initializeVelocities(T0);

	// run simulation 
	// packingObject.msFire(deltaPhi0,deltaPhiJ,phiJGuess,kineticTol,forceTol);
	// packingObject.msDamping(1e-3,kineticTol,forceTol,dampingParameter);

	// initially relax shape using ramp
	// packingObject.shapeRelaxRamp(asphericity,dCalA,kineticTol,potentialTol);

	// compress system until jamming
	packingObject.jammingFireRamp(deltaPhi,dCalA,asphericity,kb,kineticTol,potentialTol);

	// end code
	cout << "Finished running simulation in packingRamp.cpp" << endl;
	return 0;
}


















