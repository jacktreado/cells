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
	string positionsStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/packingFIRE/packingFIREPositions.dat";
	string energyStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/packingFIRE/packingFIREEnergy.dat";

	// input params
	int NCELLS 	= 8;
	int NT 		= 1e6;
	int NPRINT 	= 1e2;
	double L 	= 10.0*NCELLS;

	// instantiate object
	cout << "	** Instantiating object" << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,L,1);

	// number of vertices on small cells
	int NV = 20;

	// Basic cell variables (cells will be identical)
	double kl,ka,gam,kb,kint,del,C,l;
	double asphericity,T0,sizeRatio;

	// force and rest params
	kl 			= 2.0;
	ka 			= 2.0;
	gam 		= 0.0;
	kb 			= 0.0;
	kint 		= 1.0;
	del 		= 0.1;
	C 			= 0.0;
	l 			= 0.0;

	// asphericity
	asphericity = 1.08;

	// temperature
	T0 = 1e-4;

	// size ratio
	sizeRatio = 1.4;

	// set values in setting up packing
	cout << "	** Initializing particle quantities" << endl;
	packingObject.initializeBidisperse(NV,sizeRatio,asphericity);
	packingObject.initializeForceConstants(kl,ka,gam,kb,kint);
	packingObject.initializeInteractionParams(del,C,l);

	// initialize particle positions
	cout << "	** Initializing particle positions on a square lattice" << endl;
	packingObject.squareLattice();

	cout << "	** Initializing particle velocities with temperature T0 = " << T0 << endl;
	packingObject.initializeVelocities(1e-4);

	// open print objects
	packingObject.openPackingObject(positionsStr);
	packingObject.openEnergyObject(energyStr);
	// packingObject.printSystemPositions(0);
	// return 0;


	/****************

		Simulation

	*****************/

	// simulation params
	double timeStepMag 			= 0.1;
	double initialPhi 			= 0.5;
	double deltaPhi0 			= 0.015;
	double deltaPhiJ			= 0.005;
	double phiJGuess 			= 0.99;
	double kineticTol 			= 1e-16;
	double forceTol 			= 1e-6;
	double dampingParameter 	= 1.9;
	double fixedTemperature		= 1e-2;

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

	packingObject.setT(fixedTemperature);

	// run simulation 
	packingObject.msFire(deltaPhi0,deltaPhiJ,phiJGuess,kineticTol,forceTol);
	// packingObject.msDamping(1e-3,kineticTol,forceTol,dampingParameter);

	// end code
	cout << "Finished running simulation in packingFIRE.cpp" << endl;
	return 0;
}



