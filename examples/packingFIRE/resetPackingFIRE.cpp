/*

	Test .cpp file to test DP packing
	with input from previous step in packing sim

*/

// include file
#include "cellPacking2D.h"
#include "deformableParticles2D.h"

// use std name space
using namespace std;

// define PI
const double PI = 4.0*atan(1);

int main(){
	// local string variables
	string inputStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/packingFIRE/input.dat";
	string positionsStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/packingFIRE/packingFIREPositions.dat";
	string energyStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/packingFIRE/packingFIREEnergy.dat";

	// print stuff
	int NT = 1e6;
	int NPRINT = 100;

	// asphericity
	double asphericity = 1.04;

	// print initial greeting
	cout << "	** BEGINNING main code for jamming simulation, initial asphericity = " << asphericity << endl;

	// open input object
	cout << "	** Opening file object from string " << inputStr << endl;
	ifstream fileInputObject(inputStr.c_str());
	if (!fileInputObject.is_open()){
		cout << "	** ERROR: could not open file " << inputStr << ", ending." << endl;
		return 1;
	}

	// instantiate object
	cout << "	** Instantiating object for jamming simulation" << endl;
	cellPacking2D packingObject(fileInputObject, asphericity);

	// get number of particles
	int NCELLS = packingObject.getNCELLS();

	// set packing print info
	packingObject.setNT(NT);
	packingObject.setNPRINT(NPRINT);

	// Basic cell variables (cells will be identical)
	double kl,ka,gam,kb,kint,del,a;

	// force and rest params
	kl 			= 1.0;
	ka 			= 1.0;
	gam 		= 0.0;
	kb 			= 0.0;
	kint 		= 2.0;
	del 		= 0.1;
	a 			= 0.0;

	// set values in setting up packing
	cout << "	** Initializing particle quantities" << endl;
	packingObject.initializeForceConstants(kl,ka,gam,kb,kint);
	packingObject.initializeInteractionParams(del,a);

	// open print objects
	cout << "	** Opening print objects" << endl;
	packingObject.openPackingObject(positionsStr);
	packingObject.openEnergyObject(energyStr);

	// cout << "	** Printing particle positions" << endl;
	// packingObject.printSystemPositions(0);
	// return 0;

	/****************

		Simulation

	*****************/

	// parameters
	double timeStepMag 		= 0.05;
	double relaxTimeStepMag = 0.005;
	double deltaPhi0 		= 0.001;
	double deltaPhiJ 		= 0.0001;
	double phiJGuess 		= 0.99;
	double kineticTol 		= 1e-16;
	double forceTol 		= 1e-10;

	// get initial phi
	double phi = packingObject.packingFraction();
	packingObject.setPackingFraction(phi);

	// initial relaxation
	packingObject.setdt(relaxTimeStepMag);
	packingObject.fireRelax(forceTol,0);
	cout << "	** Energy sufficiently relaxed, beginning reset compression" << endl;

	// get updated phi, just in case
	phi = packingObject.packingFraction();
	packingObject.setPackingFraction(phi);

	// set time step
	packingObject.setdt(timeStepMag);

	// print simulation header
	cout << endl << endl;
	cout << "================================================" << endl << endl;
	cout << "	Beginning compression algorithmfrom reset using velocity-Verlet and FIRE" << endl;
	cout << " 		* NCELLS = " << NCELLS << endl;
	cout << "		* NT total = " << NT << endl;
	cout << "		* NPRINT = " << NPRINT << endl;
	cout << "		* asphericity = " << asphericity << endl;
	cout << " 		* timeScale = " << packingObject.timeScale() << endl;
	cout << "		* timeStepMag = " << timeStepMag << endl;
	cout << "		* initial phi = " << phi << endl;
	cout << "		* initial U = " << packingObject.totalPotentialEnergy() << endl;
	cout << "================================================" << endl;
	cout << endl << endl;

	// run simulation 
	packingObject.msFire(deltaPhi0,deltaPhiJ,phiJGuess,kineticTol,forceTol);
	// packingObject.msDamping(deltaPhi,kineticTol,forceTol,dampingParameter);

	// end code
	cout << "Finished running simulation in packingFIRE.cpp" << endl;
	return 0;
}