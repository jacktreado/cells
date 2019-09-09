/*

	Test .cpp file to test DP extensional
	simulation

*/

// include file
#include "cellPacking2D.h"
#include "deformableParticles2D.h"

// use std name space
using namespace std;

// define PI
const double PI = 4.0*atan(1);

// simulation constants
const int NT 					= 1e6;			// number of time steps
const int NPRINT 				= 500;			// number of time steps between prints
const double T0 				= 0.0;			// temperature scale
const double timeStepMag 		= 0.01;		// time step in MD units
const double dampingParameter 	= 0.9;			// damping parameter
const double phiTarget 			= 0.6;			// target packing fraction
const double dphi				= 0.001;		// decrease in packing fraction increment

// force parameters
const double kl 			= 1.0;			// perimeter force constant
const double ka 			= 1.0;			// area force constant
const double gam 			= 0.0;			// surface tension force constant
const double kb 			= 0.0;			// bending energy constant
const double kint 			= 5.0;			// interaction energy constant
const double del 			= 0.1;			// width of circulo lines
const double C 				= 0.01;			// attraction parameter (strength)
const double l 				= 2.0;			// attraction parameter (distance)

int main(){
	// file strings
	string inputStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/gel/input.dat";
	string positionsStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/gel/gelPositions.dat";
	string energyStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/gel/gelEnergy.dat";

	// input simulation params
	double asphericity 	= 1.09;

	// print initial greeting
	cout << "	** BEGINNING main code for gel simulation, initial asphericity = " << asphericity << endl;

	// open input object
	cout << "	** Opening file object from string " << inputStr << endl;
	ifstream fileInputObject(inputStr.c_str());
	if (!fileInputObject.is_open()){
		cout << "	** ERROR: could not open file " << inputStr << ", ending." << endl;
		return 1;
	}

	// instantiate object
	cout << "	** Instantiating object for gel simulation" << endl;
	cellPacking2D gelObject(fileInputObject, asphericity, 1);

	// get number of particles
	int NCELLS = gelObject.getNCELLS();

	// set values in setting up packing
	cout << "	** Initializing particle quantities" << endl;
	gelObject.initializeForceConstants(kl,ka,gam,kb,kint);
	gelObject.initializeInteractionParams(del,C,l);

	cout << "	** Initializing particle velocities with temperature T0 = " << T0 << endl;
	gelObject.initializeVelocities(0.0);

	// open print objects
	gelObject.openPackingObject(positionsStr);
	gelObject.openEnergyObject(energyStr);
	// gelObject.printSystemPositions(0);
	// return 0;

	/****************

		Simulation

	*****************/

	// get initial phi
	// double newPhi = 0.97;
	// gelObject.setPackingFraction(newPhi);
	double phi = gelObject.packingFraction();
	gelObject.setPackingFraction(phi);

	// set time step
	gelObject.setdt(0.1*timeStepMag);

	// print simulation header
	cout << endl << endl;
	cout << "================================================" << endl << endl;
	cout << "	Beginning gelation algorithm using velocity-Verlet and damping" << endl;
	cout << " 		* NCELLS = " << NCELLS << endl;
	cout << "		* NPRINT = " << NPRINT << endl;
	cout << "		* asphericity = " << asphericity << endl;
	cout << " 		* timeScale = " << gelObject.timeScale() << endl;
	cout << "		* timeStepMag = " << timeStepMag << endl;
	cout << " 		* initial phi = " << phi << endl;
	cout << "		* target phi = " << phiTarget << endl;
	cout << "		* dphi = " << dphi << endl;
	cout << "		* initial potential energy = " << gelObject.totalPotentialEnergy() << endl;
	cout << "================================================" << endl;
	cout << endl << endl;

	// first compress to above jamming, without attraction
	// gelObject.setNPRINT(100);
	// gelObject.setNT(NT);
	// gelObject.msFire(0.005, 0.0025, 0.96, 1e-12, 1e-6);

	// relax IC shapes
	gelObject.setNPRINT(NPRINT);
	gelObject.setNT(NT);
	int frameCount = 0;
	int plotRelax = 0;
	double Ftol = 1e-8;
	int isRelaxed = gelObject.fireRelax(Ftol,dampingParameter,plotRelax,frameCount);
	// return 0;
	// int isRelaxed = 1;

	// run simulation if initial conditions can relax
	if (isRelaxed == 1){
		cout << "	** Initial condition sufficiently relaxed!" << endl;
		cout << "	** Initializing particle velocities with temperature T0 = " << T0 << endl;
		gelObject.initializeVelocities(T0);


		cout << "	** Running quasistatic isotropic extension simulation" << endl;
		gelObject.setNPRINT(NPRINT);
		gelObject.setdt(timeStepMag);
		gelObject.initializeForceConstants(0.0,ka,gam,kb,kint);
		gelObject.initializeInteractionParams(del,C,l);
		gelObject.isoExtensionQS(phiTarget, dphi, dampingParameter);

		// close file input object
		fileInputObject.close();
	}
	else{
		cout << "	** ERROR: initial conditions could not relax in alloted time NT = " << NT << ", ending." << endl;
		return 1;
	}

	// end main successfully
	return 0;
}

