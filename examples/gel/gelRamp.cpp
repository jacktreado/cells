/*

	Test .cpp file to test DP extensional
	simulation

*/

// include file
#include "cellPacking2D.h"
#include "deformableParticles2D.h"
#include <sstream>

// use std name space
using namespace std;

// define PI
const double PI = 4.0*atan(1);

// simulation constants
const int NT 				= 5e7; 			// number of time steps
const double T0 			= 1e-2;			// temperature scale
const double sizeRatio 		= 1.4;			// size ratio between large and small particles
const double timeStepMag 	= 0.001;		// time step in MD units
const double deltaPhi 		= 0.002;		// packing fraction step
const double deltaCalA		= 0.0005;		// asphericity increase step
const double kineticTol 	= 1e-24;		// kinetic energy tolerance
const double potentialTol 	= 1e-16;		// potential energy tolerance
const double initialPhi 	= 0.86;			// initial packing fraction

// force parameters
const double kl 			= 1.0;			// perimeter force constant
const double ka 			= 1.0;			// area force constant
const double gam 			= 0.0;			// surface tension force constant
const double kb 			= 0.0;			// bending energy constant
const double kint 			= 2.0;			// interaction energy constant
const double del 			= 1.0;			// width of vertices (WHEN = 1, ITS A VERTEX FORCE!)
const double a 				= 0.0;			// attraction parameter

// int main
int main(int argc, char const *argv[])
{
	// local variables
	int NCELLS, NPRINT, NV;
	double asphericity, seed, NPRINTtmp;

	// inputs from command line
	string NPRINT_str 			= argv[1];
	string asphericity_str 		= argv[2];
	string positionFile			= argv[3];
	string energyFile 			= argv[4];

	// input file
	string inputFile = "gelRamp.input";
	ifstream inputFileObj(inputFile.c_str());
	if (!inputFileObj.is_open()){
		cout << "	** ERROR: input file cannot be opened" << endl;
		return 1;
	}

	// load strings into sstream
	stringstream NPRINTss(NPRINT_str);
	stringstream asphericityss(asphericity_str);

	// parse values from strings
	NPRINTss 		>> NPRINTtmp;
	asphericityss 	>> asphericity;

	// set NPRINT to int
	NPRINT = (int)NPRINTtmp;

	// set seed
	seed = 1;

	// instantiate object
	cout << "	** Instantiating object" << endl;
	cellPacking2D packingObject(inputFileObj,asphericity,seed);

	// get values from reading in packing fraction
	NCELLS = packingObject.getNCELLS();
	NV = packingObject.cell(0).getNV();

	// set packing print info
	packingObject.setNT(NT);
	packingObject.setNPRINT(NPRINT);

	// set values in setting up packing
	cout << "	** Initializing particle quantities" << endl;
	packingObject.initializeForceConstants(kl,ka,gam,kb,kint);
	packingObject.initializeInteractionParams(del,a);

	// open print objects
	cout << "	** Opening print objects" << endl;
	packingObject.openPackingObject(positionFile);
	packingObject.openEnergyObject(energyFile);

	/****************

		Simulation

	*****************/

	// set time step
	packingObject.setdt(timeStepMag);

	// set initial packing fraction to be starter - dphi
	packingObject.setPackingFraction(initialPhi);

	// print simulation header
	cout << endl << endl;
	cout << "================================================" << endl << endl;
	cout << "	Beginning compression algorithm using velocity-Verlet and FIRE" << endl;
	cout << " 		* NCELLS 				= " 	<< NCELLS << " and smaller cells have NV = " << NV << " vertices" << endl;
	cout << "		* NT total 				= " 	<< NT << endl;
	cout << "		* NPRINT 				= " 	<< NPRINT << endl;
	cout << "		* starting phi 			= " 	<< initialPhi << endl;
	cout << "		* target asphericity 	= " 	<< asphericity << endl;
	cout << " 		* timeScale 			= " 	<< packingObject.timeScale() << endl;
	cout << "		* timeStepMag 			= " 	<< timeStepMag << endl;
	cout << "		* delta phi 			= "	 	<< deltaPhi << endl;
	cout << "		* delta calA 			= " 	<< deltaCalA << endl;
	cout << "================================================" << endl;
	cout << endl << endl;

	// run simulation 
	// cout << "	** Compressing to a jammed state" << endl;
	// packingObject.jammingFireRamp(deltaPhi,deltaCalA,asphericity,kb,kineticTol,potentialTol);

	// initialize velocities
	cout << "	** Initializing velocities to temp scale " << T0 << endl;
	packingObject.initializeVelocities(T0);

	// ramp to target shape parameter
	double calATarget = 1.5;
	double kbTarget = 0.01;
	double dCalA = 0.001;
	double dkb = 0.001;
	cout << "	** QS ramp to cal A, target cal A= " << calATarget << endl;
	packingObject.shapeRamp(initialPhi,calATarget,dCalA,kbTarget,dkb);

	// ramp to fixed attraction and bending energy
	// double aTarget = 0.5;
	// double da = 0.001;
	// cout << "	** QS ramp to attraction, target a = " << aTarget << endl;
	// packingObject.attractionRamp(aTarget, da);



	// end main successfully
	return 0;
}











