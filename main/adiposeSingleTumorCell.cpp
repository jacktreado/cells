/*

	Main file to run background of NON-ATTRACTIVE adipose tissue
	with constantly drift forced tumor cell of a certain deformability

	* * * no attraction (C,l = 0) or bending (kb = 0) * * *

	Input parameters:
		-- NCELLS: 			number of cells (single tumor cell)
		-- NV: 				number of vertices on adipose particles (tumor = 1/2)
		-- tumor calA: 		tumor asphericity parameter
		-- forceScale:		forcing on tumor cell
		-- seed: 			initial seed for the simulation

	Files to write to
		-- positionFile: 	configuration during invasion simulation
		-- energyFile:		particle energies during invasion simulation

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
const int NTUMORCELLS 			= 1;			// only 1 tumor cell
const int NT 					= 1e6; 			// number of time steps
const int NPRINT 				= 5e2;			// number of steps between printing
const double T0 				= 1.0;			// temperature of background
const double adiposeCalA		= 1.01;			// adipose deformability
const double timeStepMag 		= 0.01;			// time step in MD units
const double initialPhi			= 0.2;			// initial packing fraction
const double targetPhi			= 0.95;			// target packing fraction
const double dphi 				= 1e-5;			// packing fraction step
const double adiposeDamping		= 0.025;		// damping on adipose tissue

// int main
int main(int argc, char const *argv[])
{
	// local variables
	int NCELLS, adiposeNV, tumorNV;
	double tumorCalA, forceScale, seed;

	// inputs from command line
	string NCELLS_str 			= argv[1];
	string adiposeNV_str 		= argv[2];
	string tumorCalA_str 		= argv[3];
	string forceScale_str 		= argv[4];
	string seed_str				= argv[5];
	string positionFile			= argv[6];
	string energyFile 			= argv[7];

	// load strings into sstream
	stringstream NCELLSss(NCELLS_str);
	stringstream adiposeNVss(adiposeNV_str);
	stringstream tumorCalAss(tumorCalA_str);
	stringstream forceScaless(forceScale_str);
	stringstream seedss(seed_str);

	// parse values from strings
	NCELLSss 		>> NCELLS;
	adiposeNVss 	>> adiposeNV;
	tumorCalAss		>> tumorCalA;
	forceScaless 	>> forceScale;
	seedss 			>> seed;

	// box size
	double L = 10.0*NCELLS;

	// tumor NV
	tumorNV = round(0.75*adiposeNV);

	// instantiate object
	cout << "	** Instantiating object" << endl;
	cellPacking2D tissueObject(NCELLS,NTUMORCELLS,tumorNV,adiposeNV,tumorCalA,adiposeCalA,seed);

	// initialize particle velocities
	cout << "	** Initializing particle velocities with temperature T0 = " << T0 << endl;
	tissueObject.initializeVelocities(T0);

	// open print objects
	cout << "	** Opening print objects" << endl;
	tissueObject.openPackingObject(positionFile);
	tissueObject.openEnergyObject(energyFile);

	/****************

		Simulation

	*****************/

	// setup sim parameters
	tissueObject.setNT(NT);
	tissueObject.setNPRINT(NPRINT);
	tissueObject.setdt(timeStepMag);
	tissueObject.setPackingFraction(initialPhi);

	// perform initial overlap relief
	cout << "	** Performing initial overlap relief" << endl;
	tissueObject.overlapRelief();

	// print simulation header
	cout << endl << endl;
	cout << "================================================" << endl << endl;
	cout << "	Beginning adipose invasion sim at phi = " << initialPhi << endl;
	cout << " 		* NCELLS 				= " 	<< NCELLS << " and adipose cells have NV = " << adiposeNV << " vertices" << endl;
	cout << "		* NT total 				= " 	<< NT << endl;
	cout << "		* NPRINT 				= " 	<< NPRINT << endl;
	cout << "		* tumorCalA 			= " 	<< tumorCalA << endl;
	cout << " 		* timeScale 			= " 	<< tissueObject.timeScale() << endl;
	cout << "		* timeStepMag 			= " 	<< timeStepMag << endl;
	cout << "================================================" << endl;
	cout << endl << endl;

	// run simulation 
	cout << "	** Quasi-static compression to target packing fraction" << endl;
	tissueObject.setNPRINT(5e1);
	tissueObject.rateCompression(targetPhi, dphi, adiposeDamping);

	// cout << "	** Forcing tumor cell through tissue to a jammed state" << endl;
	tissueObject.setNPRINT(NPRINT);
	// tissueObject.tumorForce(NTUMORCELLS, forceScale, adiposeDamping);

	return 0;
}
