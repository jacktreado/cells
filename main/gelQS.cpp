/*

	Main file to take an initial dense state ("input.dat")
	and quasistatically, isotropically extend boundaries
	by shrinking particles. Attraction should be non zero

	* * * Force parameters based on gel.cpp example file * * * 

	Input parameters:
		-- NCELLS: 				number of cells
		-- asphericity: 		p^2/4*PI*a, deformability after relaxation
		-- a: 					attraction parameter, defines strength & size of attractive shell in units of del
		-- seed: 				initial seed for the simulation

	Files to write to
		-- positionFile: 	configuration during packing simulation
		-- energyFile:		particle energies during packing simulation

	** In later versions, will add option to input desired input file

	System setup is determined by input file

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
const int NT 					= 5e7; 			// number of time steps
const double T0 				= 1e-4;			// temperature scale
const double sizeRatio 			= 1.4;			// size ratio between large and small particles
const int NPRINT 				= 500;			// number of time steps between prints
const double timeStepMag 		= 0.005;			// time step in MD unit
const double deltaPhi 			= 0.001;		// packing fraction step
const double kineticTol 		= 1e-12;		// kinetic energy tolerance
const double pressureTol 		= 1e-6;			// potential energy tolerance
const double initialPhi 		= 0.6;			// initial packing fraction
const double repPhi 			= 0.79;			// packing fraction of repulsive cells before gelation
const double initialCalA 		= 1.05;			// cal A at end of compression sim
const double deltaA 			= 0.001;		// stepping in a during attraction ramp
const double gelPhi 			= 0.02;			// final phi of gel phase

// force parameters
const double kl 			= 0.2;			// perimeter force constant
const double ka 			= 1.0;			// area force constant
const double gam 			= 0.0;			// surface tension force constant
const double kb 			= 0.0;			// bending energy constant
const double kint 			= 1.0;			// interaction energy constant
const double del 			= 1.0;			// width of vertices (WHEN = 1, ITS A VERTEX FORCE!)
const double aInitial 		= 0.0;			// attraction parameter to start

// int main
int main(int argc, char const *argv[])
{
	// local variables
	int NCELLS, NV, seed, plotIt;
	double asphericity, a;

	// inputs from command line
	string NCELLS_str 			= argv[1];
	string NV_str 				= argv[2];
	string asphericity_str 		= argv[3];
	string a_str 				= argv[4];
	string seed_str				= argv[5];
	string positionFile			= argv[6];
	string energyFile 			= argv[7];

	// load strings into sstream
	stringstream NCELLSss(NCELLS_str);
	stringstream NVss(NV_str);
	stringstream asphericityss(asphericity_str);
	stringstream ass(a_str);
	stringstream seedss(seed_str);

	// parse values from strings
	NCELLSss 		>> NCELLS;
	NVss 			>> NV;
	asphericityss 	>> asphericity;
	ass	 			>> a;
	seedss 			>> seed;

	// box size
	double L = 10.0*NCELLS;

	// instantiate object
	cout << "	** Instantiating object" << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,L,seed);

	// print initial greeting
	cout << "	** BEGINNING main code for gel simulation, initial asphericity = " << asphericity << endl;

	// set values in setting up packing
	cout << "	** Initializing particle quantities, particles are initially regular polygons" << endl;
	packingObject.initializeBidisperse(NV,sizeRatio);
	packingObject.initializeForceConstants(kl,ka,gam,kb,kint);
	packingObject.initializeInteractionParams(del,aInitial);

	// set time step
	packingObject.setdt(timeStepMag);

	// initialize particle positions by BIDISPERSE DISK PACKING
	double initPhiTarget = 0.6;
	cout << "	** Initializing particle positions as disk packing to packing fraction = " << initPhiTarget << endl;
	packingObject.bidisperseDisks(sizeRatio,initPhiTarget);

	// open output files
	cout << "	** Opening print objects" << endl;
	packingObject.openPackingObject(positionFile);
	packingObject.openEnergyObject(energyFile);




	/****************

		Simulation

	*****************/

	// REset time step (in case initial conditions changed value)
	packingObject.setdt(timeStepMag);


	// print simulation header for initial compression
	cout << endl << endl;
	cout << "================================================" << endl << endl;
	cout << "	Beginning gel algorithm using velocity-Verlet and FIRE" << endl;
	cout << " 		* NCELLS 				= " 	<< NCELLS << " and smaller cells have NV = " << NV << " vertices" << endl;
	cout << "		* NPRINT 				= " 	<< NPRINT << endl;
	cout << "		* target asphericity 	= " 	<< asphericity << endl;
	cout << " 		* timeScale 			= " 	<< packingObject.timeScale() << endl;
	cout << "		* timeStepMag 			= " 	<< timeStepMag << endl;
	cout << "		* delta phi 			= "	 	<< deltaPhi << endl;
	cout << "================================================" << endl;
	cout << endl << endl;

	// run simulation 
	int frameCount = 0;
	cout << "	** Compressing to a target phi = " << repPhi << endl;
	packingObject.compressToTarget(deltaPhi,repPhi,initialCalA,kineticTol,pressureTol,1,frameCount);

	// relax shape and attraction
	cout << "	** Relaxing attraction and shape" << endl;
	packingObject.attractionRamp(a, deltaA, asphericity, 1, frameCount);

	// run gel simulation
	cout << "	** Performing isostatic gelation sim" << endl;
	packingObject.isoExtensionQS(1, frameCount, gelPhi, deltaPhi);

	// end main successfully
	return 0;
}

