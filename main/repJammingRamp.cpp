/*

	Main file to compress initially dilute DPM particles to
	a mechanically stable (MS) state. Will ramp asphericity from
	minimum asphericity (of regular polygons) to target asphericity

	* * * no attraction (C,l = 0) or bending (kb = 0) * * *

	Input parameters:
		-- NCELLS: 			number of cells
		-- NPRINT: 			number of time steps between printing
		-- NV: 				number of vertices on smaller particles
		-- asphericity: 	p^2/4*PI*a, which defines deformability
		-- seed: 			initial seed for the simulation

	Files to write to
		-- positionFile: 	configuration during packing simulation
		-- energyFile:		particle energies during packing simulation
		-- statFile: 		packing statistics (final contact network, etc) after packing simulation

	Systems are bidisperse, 50:50 mixture with 1.4 size ratio

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
const double T0 			= 1e-4;			// temperature scale
const double sizeRatio 		= 1.0;			// size ratio between large and small particles
const double timeStepMag 	= 0.02;		// time step in MD units
const double initialPhi		= 0.7;			// initial packing fraction
const double deltaPhi 		= 0.001;		// packing fraction step
const double deltaCalA		= 0.001;		// asphericity increase step
const double kineticTol 	= 1e-24;		// kinetic energy tolerance
const double potentialTol 	= 1e-16;		// potential energy tolerance

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
	string NCELLS_str 			= argv[1];
	string NPRINT_str 			= argv[2];
	string NV_str 				= argv[3];
	string asphericity_str 		= argv[4];
	string seed_str				= argv[5];
	string positionFile			= argv[6];
	string energyFile 			= argv[7];
	string statFile 			= argv[8];

	// load strings into sstream
	stringstream NCELLSss(NCELLS_str);
	stringstream NPRINTss(NPRINT_str);
	stringstream NVss(NV_str);
	stringstream asphericityss(asphericity_str);
	stringstream seedss(seed_str);

	// parse values from strings
	NCELLSss 		>> NCELLS;
	NPRINTss 		>> NPRINTtmp;
	NVss 			>> NV;
	asphericityss 	>> asphericity;
	seedss 			>> seed;

	// set NPRINT to int
	NPRINT = (int)NPRINTtmp;

	// box size
	double L = 10.0*NCELLS;

	// instantiate object
	cout << "	** Instantiating object" << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,L,seed);

	// set values in setting up packing
	cout << "	** Initializing particle quantities, particles are initially regular polygons" << endl;
	packingObject.initializeBidisperse(NV,sizeRatio);
	packingObject.initializeForceConstants(kl,ka,gam,kb,kint);
	packingObject.initializeInteractionParams(del,a);

	// initialize particle positions
	cout << "	** Initializing particle positions on a square lattice" << endl;
	packingObject.squareLattice();

	// initialize particle velocities
	cout << "	** Initializing particle velocities with temperature T0 = " << T0 << endl;
	packingObject.initializeVelocities(T0);

	// open print objects
	cout << "	** Opening print objects" << endl;
	packingObject.openPackingObject(positionFile);
	packingObject.openEnergyObject(energyFile);
	packingObject.openStatObject(statFile);

	/****************

		Simulation

	*****************/

	// set time step
	packingObject.setdt(timeStepMag);

	// set initial packing fraction
	packingObject.setPackingFraction(initialPhi);

	// print simulation header
	cout << endl << endl;
	cout << "================================================" << endl << endl;
	cout << "	Beginning compression algorithm using velocity-Verlet and FIRE" << endl;
	cout << " 		* NCELLS 				= " 	<< NCELLS << " and smaller cells have NV = " << NV << " vertices" << endl;
	cout << "		* NT total 				= " 	<< NT << endl;
	cout << "		* NPRINT 				= " 	<< NPRINT << endl;
	cout << "		* target asphericity 	= " 	<< asphericity << endl;
	cout << " 		* timeScale 			= " 	<< packingObject.timeScale() << endl;
	cout << "		* timeStepMag 			= " 	<< timeStepMag << endl;
	cout << "		* delta phi 			= "	 	<< deltaPhi << endl;
	cout << "		* delta calA 			= " 	<< deltaCalA << endl;
	cout << "================================================" << endl;
	cout << endl << endl;

	// perform initial overlap relief
	cout << "	** Peforming a wee relief of any initial overlaps" << endl;
	packingObject.overlapRelief();

	// run simulation 
	double phiTarget = 1.1;
	cout << "	** Compressing to a jammed state" << endl;
	packingObject.jammingFireRamp(deltaPhi,deltaCalA,asphericity,kb,phiTarget,kineticTol,potentialTol,1);

	// once completed, print stats
	cout << "	** Compression protocol completed, printing stats to file." << endl;
	packingObject.printSystemStats();

	return 0;
}
