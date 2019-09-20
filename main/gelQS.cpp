/*

	Main file to take an initial dense state ("input.dat")
	and quasistatically, isotropically extend boundaries
	by shrinking particles. Attraction should be non zero

	* * * Force parameters based on gel.cpp example file * * * 

	Input parameters:
		-- NCELLS: 				number of cells
		-- asphericity: 		p^2/4*PI*a, deformability after relaxation
		-- C: 					attraction parameter, defines strength of attractive shell in units of del
		-- l: 					attraction parameter, defines distance of attractive shell in units of del
		-- seed: 				initial seed for the simulation

	Files to write to
		-- positionFile: 	configuration during packing simulation
		-- energyFile:		particle energies during packing simulation
		-- statFile: 		packing statistics (final contact network, etc) after packing simulation

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
const double sizeRatio 			= 1.0;			// size ratio between large and small particles
const int NPRINT 				= 500;			// number of time steps between prints
const double timeStepMag 		= 0.02;			// time step in MD unit
const double deltaPhi 			= 0.001;		// packing fraction step
const double deltaCalA			= 0.003;		// asphericity increase step
const double kineticTol 		= 1e-24;		// kinetic energy tolerance
const double potentialTol 		= 1e-16;		// potential energy tolerance
const double initialPhi 		= 0.7;			// initial packing fraction

// force parameters
const double kl 			= 1.0;			// perimeter force constant
const double ka 			= 1.0;			// area force constant
const double gam 			= 0.0;			// surface tension force constant
const double kb 			= 0.0;			// bending energy constant
const double kint 			= 2.0;			// interaction energy constant
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




	/****************

		Simulation

	*****************/

	// set time step
	packingObject.setdt(timeStepMag);

	// set initial packing fraction
	packingObject.setPackingFraction(initialPhi);

	// print simulation header for initial compression
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
	double phiTarget = 0.75;
	plotIt = 0;
	cout << "	** Compressing to a jammed state, do not plot yet" << endl;
	packingObject.jammingFireRamp(deltaPhi,deltaCalA,asphericity,kb,phiTarget,kineticTol,potentialTol,plotIt);

	// ramp to fixed attraction and bending energy
	int frameCount = 0;
	double aTarget = a;
	double da = 0.001;
	cout << "	** QS ramp to attraction, target a = " << aTarget << endl;
	packingObject.attractionRamp(aTarget, da, plotIt, frameCount);

	// run isotropic extension
	frameCount = 0;
	phiTarget = 0.2;
	plotIt = 1;
	double dphi = 0.001;
	cout << "	** QS isotropic extension to form cell gel, starting to plot! aTarget = " << aTarget << endl;
	packingObject.isoExtensionQS(plotIt, frameCount, phiTarget, dphi);

	// end main successfully
	return 0;
}

