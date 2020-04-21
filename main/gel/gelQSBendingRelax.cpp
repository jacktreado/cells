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
const int NT 					= 1e6;			// number of time steps
const int NPRINT 				= 500;			// number of time steps between prints
const double T0 				= 0.0;			// temperature scale
const double timeStepMag 		= 0.0075;		// time step in MD units
const double dampingParameter 	= 1.9;			// damping parameter
const double phiTarget 			= 0.6;			// target packing fraction
const double dphi				= 0.001;		// decrease in packing fraction increment

// force parameters
const double kl 			= 1.0;			// perimeter force constant
const double ka 			= 1.0;			// area force constant
const double gam 			= 0.0;			// surface tension force constant
const double kint 			= 5.0;			// interaction energy constant
const double del 			= 0.1;			// width of circulo lines

// int main
int main(int argc, char const *argv[])
{
	// local variables
	int NCELLS;
	double asphericity, C, l, kb, seed;

	// inputs from command line
	string NCELLS_str 			= argv[1];
	string asphericity_str 		= argv[2];
	string C_str 				= argv[3];
	string l_str 				= argv[4];
	string kb_str 				= argv[5];
	string seed_str				= argv[6];
	string positionFile			= argv[7];

	// load strings into sstream
	stringstream NCELLSss(NCELLS_str);
	stringstream asphericityss(asphericity_str);
	stringstream Css(C_str);
	stringstream lss(l_str);
	stringstream kbss(kb_str);
	stringstream seedss(seed_str);

	// parse values from strings
	NCELLSss 		>> NCELLS;
	asphericityss 	>> asphericity;
	Css	 			>> C;
	lss 			>> l;
	kbss 			>> kb;
	seedss 			>> seed;

	// box size
	double L = 10.0*NCELLS;

	// open input object
	string inputStr = "input.dat";
	cout << "	** Opening file object from string " << inputStr << endl;
	ifstream fileInputObject(inputStr.c_str());
	if (!fileInputObject.is_open()){
		cout << "	** ERROR: could not open file " << inputStr << ", ending." << endl;
		return 1;
	}

	// instantiate object
	cout << "	** Instantiating object for gel simulation" << endl;
	cellPacking2D gelObject(fileInputObject, asphericity, 1);

	// print initial greeting
	cout << "	** BEGINNING main code for gel simulation, initial asphericity = " << asphericity << endl;

	// set values in setting up packing
	cout << "	** Initializing particle quantities * WITH BENDING ENERGY = " << kb << " *" << endl;
	gelObject.initializeForceConstants(kl,ka,gam,kb,kint);
	gelObject.initializeInteractionParams(del,C,l);

	cout << "	** Initializing particle velocities with temperature T0 = " << T0 << endl;
	gelObject.initializeVelocities(0.0);

	// open print objects
	gelObject.openPackingObject(positionFile);

	/****************

		Simulation

	*****************/

	// get initial phi from input.dat
	double phi = gelObject.packingFraction();
	gelObject.setPackingFraction(phi);

	// set time step
	gelObject.setdt(0.5*timeStepMag);

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

	// relax IC shapes
	gelObject.setNPRINT(50*NPRINT);
	gelObject.setNT(NT);
	int frameCount = 0;
	int plotRelax = 0;
	double Ftol = 1e-8;
	int isRelaxed = gelObject.fireRelax(Ftol,dampingParameter,plotRelax,frameCount);

	// run simulation if initial conditions can relax
	if (isRelaxed == 1){
		cout << "	** Initial condition sufficiently relaxed!" << endl;
		cout << "	** Initializing particle velocities with temperature T0 = " << T0 << endl;
		gelObject.initializeVelocities(T0);

		cout << "	** Running quasistatic isotropic extension simulation, now kb = " << kb << endl;
		gelObject.setNPRINT(NPRINT);
		gelObject.setdt(timeStepMag);
		gelObject.initializeForceConstants(kl,ka,gam,kb,kint);
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

