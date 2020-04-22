/*

	Main .cpp file to compress NCELLS adipose cells
	to target phi, and then input a motile tumor cell which
	tries to migrate in the system

	Simulation:

		-- NCELLS initially packed to target packing fraction
			** uniform size, with shape parameter adiposeCalA0

		-- particle 1 selected as tumor cell
			** shrunk to half area
			** given new shape parameter tumorCalA0

		-- given initial motility, allowed to diffuse within tissue with Brownian dynamics
			** tumor cell active brownian cell, tissue athermal brownian particles

		-- compress initial cells without adhesion, turn on after particle 1 turned into a tumor cell
*/

// include files
#include "cellPacking2D.h"
#include "deformableParticles2D.h"
#include <sstream>

// use std name space
using namespace std;

// define PI
const double PI = 4.0*atan(1);

// simulation constants
const int NT 					= 1e7; 			// number of time steps
const int NPRINT 				= 2e3;			// number of time steps between prints
const double timeStepMag 		= 0.05;			// time step in MD unit
const double phiTarget 			= 0.95;			// target phi of initial packing
const double dphi 				= 1e-3;			// delta phi during packing phase
const double sizeRatio 			= 1.4;			// ratio between small and large particles
const double sizeFraction		= 0.5;			// fraction of small particles

// force parameters
const double ka 				= 1.0;			// area force constant (should be = 1)
const double gam 				= 0.0;			// surface tension force constant
const double kint 				= 0.1;			// interaction energy constant
const double del 				= 1.0;			// radius of vertices in units of l0

// tolerances
const double Ftol 				= 1e-6;			// force tolerance (for FIRE min)
const double Ktol 				= 1e-8;		// kinetic energy tolerance

// int main
int main(int argc, char const *argv[])
{
	// local variables
	int NCELLS, NV, seed, plotIt;
	double a, sizeDisp, calA0, kl, ka, kb;

	// inputs from command line
	string NCELLS_str 			= argv[1];
	string NV_str 				= argv[2];
	string tumorCalA0_str 		= argv[3];
	string adiposeCalA0_str 	= argv[3];
	string tumorMotility_str 	= argv[4];
	string kl_str 				= argv[4];
	string kb_str 				= argv[5];
	string seed_str				= argv[6];
	string jammingFile 			= argv[7];

	// load strings into sstream
	stringstream NCELLSss(NCELLS_str);
	stringstream NVss(NV_str);
	stringstream calA0ss(calA0_str);
	stringstream klss(kl_str);
	stringstream kbss(kb_str);
	stringstream seedss(seed_str);

	// parse values from strings
	NCELLSss 		>> NCELLS;
	NVss 			>> NV;
	calA0ss 		>> calA0;
	klss 			>> kl;
	kbss 			>> kb;
	seedss 			>> seed;

	// temporary box length; will be modified in initialization
	double Ltmp 	= 1.0;

	// instantiate main packing object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed);

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initializing gel at phiDisk = " << phiDisk << " using SP model" << endl;
	packingObject.initializeBidisperse(NV, phiDisk, sizeRatio, sizeFraction, del);

	// set deformability, force values
	packingObject.forceVals(calA0,ka,kl,gam,kb,kint,del,a);

	// update time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// open position output file
	packingObject.openJamObject(jammingFile);

	// compress to set packing fraction using FIRE, pressure relaxation
	cout << "	** jamming protocol with Ptol = " << Ptol << " and Ktol = " << Ktol << endl;
	packingObject.findJamming(deltaPhi0, Ktol, Ftol, Ptol);

	cout << "	** FINISHED COMPRESSING TO JAMMING, ENDING MAIN FILE" << endl;
	return 0;
}