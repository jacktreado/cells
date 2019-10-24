/*

	Main file to take an initial dense state
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
const int NPRINT 				= 5e2;			// number of time steps between prints
const double timeStepMag 		= 0.02;			// time step in MD unit
const double deltaPhi 			= 0.001;		// packing fraction step
const double deltaA 			= 0.001;		// stepping in a during attraction ramp
const double phiDisk 			= 0.75;			// initial phi of SP disks
const double phiGel 			= 0.02;			// final phi of gel phase

// force parameters
const double gam 			= 0.0;			// surface tension force constant
const double kint 			= 1.0;			// interaction energy constant
const double aInitial 		= 0.0;			// attraction parameter to start

// int main
int main(int argc, char const *argv[])
{
	// local variables
	int NCELLS, NV, seed, plotIt;
	double asphericity, a, sizeDisp, phiTarget, gelRate, kl, ka, kb, del;

	// inputs from command line
	string NCELLS_str 			= argv[1];
	string NV_str 				= argv[2];
	string sizeDisp_str 		= argv[3];
	string phiTarget_str 		= argv[4];
	string gelRate_str 			= argv[5];
	string kl_str 				= argv[6];
	string ka_str 				= argv[7];
	string kb_str 				= argv[8];
	string del_str 				= argv[9];
	string asphericity_str 		= argv[10];
	string a_str 				= argv[11];
	string seed_str				= argv[12];
	string positionFile			= argv[13];
	string energyFile 			= argv[14];
	string contactFile 			= argv[15];

	// load strings into sstream
	stringstream NCELLSss(NCELLS_str);
	stringstream NVss(NV_str);
	stringstream sizeDispss(sizeDisp_str);
	stringstream phiTargetss(phiTarget_str);
	stringstream gelRatess(gelRate_str);
	stringstream klss(kl_str);
	stringstream kass(ka_str);
	stringstream kbss(kb_str);
	stringstream delss(del_str);
	stringstream asphericityss(asphericity_str);
	stringstream ass(a_str);
	stringstream seedss(seed_str);

	// parse values from strings
	NCELLSss 		>> NCELLS;
	NVss 			>> NV;
	sizeDispss 		>> sizeDisp;
	phiTargetss 	>> phiTarget;
	gelRatess 		>> gelRate;
	klss 			>> kl;
	kass 			>> ka;
	kbss 			>> kb;
	delss 			>> del;
	asphericityss 	>> asphericity;
	ass	 			>> a;
	seedss 			>> seed;

	// temporary box length; will be modified in initialization
	double Ltmp 	= 1.0;

	// instantiate object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed); 	// NOTE: NEED TO MAKE NEW CONSTRUCTOR, EVERYTHING ELSE DONE IN initializeGel AND regularPolygon FUNCTIONS

	// open position output file
	packingObject.openPackingObject(positionFile);
	packingObject.openEnergyObject(energyFile);
	packingObject.openStatObject(contactFile);

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initializing gel at phiDisk = " << phiDisk << " using SP model" << endl;
	packingObject.initializeGel(NV,phiDisk,sizeDisp,del);

	// set deformability, force values
	packingObject.gelForceVals(asphericity,kl,ka,gam,kb,kint,del,aInitial);

	// update time scale for compresion protocol
	packingObject.setdt(0.2*timeStepMag);

	// compress to set packing fraction using FIRE, pressure relaxation
	cout << "	** QS compresison protocol to phiTarget = " << phiTarget << endl;

	// compress to target phi
	packingObject.qsIsoCompression(phiTarget,deltaPhi);

	// -- ramp attraction
	cout << "	** Ramping attraction to a = " << a << endl;
	packingObject.attractionRamp(a,deltaA);

	// -- decrease phi as if boundary was growing: phi(t) = phi(0)/(1 + a*t)
	cout << "	** Running gel extension simulation with gelRate = " << gelRate << ", phiGel = " << phiGel << endl;
	packingObject.gelRateExtension(phiGel,gelRate,0.5*timeStepMag);

	// end main successfully
	return 0;
}

