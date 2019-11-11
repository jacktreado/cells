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
const int NT 					= 5e7; 			// number of time steps
const int NPRINT 				= 5e3;			// number of time steps between prints
const double timeStepMag 		= 0.005;		// time step in MD unit
const double phiDisk 			= 0.6;			// initial phi of SP disks

// force parameters
const double gam 			= 0.0;			// surface tension force constant
const double kint 			= 1.0;			// interaction energy constant

/// int main
int main(int argc, char const *argv[])
{
	// local variables
	int NCELLS, NV, seed, plotIt;
	double asphericity, a, sizeDisp, phiTarget, dphi, kl, ka, kb, del;
	double Ltmp = 1.0;

	// inputs from command line
	string NCELLS_str 			= argv[1];
	string NV_str 				= argv[2];
	string sizeDisp_str 		= argv[3];
	string phiTarget_str 		= argv[4];
	string dphi_str 			= argv[5];
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
	stringstream dphiss(dphi_str);
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
	dphiss 			>> dphi;
	klss 			>> kl;
	kass 			>> ka;
	kbss 			>> kb;
	delss 			>> del;
	asphericityss 	>> asphericity;
	ass	 			>> a;
	seedss 			>> seed;

	// instantiate object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed); 	// NOTE: NEED TO MAKE NEW CONSTRUCTOR, EVERYTHING ELSE DONE IN initializeGel AND regularPolygon FUNCTIONS

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initializing particles at phiDisk = " << phiDisk << " using SP model" << endl;
	packingObject.initializeGel(NV, phiDisk, sizeDisp, del);

	// set deformability, force values
	packingObject.gelForceVals(asphericity,kl,ka,gam,kb,kint,del,a);

	// update time scale
	packingObject.setdt(timeStepMag);

		// open position output file
	packingObject.openPackingObject(positionFile);
	packingObject.openEnergyObject(energyFile);
	packingObject.openStatObject(contactFile);

	// compress to set packing fraction using FIRE, pressure relaxation
	cout << "	** QS compresison protocol to phiTarget = " << phiTarget << endl;
	packingObject.qsIsoCompression(phiTarget,dphi);

	cout << endl << endl << endl;
	cout << "	** QS compression protocol completed, ending. " << endl;
	return 0;
}
