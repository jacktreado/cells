/*

	Main .cpp file to compress NCELLS to jamming
	print jammed configuration, and then
	compress to confluency (phi = 1.03) and record steps in
	between to monitor the pressure required to deform cells 
	into confluence

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
const double timeStepMag 		= 0.005;		// time step in MD unit
const double phiDisk 			= 0.7;			// initial phi of SP disks
const double phiTarget 			= 1.0;			// confluent packing fraction target
const double deltaPhi0 			= 1e-3;			// initial delta phi

// force parameters
const double gam 			= 0.0;			// surface tension force constant
const double kint 			= 0.5;			// interaction energy constant
const double aInitial 		= 0.0;			// attraction parameter to start
const double del 			= 1.0;			// radius of vertices in units of l0

// tolerances
const double Ftol 			= 1e-9;			// force tolerance (for FIRE min)
const double Ktol 			= 1e-11;		// kinetic energy tolerance
const double Ptol 			= 1e-7;			// pressure tolerance

// int main
int main(int argc, char const *argv[])
{
	// local variables
	int NCELLS, NV, seed, plotIt;
	double a, sizeDisp, calA0, kl, ka, kb;

	// inputs from command line
	string NCELLS_str 			= argv[1];
	string NV_str 				= argv[2];
	string sizeDisp_str 		= argv[3];
	string calA0_str 			= argv[4];
	string kl_str 				= argv[5];
	string ka_str 				= argv[6];						
	string kb_str 				= argv[7];
	string seed_str				= argv[8];
	string energyFile 			= argv[9];
	string jammingFile 			= argv[10];

	// load strings into sstream
	stringstream NCELLSss(NCELLS_str);
	stringstream NVss(NV_str);
	stringstream sizeDispss(sizeDisp_str);
	stringstream calA0ss(calA0_str);
	stringstream klss(kl_str);
	stringstream kass(ka_str);
	stringstream kbss(kb_str);
	stringstream seedss(seed_str);

	// parse values from strings
	NCELLSss 		>> NCELLS;
	NVss 			>> NV;
	sizeDispss 		>> sizeDisp;
	calA0ss 		>> calA0;
	klss 			>> kl;
	kass 			>> ka;
	kbss 			>> kb;
	seedss 			>> seed;

	// temporary box length; will be modified in initialization
	double Ltmp 	= 1.0;

	// instantiate main packing object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed); 	// NOTE: NEED TO MAKE NEW CONSTRUCTOR, EVERYTHING ELSE DONE IN initializeGel AND regularPolygon FUNCTIONS

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initializing gel at phiDisk = " << phiDisk << " using SP model" << endl;
	packingObject.initializeGel(NV, phiDisk, sizeDisp, del, ka);

	// set deformability, force values
	packingObject.gelForceVals(calA0,kl,ka,gam,kb,kint,del,aInitial);

	// update time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// open position output file
	packingObject.openJamObject(jammingFile);

	// compress to set packing fraction using FIRE, pressure relaxation
	cout << "	** jamming protocol with Ptol = " << Ptol << " and Ktol = " << Ktol << endl;
	packingObject.findJamming(deltaPhi0, Ktol, Ftol, Ptol);

	// get packing fraction, test to see if we should keep compressing
	double phiJ = packingObject.packingFraction();
	double phiTmp, phiTargetTmp, deltaPhiTmp;

	// open energy and position file
	packingObject.openEnergyObject(energyFile);

	// create vector of logarithmically spaced points
	int nlogpts = 30;
	double p0 = 1e-6;
	double p1 = 1e-2;
	vector<double> dphiLogSpace(nlogpts,0.0);
	double dp = (log10(p1) - log10(p0))/(nlogpts - 1);
	double logp, linp;
	logp = log10(p0);
	dphiLogSpace.at(0) = p0;

	// loop over vector, populate with points
	for (int i=1; i<nlogpts; i++){
		logp = logp + dp;
		linp = pow(10.0,logp);
		dphiLogSpace.at(i) = linp;
	}

	// check vs phiTarget
	if (phiJ < phiTarget){
		// compress packing in nlogpts logarithmically spaced steps to 1e-2 over phiJ
		for (int i=0; i<nlogpts; i++){
			// get current packing fraction
			phiTmp = packingObject.packingFraction();

			// get log spaced delta phi
			deltaPhiTmp = dphiLogSpace.at(i);

			// get new phi target
			phiTargetTmp = phiTmp + deltaPhiTmp;

			// compress to new phiTarget by new deltaPhi
			cout << "	** QS compression protocol i = " << i << " from phi = " << phiTmp << " to phiTarget = " << phiTargetTmp << " by dphi = " << deltaPhiTmp << endl;
			packingObject.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftol, Ktol);
		}

		// after logarithmic compression, print to jam file
		phiTmp = packingObject.packingFraction();
		cout << "	** Printing compressed state at phi = " << phiTmp << " to jam file" << endl;
		packingObject.printJammedConfig();

		// if still not confluent, compress to confluency
		if (phiTmp < phiTarget){
			phiTargetTmp = phiTarget;
			deltaPhiTmp = 1e-3;
			cout << "	** QS compresison protocol to final phiTarget = " << phiTargetTmp << endl;
			packingObject.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftol, Ktol);
		}
	}

	// Print final confluent config to jammed file
	cout << "	** Printing final confluent state at phi = " << packingObject.packingFraction() << " to jam file" << endl;
	packingObject.printJammedConfig();

	cout << "	** FINISHED COMPRESSING ABOVE JAMMING, ENDING MAIN FILE" << endl;
	return 0;
}