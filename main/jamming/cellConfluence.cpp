/*

	Main .cpp file to read in existing JAM 
	file, then compress to confluency (nominally phi = 1.0)

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
const double timeStepMag 		= 0.01;			// time step in MD unit
const double dphi 				= 1e-4;			// packing fraction increase
const double T0 				= 1e-8;			// initial velocities for read-in cells

// force parameters
const double ka 				= 1.0;			// area force constant (should be = 1)
const double gam 				= 0.0;			// surface tension force constant
const double kint 				= 1.0;			// interaction energy constant
const double a 					= 0.0;			// attraction parameter 
const double del 				= 1.0;			// radius of vertices in units of l0

// tolerances
const double Ftol 				= 1e-12;		// force tolerance (for FIRE min)

// int main
int main(int argc, char const *argv[])
{
	// input variables
	int NOUTPUTS, seed;
	double calA0, kl, kb, phiTarget;

	// inputs from command line
	string inputFile 			= argv[1];
	string calA0_str 			= argv[2];
	string kl_str 				= argv[3];
	string kb_str 				= argv[4];
	string NOUTPUTS_str 		= argv[5];
	string phiTarget_str 		= argv[6];
	string seed_str				= argv[7];
	string energyFile 			= argv[8];
	string jammingFile 			= argv[9];
	string vdosFile 			= argv[10];

	// load strings into sstream
	stringstream calA0ss(calA0_str);
	stringstream klss(kl_str);
	stringstream kbss(kb_str);
	stringstream NOUTPUTSss(NOUTPUTS_str);
	stringstream phiTargetss(phiTarget_str);
	stringstream seedss(seed_str);

	// parse values from strings
	calA0ss 		>> calA0;
	klss 			>> kl;
	kbss 			>> kb;
	NOUTPUTSss 		>> NOUTPUTS;
	phiTargetss		>> phiTarget;
	seedss 			>> seed;

	// instantiate main packing object
	cout << "	** Reading in from file " << inputFile << endl;
	cellPacking2D packingObject(inputFile,T0,seed);

	// set deformability, force values
	packingObject.forceVals(calA0,ka,kl,gam,kb,kint,del,a);

	// update time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// open position output file
	packingObject.openJamObject(jammingFile);
	packingObject.openEnergyObject(energyFile);
	packingObject.openStatObject(vdosFile);

	// set NT and NPRINT
	packingObject.setNT(NT);
	packingObject.setNPRINT(NPRINT);

	// compress to set packing fraction using FIRE, pressure relaxation
	double Fcheck, Kcheck;
	cout << "	** relaxing system with Ftol = " << Ftol << endl;
	packingObject.fireMinimizeF(Ftol, Fcheck, Kcheck);
	cout << "	** Fcheck = " << Fcheck << endl;
	cout << "	** Kcheck = " << Kcheck << endl;

	// compute pressure
	double Pcheck = 0.5*(packingObject.getSigmaXX() + packingObject.getSigmaYY())/(packingObject.getNCELLS()*packingObject.getL(0)*packingObject.getL(0));
	cout << "	** Pcheck = " << Pcheck << endl << endl;

	// print initial configuration and compute VDOS
	cout << "	** computing VDOS, printing to " << vdosFile << endl << endl;
	packingObject.printJammedConfig();
	packingObject.vdos();
	
	// grow cells to confluence, only print NOUTPUTS frames during compression
	double phi = packingObject.packingFraction();
	double rscale = sqrt((phi + dphi)/phi);

	// determine number of frames to skip based on current packing fraction
	int NSTEPS = round(abs(phiTarget - phi)/dphi);
	int PLOTSKIP = NSTEPS/NOUTPUTS;

	// output to console
	cout << "	** Compressing to confluence over " << NSTEPS << " steps, outputting over " << NOUTPUTS << " frames every " << PLOTSKIP << " steps. " << endl;

	// compress to jamming
	int itmax = 1e3;
	int it = 0;
	while(phi < phiTarget && it < itmax){
		// iterator
		it++;

		// decrease by packing fraction
		cout << "	** Compression protocol it = " << it << " from phi = " << phi << " to phi + dphi = " << phi + dphi << " by dphi = " << dphi << endl;
		packingObject.scaleLengths(rscale);

		// minimize
		packingObject.fireMinimizeF(Ftol, Fcheck, Kcheck);

		// updated packing fraction
		phi = packingObject.packingFraction();
		packingObject.updatePackingFraction();

		// updated pressure/cell
		Pcheck = 0.5*(packingObject.getSigmaXX() + packingObject.getSigmaYY())/(packingObject.getNCELLS()*packingObject.getL(0)*packingObject.getL(0));

		// print vdos and config
		if (it % PLOTSKIP == 0){
			cout << "	** at it = " << it << ", outputting vdos and config to files..." << endl;
			packingObject.vdos();
			packingObject.printJammedConfig();
		}
	}

	cout << "	** FINISHED COMPRESSING ABOVE JAMMING, ENDING MAIN FILE" << endl;
	return 0;
}