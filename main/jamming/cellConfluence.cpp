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
const double timeStepMag 		= 0.03;			// time step in MD unit
const double dphi 				= 5e-4;			// initial packing fraction increase
const double dphiFast			= 10*dphi; 		// second packing fraction incremenet, to get to larger pressures faster
const double T0 				= 1e-8;			// initial velocities for read-in cells

// force parameters
const double ka 				= 1.0;			// area force constant (should be = 1)
const double gam 				= 0.0;			// surface tension force constant
const double kint 				= 0.5;			// interaction energy constant
const double a 					= 0.0;			// attraction parameter 
const double del 				= 1.0;			// radius of vertices in units of l0

// tolerances
const double Ftol 				= 1e-12;		// force tolerance (for FIRE min)
const double Ptol 				= 1e-8;			// pressure tolerance

// ouputs
const double dlogpTol 			= 0.1;			// log difference between pressures of different frames

// int main
int main(int argc, char const *argv[])
{
	// input variables
	int seed;
	double calA0, kl, kb, pTarget;

	// inputs from command line
	string inputFile 			= argv[1];
	string calA0_str 			= argv[2];
	string kl_str 				= argv[3];
	string kb_str 				= argv[4];
	string pTarget_str 			= argv[5];
	string seed_str				= argv[6];
	string energyFile 			= argv[7];
	string jammingFile 			= argv[8];
	string vdosFile 			= argv[9];

	// load strings into sstream
	stringstream calA0ss(calA0_str);
	stringstream klss(kl_str);
	stringstream kbss(kb_str);
	stringstream pTargetss(pTarget_str);
	stringstream seedss(seed_str);

	// parse values from strings
	calA0ss 		>> calA0;
	klss 			>> kl;
	kbss 			>> kb;
	pTargetss		>> pTarget;
	seedss 			>> seed;

	// instantiate main packing object
	cout << "	** Reading in from file " << inputFile << endl;
	cellPacking2D packingObject(inputFile,T0,seed);

	// set deformability, force values
	packingObject.forceVals(calA0,ka,kl,gam,kb,kint,del,a);

	// update time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// open jamming and vdos output file
	packingObject.openJamObject(jammingFile);
	packingObject.openStatObject(vdosFile);

	// set NT and NPRINT
	packingObject.setNT(NT);
	packingObject.setNPRINT(NPRINT);

	// compress to set packing fraction using FIRE, pressure relaxation
	double Fcheck, Kcheck;
	cout << "	** relaxing system Enthalpy H = U + PV with Ftol = " << Ftol << endl;
	packingObject.enthalpyMin(dphi, Ftol, Ptol);

	// compute pressure in initial jammed state
	double Pcheck = 0.5*(packingObject.getSigmaXX() + packingObject.getSigmaYY())/(packingObject.getNCELLS()*packingObject.getL(0)*packingObject.getL(0));
	cout << "	** Pcheck = " << Pcheck << endl << endl;

	// print initial configuration and compute VDOS
	cout << "	** computing VDOS, printing to " << vdosFile << endl << endl;
	packingObject.vdos();

	// IDEA: compress by small packing fraction steps, 
	// only output when difference between pressures exceeds increasing scale

	// open energy output file
	packingObject.openEnergyObject(energyFile);

	// print first line of energy relaxation (to match with first jamming frame)
	packingObject.printSystemEnergy(1);

	// determine thresholds for output based on difference in pressure between frames
	double pCurrent = Pcheck;
	double pLast = pCurrent;
	double dlogp = 0.0;
	int NFRAMES = ceil((log10(pTarget) - log10(pCurrent))/dlogpTol);

	// output to console
	cout << "	** Compressing to confluence, printing " << NFRAMES << " frames, starting from pressure " << pCurrent << " to final pressure " << pTarget << endl;

	// compute degree to which size should be increased
	double phi = packingObject.packingFraction();
	double rscale = sqrt((phi + dphi)/phi);

	// compress to target pressure
	int itmax = 1e3;
	int it = 0;
	while(Pcheck < pTarget && it < itmax){
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
		pCurrent = Pcheck;

		// compute difference in pressure since last output
		dlogp = log10(pCurrent) - log10(pLast);

		// print vdos and config when difference in pressures is large enough
		if (dlogp > dlogpTol){
			cout << "	** at it = " << it << ", dlogp = " << dlogp << " which is > " << dlogpTol << ", so outputting vdos and config to files..." << endl;
			packingObject.vdos();
			packingObject.printJammedConfig();

			// reset pLast
			pLast = pCurrent;

			// if pressure above fast pressure, increase compression step size
			if (pCurrent > 1e-2*pTarget)
				rscale = sqrt((phi + dphiFast)/phi);
		}
	}

	// print final configuration
	cout << "	** at it = " << it << ", outputting FINAL vdos, config and energy to files..." << endl;
	packingObject.vdos();
	packingObject.printJammedConfig();
	packingObject.printSystemEnergy(1);

	cout << "	** FINISHED COMPRESSING ABOVE JAMMING, ENDING MAIN FILE" << endl;
	return 0;
}