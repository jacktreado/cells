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
const double timeStepMag 		= 0.02;			// time step in MD unit
const double phiDisk 			= 0.5;			// initial phi of SP disks
const double deltaPhi0 			= 5e-4;			// initial delta phi
const double sizeRatio 			= 1.4;			// ratio between small and large particles
const double sizeFraction		= 0.5;			// fraction of small particles

// target packing fraction (confluence)
const double phiTarget 			= 1.0;

// force parameters
const double ka 				= 1.0;			// area force constant (should be = 1)
const double gam 				= 0.0;			// surface tension force constant
const double kint 				= 0.05;			// interaction energy constant
const double a 					= 0.0;			// attraction parameter 
const double del 				= 1.0;			// radius of vertices in units of l0

// tolerances
const double Ftol 				= 1e-10;		// force tolerance (for FIRE min)
const double Ptol 				= 1e-8;			// pressure tolerance

// int main
int main(int argc, char const *argv[])
{
	// local variables
	int NCELLS, NV, seed;
	double calA0, kl, kb;

	// inputs from command line
	string NCELLS_str 			= argv[1];
	string NV_str 				= argv[2];
	string calA0_str 			= argv[3];
	string kl_str 				= argv[4];
	string kb_str 				= argv[5];
	string seed_str				= argv[6];
	string energyFile 			= argv[7];
	string jammingFile 			= argv[8];
	string vdosFile 			= argv[9];

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

	// instantiate main packing object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,1.0,seed);

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
	packingObject.findJamming(deltaPhi0, Ftol, Ptol);

	// get packing fraction, test to see if we should keep compressing
	double phiJ = packingObject.packingFraction();
	double phiTmp, phiTargetTmp, deltaPhiTmp;

	// open energy and vdos file
	packingObject.openEnergyObject(energyFile);
	packingObject.openStatObject(vdosFile);

	// compute initial vdos
	packingObject.vdos();

	// step 1: compress by dphi to dphi*10, repeat for dphi = 1e-8, 1e-7, 1e-6, 1e-5, 1e-4
	// 	-- only save jammed configuration and vdos after compression (5 configs, 5 sets fo evals)
	// 
	// step 2: compress by dphi = 1e-3 to phiTarget, saving every NSKIP configs (based on distance of phiJ to confluence)
	// step 3: save configurationa and vdos of configuration at phi target

	// variables for compression to confluence
	int NSTEPS;
	double phiRange;
	int PRINTSTEPS = 20;
	int NSKIP;


	// check vs phiTarget
	if (phiJ < phiTarget){

		// initial dphi
		deltaPhiTmp = 1e-8;

		// initial jammed packing fraction
		phiTmp = phiJ;

		// print to console
		cout << "	** Beginning compression protocol, first compressing by dphi = 1e-8 from phiJ = " << phiTmp << endl; 

		// compress packing in nlogpts logarithmically spaced steps to 1e-2 over phiJ
		for (int i=0; i<5; i++){
			// get new phi target
			phiTargetTmp = phiTmp + 10*deltaPhiTmp;

			// compress to new phiTarget by new deltaPhi
			cout << "	** QS compression protocol i = " << i << " from phi = " << phiTmp << " to phiTarget = " << phiTargetTmp << " by dphi = " << deltaPhiTmp << endl;
			packingObject.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftol);

			// increment dphi by a factor of 10
			deltaPhiTmp *= 10;

			// print jammed configuration and vdos
			phiTmp = packingObject.packingFraction();
			cout << "	** Printing compressed state at i = " << i << ", phi = " << phiTmp << " to jam file" << endl;
			packingObject.printJammedConfig();
			packingObject.vdos();
		}

		// if still not confluent, compress to confluency
		if (phiTmp < phiTarget){
			// distance from current packing fraction to target
			phiRange = phiTarget - phiTmp;

			// determine number of compressive steps to target confluency
			NSTEPS = ceil(phiRange/deltaPhiTmp);

			// recompute dphi
			deltaPhiTmp = phiRange/NSTEPS;

			// choose to print during compression or not
			if (NSTEPS < 2*PRINTSTEPS){
				// compress to target packing fraction, only store configurations at the end
				cout << "	** QS compression protocol to final phiTarget = " << phiTarget << endl;
				packingObject.qsIsoCompression(phiTarget, deltaPhiTmp, Ftol);
			}
			else{
				// compress in steps, output every NSKIP steps
				NSKIP = NSTEPS/PRINTSTEPS;

				// loop over compression steps until target reached
				while (phiTmp < phiTarget){
					// temporary target is NSKIPS ahead
					phiTargetTmp = phiTmp + NSKIP*deltaPhiTmp;

					// compress
					cout << "	** QS compression protocol to phiTarget = " << phiTargetTmp << endl;
					packingObject.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftol);

					// check if over true phiTarget, if so then break
					phiTmp = packingObject.packingFraction();
					if (phiTmp > phiTarget)
						break;
					else{
						cout << "	** Printing compressed state at phi = " << phiTmp << " to jam file" << endl;
						packingObject.printJammedConfig();
						packingObject.vdos();
					}
				}
			}
		}
	}

	// Print final confluent config to jammed file
	cout << "	** Printing final confluent state at phiC = " << packingObject.packingFraction() << " to jam file" << endl;
	packingObject.printJammedConfig();
	packingObject.vdos();

	cout << "	** FINISHED COMPRESSING ABOVE JAMMING, ENDING MAIN FILE" << endl;
	return 0;
}