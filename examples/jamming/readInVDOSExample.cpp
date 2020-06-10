/*

	Example .cpp file to compress NCELLS to jamming
	and to compute shear modulus

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
const double deltaPhi0 			= 1e-4;			// initial delta phi

// force parameters
const double ka 			= 1.0;				// area force constant (should be = 1)
const double gam 			= 0.0;				// surface tension force constant
const double kint 			= 0.1;				// interaction energy constant
const double a 				= 0.0;				// attraction parameter 
const double del 			= 1.0;				// radius of vertices in units of l0

// tolerances
const double Ftol 			= 1e-12;			// force tolerance (for FIRE min)
const double Ptol 			= 1e-8;			// pressure tolerance

// main function
int main()
{
	// local variables

	// input file
	string inputFile = "/Users/JackTreado/Jamming/CellSim/viz/jamming/data/bidcells_N16_NV16_calA1.02_kl1.0_kb0_seed1.jam";

	// output files
	string vdosFile = "vdos_comp.test";
	string jamFile = "jam_comp.test";
	string enFile = "en_comp.test";

	// system details
	double seed 	= 1;
	double T0 		= 1e-6;

	// mechanical parameters
	double kl = 1.0;
	double kb = 0;
	double calA0 = 1.02;
	double Fcheck, Kcheck;

	// instantiate main packing object
	cout << "	** Reading in from file " << inputFile << endl;
	cellPacking2D packingObject(inputFile,T0,seed);

	// set deformability, force values
	packingObject.forceVals(calA0,ka,kl,gam,kb,kint,del,a);

	// update time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// open files
	packingObject.openJamObject(jamFile);
	packingObject.openStatObject(vdosFile);
	packingObject.openEnergyObject(enFile);

	// set NT and NPRINT
	packingObject.setNT(NT);
	packingObject.setNPRINT(NPRINT);

	// compress to set packing fraction using FIRE, pressure relaxation
	cout << "	** relaxing system with Ftol = " << Ftol << endl;
	packingObject.fireMinimizeF(Ftol, Fcheck, Kcheck);
	cout << "Fcheck = " << Fcheck << endl;
	cout << "Kcheck = " << Kcheck << endl;

	// compute pressure
	double Pcheck = 0.5*(packingObject.getSigmaXX() + packingObject.getSigmaYY())/(packingObject.getNCELLS()*packingObject.getL(0)*packingObject.getL(0));
	cout << "Pcheck = " << Pcheck << endl;

	// if Pcheck > 2*Ptol, compute VDOS, else, find nearest jammed state
	double dphi = deltaPhi0;
	packingObject.printJammedConfig();

	cout << "	** computing VDOS, printing to " << vdosFile << endl;
	packingObject.vdos();

	// decrease packing fraction by small steps until unjamming
	double phi = packingObject.packingFraction();
	double rscale = sqrt((phi + deltaPhi0)/phi);
	int itmax = 1e3;
	int it = 0;
	while(Pcheck < 1e-5 && it < itmax){
		// iterator
		it++;

		// decrease by packing fraction
		cout << "	** Compression protocol it = " << it << " from phi = " << phi << " to phi + dphi = " << phi + deltaPhi0 << " by dphi = " << deltaPhi0 << endl;
		packingObject.scaleLengths(rscale);

		// minimize
		packingObject.fireMinimizeF(Ftol, Fcheck, Kcheck);

		// updated packing fraction
		phi = packingObject.packingFraction();

		// updated pressure/cell
		Pcheck = 0.5*(packingObject.getSigmaXX() + packingObject.getSigmaYY())/(packingObject.getNCELLS()*packingObject.getL(0)*packingObject.getL(0));

		// print vdos and config
		packingObject.vdos();
		packingObject.printJammedConfig();
	}

	return 0;
}
























