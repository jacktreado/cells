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
const double T0 				= 1e-8;			// initial velocities for read-in cells

// force parameters
const double ka 				= 1.0;			// area force constant (should be = 1)
const double gam 				= 0.0;			// surface tension force constant
const double kint 				= 0.5;			// interaction energy constant
const double a 					= 0.0;			// attraction parameter 
const double del 				= 1.0;			// radius of vertices in units of l0

// tolerances
const double Ftol 				= 1e-15;		// force tolerance (for FIRE min)
const double pscale 			= 0.8;			// scale to get new pressure

// ouputs
const double dlogpTol 			= 0.1;			// log difference between pressures of different frames

// int main
int main()
{
	// input variables
	int seed;
	double calA0, kl, kb, pTarget, Ptol, dphi;

	// inputs from command line
	string inputFile 			= "examples/jamming/dpmb.pos";
	string energyFile 			= "zeropEn.test";
	string jammingFile 			= "zeropPos.test";
	string vdosFile 			= "zeropVDOS.test";

	// simulation values
	calA0 = 1.02;
	kl = 1.0;
	kb = 0.01;
	pTarget = 1e-13;
	seed = 1;

	// initial pressure tolerance
	Ptol = 9e-9;

	// initial packing fraction increment
	dphi = 5e-5;

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
	cout << "	** relaxing system Enthalpy H = U + PV to Ptol = " << Ptol << " and Ftol = " << Ftol << endl;
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

	// relax to smaller pressures
	int i, NSTEPS;
	NSTEPS = round((log10(pTarget)-log10(Ptol))/log10(pscale));
	cout << "	** Beginning decompression protocol, scaling target pressure by " << pscale << " each time, target NSTEPS = " << NSTEPS << endl;

	// loop over steps, scale Ptol each time
	for(i=0; i<NSTEPS; i++){
		// scale tolerance
		Ptol *= pscale;

		// minimize enthalpy to target pressure
		cout << "	** Minimize enthalpy to Ptol = " << Ptol << endl;
		packingObject.enthalpyMin(dphi, Ftol, Ptol);
		packingObject.vdos();
	}

	cout << "	** FINISHED RELAXATION TO TARGET PRESSURE, ENDING MAIN FILE" << endl;
	return 0;
}