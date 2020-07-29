/*

	Main file to start flow of bidisperse 
	repulsive DPM particles in hopper

	Input parameters:
		-- NCELLS: 				number of cells
		-- NT:					length of simulation
		-- sizeDispersion: 		std dev of particles sizes
		-- g: 					scale of body forces (relative to spring force f_0 = 1)
		-- w0:					width of the hopper reservoir
		-- w: 					width of the hopper orifice
		-- seed: 				initial seed 


	Files to write to
		-- positionFile: 		configuration/force information during hopper flow simulation

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
const int NPRINT 				= 1e3;			// number of steps between printing
const double smallRadius 		= 0.5;			// radius fo smaller particles (diameter is length unit)
const double sizeRatio 			= 1.4;			// ratio of small diameter to large diameter
const double th 				= PI/4.0;		// hopper angle (pi - th = deflection angle from horizontal)
const double b 					= 0.1;			// damping coefficient
const double timeStepMag 		= 0.05;			// time step

// force parameters
const double ka 			= 1.0;				// area force constant (should be = 1)
const double kint 			= 1.0;				// interaction energy constant
const double a 				= 0.0;				// attraction parameter 
const double del 			= 1.0;				// radius of vertices in units of l0
const double calA0 			= 1.0;				// target shape parameter, set to 1 for hopper simulations

// int main
int main(int argc, char const *argv[])
{
	// local variables
	int ci;
	double Lmin;

	// inputs from command line
	string NCELLS_str 			= argv[1];
	string NV_str 				= argv[2];
	string NT_str 				= argv[3];
	string g_str 				= argv[4];
	string w0_str				= argv[5];
	string w_str 				= argv[6];
	string kl_str 				= argv[7];
	string gam_str 				= argv[8];
	string kb_str 				= argv[9];
	string seed_str 			= argv[10];
	string positionFile			= argv[11];

	// load strings into sstream
	stringstream NCELLSss(NCELLS_str);
	stringstream NVss(NV_str);
	stringstream NTss(NT_str);
	stringstream gss(g_str);
	stringstream w0ss(w0_str);
	stringstream wss(w_str);
	stringstream klss(kl_str);
	stringstream gamss(gam_str);
	stringstream kbss(kb_str);
	stringstream seedss(seed_str);

	// read-in variables
	int NCELLS, NV, NT, seed;
	double NT_dbl, L, g, w0, w, kl, gam, kb;

	// parse values from strings
	NCELLSss 		>> NCELLS;
	NVss 			>> NV;
	NTss 			>> NT_dbl;
	gss				>> g;
	w0ss 			>> w0;
	wss 			>> w;
	klss 			>> kl;
	gamss 			>> gam;
	kbss 			>> kb;
	seedss 			>> seed;

	// get integer NT from input double value
	NT = (int)round(NT_dbl);

	// initialize radii
	vector<double> radii(NCELLS,0.0);
	for (ci=0; ci<NCELLS; ci++){
		if (ci < round(0.5*NCELLS))
			radii.at(ci) = smallRadius;
		else
			radii.at(ci) = smallRadius*sizeRatio;
	}

	// determine scale of reservoir size to make sure that phi \approx 2 given width
	Lmin = 0.2*w0;

	// instantiate object
	cout << "	** Instantiating object with NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Lmin,seed);

	// set deformability, force values
	packingObject.forceVals(calA0,ka,kl,gam,kb,kint,del,a);

	// open print objects
	cout << "	** Opening printing objects for positions and energy " << endl;
	packingObject.openPackingObject(positionFile);

	// initialize positions in hopper reservoir
	cout << "	** Initially placing particles in hopper using SP model" << endl;
	packingObject.initializeHopperDP(radii,w0,w,th,Lmin,NV);

	// update time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// run flow simulation for NT time steps
	cout << "	** Running hopper DPM flow with g = " << g << endl;
	packingObject.flowHopperDP(w0,w,th,g,b);


	// print ending statement
	cout <<" 	** FINISHED RUNNING DP HOPPER FLOW FOR NT = " << NT << " TIME STEPS, ENDING. " << endl;
	return 0;
}