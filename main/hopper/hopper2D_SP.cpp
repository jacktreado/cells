/*

	Main file to start with a random distribution of 
	particles in a hopper reservoir, and then flow
	toward orifice with body force of strength g/f_0 (f_0 = 1)

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
const int NPRINT 				= 2e2;			// number of steps between printing
const int NV 					= 12;			// number of vertices (not important for SP model)
const double meanRadius 		= 0.5;			// mean radius (diameter is length unit)
const double th 				= PI/3.0;		// hopper angle (pi - th = deflection angle from horizontal)

// int main
int main(int argc, char const *argv[])
{
	// local variables
	int ci;
	double r1, r2, grv, Lmin;

	// inputs from command line
	string NCELLS_str 			= argv[1];
	string NT_str 				= argv[2];
	string sizeDisp_str 		= argv[3];
	string g_str 				= argv[4];
	string w0_str				= argv[5];
	string w_str 				= argv[6];
	string seed_str 			= argv[7];
	string positionFile			= argv[8];

	// load strings into sstream
	stringstream NCELLSss(NCELLS_str);
	stringstream NTss(NT_str);
	stringstream sizeDispss(sizeDisp_str);
	stringstream gss(g_str);
	stringstream w0ss(w0_str);
	stringstream wss(w_str);
	stringstream seedss(seed_str);

	// read-in variables
	int NCELLS, NT, seed;
	double NT_dbl, L, sizeDispersion, g, w0, w;

	// parse values from strings
	NCELLSss 		>> NCELLS;
	NTss 			>> NT_dbl;
	sizeDispss 		>> sizeDispersion;
	gss				>> g;
	w0ss 			>> w0;
	wss 			>> w;
	seedss 			>> seed;

	// get integer NT from input double value
	NT = (int)round(NT_dbl);

	// seed random number generator
	srand48(341234018*seed);

	// determine L from hopper geometry parameters
	L = 0.5*(w0 - w)/tan(th);

	// initialize radii as gaussian random variables
	vector<double> radii(NCELLS,0.0);
	for (ci=0; ci<NCELLS; ci++){
		// generate random numbers
		r1 = drand48();
		r2 = drand48();

		// calculate gaussian random variable using Box-Muller transform
		grv = sqrt(-2.0*log(r1))*cos(2*PI*r2);

		// get radius
		radii.at(ci) = grv*sizeDispersion + meanRadius;
	}

	// determine scale of reservoir size to make sure that phi \approx 2 given width
	Lmin = (NCELLS*PI)/(L*w0*3.0);

	// instantiate object
	cout << "	** Instantiating object with NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,L,seed);

	// open print objects
	cout << "	** Opening printing objects for positions and energy " << endl;
	packingObject.openPackingObject(positionFile);

	// initialize positions in hopper reservoir
	cout << "	** Relaxing particle positions using SP model" << endl;
	packingObject.initializeHopperSP(radii,w0,w,th,Lmin);

	// flow particles through hopper with force strength g
	cout << "	** Running hopper FLOW with g = " << g << endl;
	packingObject.flowHopperSP(radii,w0,w,th,g);

	return 0;
}