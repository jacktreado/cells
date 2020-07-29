/*

	Example file for 2D soft particles
		in a simplified hopper flow

	Packs particles randomly in a reservoir, 
		lets system settle, then flow begins

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
const int NV 					= 12;			// number of vertices (not relevant for SP model)
const int seed 					= 20;			// initial seed
const int NT 					= 1e5;			// number of time steps for flow simulation
const int NPRINT 				= 2e3;			// number of steps between printing
const double sizeDispersion 	= 0.1;			// std deviation of polydispersity
const double meanRadius 		= 0.5;			// mean radius (diameter is length unit)
const double w0 				= 10.0;			// width of hopper reservoir (in units of mean diameter sigma)
const double w 					= 2.01;			// orifice width (in units of mean diameter sigma)
const double th 				= PI/3.0;		// hopper angle (pi - th = deflection angle from horizontal)
const double phi0 				= 0.4;			// initial packing fraction

// main function
int main()
{
	// local variables
	int ci;
	double r1, r2, g1;

	// seed random number generator
	srand48(4809*seed);

	// output files
	string posFile = "pos.test";

	// determine L from w0, w, th
	double L = 0.5*(w0 - w)/tan(th);

	// determine number of particles given initial geometry, packing fraction
	double Ndbl = 2.0*phi0*(3.0*w0)*(w0-w)/(PI*tan(th));
	int NCELLS = round(Ndbl);

	// determine scale of reservoir size to make sure that phi \approx 2 given width
	double Lmin = (NCELLS*PI)/(L*w0*3.0);

	// initialize radii as gaussian random variables
	vector<double> radii(NCELLS,0.0);
	for (ci=0; ci<NCELLS; ci++){
		// generate random numbers
		r1 = drand48();
		r2 = drand48();

		// calculate gaussian random variable using Box-Muller transform
		g1 = sqrt(-2.0*log(r1))*cos(2*PI*r2);

		// get radius
		radii.at(ci) = g1*sizeDispersion + meanRadius;
	}

	// instantiate object
	cout << "	** Instantiating object with NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,L,seed);

	// open print objects
	cout << "	** Opening printing object for positions " << endl;
	packingObject.openPackingObject(posFile);

	// initialize positions in hopper reservoir
	cout << "	** Relaxing particle positions using SP model" << endl;
	packingObject.initializeHopperSP(radii,w0,w,th,Lmin);


	// flow particles through hopper with force strength g
	double g = 1e-2;
	cout << "	** Running hopper FLOW with g = " << g << endl;
	packingObject.flowHopperSP(radii,w0,w,th,g);

	return 0;
}