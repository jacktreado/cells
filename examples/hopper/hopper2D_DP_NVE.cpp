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
const int NV 					= 12;			// number of vertices
const int seed 					= 20;			// initial seed
const int NT 					= 1e5;			// number of time steps for flow simulation
const int NPRINT 				= 2e2;			// number of steps between printing
const double smallRadius 		= 0.5;			// radius fo smaller particles (diameter is length unit)
const double sizeRatio 			= 1.4;			// ratio of small diameter to large diameter
const double w0 				= 10.0;			// width of hopper reservoir (in units of small diameter)
const double w 					= 1.5;			// orifice width (in units of small diameter)
const double th 				= PI/6.0;		// hopper angle (pi - th = deflection angle from horizontal)
const double phi0 				= 0.4;			// initial packing fraction
const double T 					= 1e-2;			// constant temperature
const double timeStepMag 		= 0.05;			// time step

// force parameters
const double ka 			= 1.0;				// area force constant (should be = 1)
const double kl 			= 1.0;				// perimeter force constant
const double kb 			= 0.0;				// bending force constant
const double gam 			= 0.0;				// surface tension force constant
const double kint 			= 0.5;				// interaction energy constant
const double a 				= 0.0;				// attraction parameter 
const double del 			= 1.0;				// radius of vertices in units of l0

// main function
int main()
{
	// local variables
	int ci;

	// particle calA0
	double calA0 = 1.0;

	// damping
	double b = 0.1;

	// seed random number generator
	srand48(4809*seed);

	// output files
	string posFile = "hopperDP_pos.test";
	string enFile = "hopperDP_en.test";
	// string posFile = "hopperSP_pos.test";

	// determine L from w0, w, th
	double L = 0.5*(w0 - w)/tan(th);

	// initial number of cells
	int NCELLS = 20;

	// determine scale of reservoir size to make sure that phi \approx 2 given width
	double Lmin = (NCELLS*PI)/(L*w0*3.0);

	// initialize radii
	vector<double> radii(NCELLS,0.0);
	for (ci=0; ci<NCELLS; ci++){
		if (ci < round(0.5*NCELLS))
			radii.at(ci) = smallRadius;
		else
			radii.at(ci) = smallRadius*sizeRatio;
	}

	// instantiate object
	cout << "	** Instantiating object with NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,L,seed);

	// set deformability, force values
	packingObject.forceVals(calA0,ka,kl,gam,kb,kint,del,a);

	// update time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// open print objects
	cout << "	** Opening printing objects for positions and energy " << endl;
	packingObject.openPackingObject(posFile);
	packingObject.openEnergyObject(enFile);

	// initialize positions in hopper reservoir
	cout << "	** Initially placing particles in hopper using SP model" << endl;
	packingObject.initializeHopperSP(radii,w0,w,th,Lmin,NV);

	// check hopper NVE
	double g = 0.1;
	cout << "	** Running hopper NVE with g = " << g << endl;
	// packingObject.hopperDPNVE(w0,w,th,g,T);
	packingObject.flowHopperDP(w0,w,th,g,b);

	return 0;
}