/*

	Example file for active
	particles in zebrafish
	boundary condition

*/

// include files
#include "cellPacking2D.h"
#include "deformableParticles2D.h"
#include <sstream>

// use std name space
using namespace std;

// define PI
const double PI = 4.0*atan(1);

// length paramaters
const int NT 					= 1e3;
const int NPRINT				= 100;

// simulation constants
const double sizeDispersion 	= 0.1;		// size dispersion (std dev of cell sizes)
const double timeStepMag		= 0.01;			// time step scale

// disk constants
const double phiInit 			= 0.75;			// initial packing fraction of disks (sets boundary)
const double phiTarget	 		= 0.99;			// target packing fraction of dp particles

// force constants
const double kl 				= 1.0;				// perimeter force constant
const double ka 				= 1.0;				// area force constant
const double gam 				= 0.0;				// surface tension force constant
const double kb 				= 0.1;				// bending energy constant
const double kint 				= 1.0;				// interaction energy constant
const double del 				= 1.0;				// width of vertices in units of l0, vertex sep on regular polygon
const double a 					= 0.0;				// attraction increment

// main function
int main()
{
	// local variables

	// output files
	string posFile = "pos.test";
	string enFile = "en.test";

	// system details
	int NCELLS 		= 64; 
	int seed 		= 3;
	int NV 			= 20;
	double Ltmp 	= 10.0;
	double Lscale 	= 2.0;
	double v0 		= 0.0;
	double vtau 	= 0.1;
	double dh 		= 0.0005;
	double Pthresh	= 0.01;
	double Dr 		= 0.5;
	double calA0 	= 1.1;
	double Ptol 	= 1e-6;
	double dphi 	= 5e-3;

	// vector of radii
	vector<double> radii(NCELLS,0.0);

	// instantiate object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed); 	// NOTE: NEED TO MAKE NEW CONSTRUCTOR, EVERYTHING ELSE DONE IN initializeGel AND regularPolygon FUNCTIONS

	// open print objects
	packingObject.openPackingObject(posFile);
	packingObject.openEnergyObject(enFile);

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initializing DP particles in zebrafish boundary with Lscale = " << Lscale << endl;
	packingObject.initializeActiveStickySP(radii,NV,phiInit,sizeDispersion,Lscale);

	// // set attraction of each particle
	// for (int ci=0; ci<NCELLS; ci++)
	// 	packingObject.cell(ci).seta(a);

	packingObject.setdt(timeStepMag);

	// set new force scales
	packingObject.gelForceVals(calA0,kl,ka,gam,kb,kint,del,a);

	// run sticky SP NVE
	// cout << "	** Running DP NVE in tailbud geometry" << endl;
	packingObject.dpPntNVE(v0);

	// set new force scales
	packingObject.gelForceVals(calA0,kl,ka,gam,1e-4*kb,kint,del,a);

	// dp zebrafish compression
	packingObject.dpPntIsoCompression(Ptol, phiTarget, dphi);

	return 0;
}
























