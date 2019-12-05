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
const int NT 					= 5e4;
const int NPRINT				= 200;

// simulation constants
const double sizeDispersion 	= 0.1;		// size dispersion (std dev of cell sizes)
const double timeStepMag		= 0.05;			// time step scale

// disk constants
const double phiDisk	 		= 0.9;			// initial packing fraction of disks (sets boundary)

// force constants
const double a 					= 0.1;			// attractive parameter

// main function
int main()
{
	// local variables

	// output files
	string posFile = "pos.test";
	string enFile = "en.test";

	// system details
	int NCELLS 		= 256; 
	int seed 		= 3;
	int NV 			= 24;
	double Ltmp 	= 10.0;
	double R0 		= 10.0;
	double v0 		= 0.05;
	double vtau 	= 0.1;
	double dh 		= 0.0005;
	double Pthresh	= 0.01;
	double Dr 		= 0.5;

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
	cout << "	** Initializing particles at phiDisk = " << phiDisk << " using SP model in zebrafish boundary with R0 = " << R0 << endl;
	packingObject.initializeActiveZebrafish(radii,NV,phiDisk,sizeDispersion,R0);

	// // set attraction of each particle
	// for (int ci=0; ci<NCELLS; ci++)
	// 	packingObject.cell(ci).seta(a);

	// run sticky SP NVE
	// cout << "	** Running active pipe flow using sticky SP model" << endl;
	// packingObject.spActivePipeFlow(radii, a, v0, Dr);
	packingObject.spActiveZebrafishVicsek(radii, a, v0, Dr, vtau, Pthresh, dh);

	return 0;
}
























