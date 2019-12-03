/*

	MAIN FILE for active brownian particles (ABPs)
	in 

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
const double timeStepMag		= 0.05;			// time step scale

// disk constants
const double phiDisk	 		= 0.9;			// initial packing fraction of disks (sets boundary)

// boundary constants
const double R0 				= 10.0;			// radius of horseshoe
const double dh					= 1e-3;			// step size for horseshoe increase
const double LINIT 				= 1.0;			// throwaway variable for initial box size

// cell constants
const int NV 					= 12;			// number of vertices per cell (does not matter for this sim)
	
int main(int argc, char const *argv[])
{
	// variables for main
	int NCELLS, NT, NPRINT, seed;
	double NTtmp, NPRINTtmp, sizeDispersion, a, v0, Dr, Pthresh; 

	// inputs from command line
	string NCELLS_str 			= argv[1];
	string NT_str 				= argv[2];
	string NPRINT_str 			= argv[3];
	string sizeDisp_str 		= argv[4];
	string v0_str 				= argv[5];
	string Dr_str 				= argv[6];
	string Pthresh_str 			= argv[7];
	string a_str 				= argv[8];
	string seed_str				= argv[9];
	string positionFile			= argv[10];
	string energyFile 			= argv[11];

	// load strings into sstream
	stringstream NCELLSss(NCELLS_str);
	stringstream NTss(NT_str);
	stringstream NPRINTss(NPRINT_str);
	stringstream sizeDispss(sizeDisp_str);
	stringstream v0ss(v0_str);
	stringstream Drss(Dr_str);
	stringstream Pthreshss(Pthresh_str);
	stringstream ass(a_str);
	stringstream seedss(seed_str);

	// parse values from strings
	NCELLSss 		>> NCELLS;
	NTss 			>> NTtmp;
	NPRINTss 		>> NPRINTtmp;
	sizeDispss 		>> sizeDispersion;
	v0ss 			>> v0;
	Drss 			>> Dr;
	Pthreshss 		>> Pthresh;
	ass	 			>> a;
	seedss 			>> seed;

	// save values to ints
	NT 		= (int)NTtmp;
	NPRINT 	= (int)NPRINTtmp;

	// vector of radii
	vector<double> radii(NCELLS,0.0);

	// instantiate object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,LINIT,seed); 	// NOTE: NEED TO MAKE NEW CONSTRUCTOR, EVERYTHING ELSE DONE IN initializeGel AND regularPolygon FUNCTIONS

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initializing particles at phiDisk = " << phiDisk << " using SP model in zebrafish boundary with R0 = " << R0 << endl;
	packingObject.initializeActiveZebrafish(radii,NV,phiDisk,sizeDispersion,R0);

	// open print objects
	packingObject.openPackingObject(positionFile);
	packingObject.openEnergyObject(energyFile);

	// run active brownian particle simulation
	cout << "	** Running particles at v0 = " << v0 << ", Dr = " << Dr << ", and a = " << a << endl;
	packingObject.spAciveZebrafishABPs(radii, a, v0, Dr, Pthresh, dh);

	// end program
	cout << endl << "	** Finishing program in main." << endl;
	return 0;
}