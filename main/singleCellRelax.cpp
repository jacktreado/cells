/*

	Code to find the ground state shape of
	a single cell by changing mechanical and
	shape parameters

	Uses FIRE algorithm to relax forces

*/

// include files
#include "cellPacking2D.h"
#include "deformableParticles2D.h"
#include <sstream>

// use std name space
using namespace std;

// define PI
const double PI = 4.0*atan(1);

// dimension of space
const int NDIM 					= 2;

// length paramaters
const int NT 					= 1e3;
const int NPRINT				= 50;

// simulation constants
const int NCELLS 				= 1;			// only one cell
const double timeStepMag		= 0.005;		// time step scale
const double gam				= 0.0;			// surface tension force scale
const double del 				= 1.0;			// vertex size (set to 1)
const double kint 				= 1.0;			// interaction spring constant (set to 1)
const double a 					= 1.0;			// attraction parameter (set to 1)

// area spring: SET TO 1, SETS ENERGY SCALE
const double ka 				= 1.0;

// disk constants
const double phiDisk	 		= 0.1;			// initial packing fraction of cell-disk (sets boundary)

// minimization tolerances
const double Ftol 				= 1e-10;		// force tolerance
const double Ktol 				= 1e-20;		// kinetic energy tolerance

// main function
int main(int argc, char const *argv[])
{
	// local variables
	int NV, seed;
	double calA0, kl, kb, th0, t0;

	// inputs from command line
	string NV_str 				= argv[1];
	string calA0_str 			= argv[2];
	string kl_str 				= argv[3];
	string kb_str 				= argv[4];
	string th0_str 				= argv[5];
	string seed_str				= argv[6];
	string positionFile			= argv[7];

	// load strings into sstream
	stringstream NVss(NV_str);
	stringstream calA0ss(calA0_str);
	stringstream klss(kl_str);
	stringstream kbss(kb_str);
	stringstream th0ss(th0_str);
	stringstream seedss(seed_str);

	// parse values from strings
	NVss 			>> NV;
	calA0ss 		>> calA0;
	klss 			>> kl;
	kbss 			>> kb;
	th0ss	 		>> th0;
	seedss 			>> seed;

	// scale theta0 by PI
	th0 *= PI;

	// temporary box length; will be modified in initialization
	double Ltmp 	= 1.0;

	// instantiate object
	cout << "	** Instantiating object for single particle to be relaxed" << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed);

	// initialize boundary for cell of packing fraction phiDisk
	Ltmp = sqrt(PI/phiDisk);
	packingObject.setL(0,Ltmp);
	packingObject.setL(1,Ltmp);

	// set up cell object
	for (int d=0; d<NDIM; d++){
		packingObject.cell(0).setL(d,Ltmp);
		packingObject.cell(0).setpbc(d,1);
	}
	
	// set number of vertices
	packingObject.cell(0).setNV(NV);

	// initialize arrays
	packingObject.cell(0).initializeVertices();
	packingObject.cell(0).initializeCell();

	// initialize asphericity
	double l0 = 2.0*sin(PI/NV);
	double a0 = 0.5*NV*sin(2.0*PI/NV);
	packingObject.cell(0).setl0(l0);
	packingObject.cell(0).seta0(a0);

	// initialize cell position
	packingObject.cell(0).setCPos(0,0.5*Ltmp);
	packingObject.cell(0).setCPos(1,0.5*Ltmp);

	// initialize cell velocities
	packingObject.cell(0).setCVel(0,0.0);
	packingObject.cell(0).setCVel(1,0.0);


	// initialize vertices as regular polygon
	packingObject.cell(0).regularPolygon();

	// set del to 1, as this is the bumpy cell model
	packingObject.cell(0).setdel(del);

	// set force values
	packingObject.gelForceVals(calA0,kl,ka,gam,kb,kint,del,a);

	// get box length
	Ltmp = packingObject.getL(0);

	// set cell position to box center, deform position
	packingObject.cell(0).setCPos(0,0.5*Ltmp);
	packingObject.cell(0).setCPos(1,0.5*Ltmp);

	// deform vertex positions slightly
	packingObject.cell(0).vertexPerturbation(1e-2);

	// print initial position to frame
	packingObject.printSystemPositions();

	// set preferred theta0
	packingObject.cell(0).setc0Angle(th0);

	// set time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// use FIRE to relax particle shape to desired force tolerance
	cout << "	** Relaxing single particle at phiDisk = " << phiDisk << " force-based FIRE algorithm" << endl;
	packingObject.fireMinimizeF(Ftol, Ktol);

	// print positions to file
	cout << "	** Printing vetex positions to file" << endl;
	packingObject.printSystemPositions();

	// print that end of main has been found
	cout << "	** Reached end of main, single cell relaxation protocol has ended." << endl;

	// end function
	return 0;
}