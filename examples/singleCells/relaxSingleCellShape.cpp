/*

	Example file to generate a gel of cells
	from an initial bidisperse sphere packing,
	and to decompress at a fixed rate
	rather than QS

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
const int NT 					= 1e1;
const int NPRINT				= 50;

// simulation constants
const double timeStepMag		= 0.01;			// time step scale

// disk constants
const double phiDisk	 		= 0.1;			// initial packing fraction of disks (sets boundary)

// main function
int main()
{
	// output files
	string posFile = "pos.test";
	string enFile = "en.test";

	// system details
	int NCELLS 		= 1; 
	int seed 		= 1;
	int NV 			= 50;
	double Ltmp 	= 10.0;

	// characteristic time scale
	double t0 		= 0.0;

	// force constants
	double calA0 	= 1.5;
	double kl		= 0.1;
	double ka 		= 1.0;
	double gam 		= 0.0;
	double kb 		= 1.0;
	double kint 	= 1.0;
	double del 		= 1.0;
	double a 		= 0.0;

	// preferred angle
	double theta0 	= PI/15.0;

	// relaxation tolerances
	double Ftol = 1e-10;
	double Ktol = 1e-20;

	// instantiate object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed); 	// NOTE: NEED TO MAKE NEW CONSTRUCTOR, EVERYTHING ELSE DONE IN initializeGel AND regularPolygon FUNCTIONS

	// open print files
	packingObject.openPackingObject(posFile);
	packingObject.openEnergyObject(enFile);

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
	packingObject.cell(0).vertexPerturbation(1e-4);

	// print initial position to frame
	packingObject.printSystemPositions();

	// set preferred theta0
	packingObject.cell(0).setc0Angle(theta0);

	// set time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// use FIRE to relax particle shape to desired force tolerance
	cout << "	** Relaxing single particle at phiDisk = " << phiDisk << " using SP model" << endl;
	packingObject.fireMinimizeF(Ftol, Ktol);

	// end function
	return 0;
}








