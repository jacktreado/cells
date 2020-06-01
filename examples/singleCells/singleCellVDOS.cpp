/*

	Example .cpp file to test single cell VDOS computation

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
const double timeStepMag 		= 0.02;			// time step in MD unit
const double phiDisk 			= 0.2;			// initial phi of SP disks

// target packing fraction (confluence)
const double phiTarget 			= 1.0;

// force parameters
const double kl 				= 1.0;			// perimeter force constant
const double kb 				= 0.01;			// bending force constant
const double ka 				= 1.0;			// area force constant (should be = 1)
const double gam 				= 0.0;			// surface tension force constant
const double kint 				= 0.05;			// interaction energy constant
const double a 					= 0.0;			// attraction parameter 
const double del 				= 1.0;			// radius of vertices in units of l0

// force tolerance
const double Ftol 				= 1e-14;		// force tolerance (for FIRE min)

int main(){
	// initial single cell with NV vertices
	int NCELLS = 1;
	double Ltmp = 1.0;
	double calA0 = 1.06;
	int seed = 1;
	int NV = 32;
	double sizeRatio = 1.0;
	double sizeFraction = 1.0;
	double T0 = 1e-6;

	// instantiate object
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed);

	// initialize single cell
	packingObject.initializeBidisperse(NV, phiDisk, sizeRatio, sizeFraction, del);

	// set force values
	packingObject.forceVals(calA0,ka,kl,gam,kb,kint,del,a);

	// set time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// use FIRE to relax particle shape to desired force tolerance
	cout << "	** Relaxing single cell" << endl;
	double Ktest, Ftest;
	packingObject.initializeVelocities(T0);
	packingObject.fireMinimizeF(Ftol,Ftest,Ktest);
	cout << "	** Particle relaxed, Ftest = " << Ftest << ", Ktest = " << Ktest << endl;

	// open file objects
	string vdosFile = "singleCellVDOS.test";
	string posFile = "singleCellPos.test";
	packingObject.openPackingObject(posFile);
	packingObject.openStatObject(vdosFile);

	// place cell in center of box
	double L = packingObject.getL(0);
	double vx, vy;
	for (int i=0; i<NV; i++){
		// get position
		vx = packingObject.cell(0).vrel(i,0);
		vy = packingObject.cell(0).vrel(i,1);

		// move to box center
		packingObject.cell(0).setVPos(i,0,0.5*L + vx);
		packingObject.cell(0).setVPos(i,1,0.5*L + vy);
	}
	packingObject.cell(0).updateCPos();
	

	// print positions of relaxed cell shape
	cout << "	** Printing single cell positions" << endl;
	packingObject.printSystemPositions();

	// compute vdos
	cout << "	** Computing single cell VDOS" << endl;
	packingObject.vdos();


	// return
	return 0;
}