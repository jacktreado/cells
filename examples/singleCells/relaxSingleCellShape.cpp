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
const int NPRINT				= 200;

// simulation constants
const double timeStepMag		= 0.03;			// time step scale

// disk constants
const double phiDisk	 		= 0.1;			// initial packing fraction of disks (sets boundary)

// main function
int main()
{
	// local variables
	double Fcheck, Kcheck;

	// output files
	string posFile = "pos.test";

	// system details
	int NCELLS 		= 1; 
	int seed 		= 1;
	int NV 			= 23;
	double Ltmp 	= 1.0;

	// force constants
	double calA0 	= 1.02;
	double kl		= 1.0;
	double ka 		= 1.0;
	double gam 		= 0.0;
	double kb 		= 0.01;
	double kint 	= 1.0;
	double del 		= 1.0;
	double a 		= 0.0;

	// relaxation tolerances
	double Ftol = 1e-15;

	// instantiate object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed); 	// NOTE: NEED TO MAKE NEW CONSTRUCTOR, EVERYTHING ELSE DONE IN initializeGel AND regularPolygon FUNCTIONS

	// open print files
	packingObject.openPackingObject(posFile);

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

	// initialize asphericity (hard code from matlab output)
	double l0 = 0.248523189510815;
	double a0 = 2.5331841512012;

	// set asphericity values
	packingObject.cell(0).seta0(a0);
	packingObject.cell(0).setl0(l0);

	// initialize relative positions (hard code from matlab output)
	packingObject.cell(0).setVPos(0,0,6.72674288957704);
	packingObject.cell(0).setVPos(1,0,6.78932241064127);
	packingObject.cell(0).setVPos(2,0,6.77019355336481);
	packingObject.cell(0).setVPos(3,0,6.6681063224058);
	packingObject.cell(0).setVPos(4,0,6.49672565716751);
	packingObject.cell(0).setVPos(5,0,6.27879279014721);
	packingObject.cell(0).setVPos(6,0,6.03709474196122);
	packingObject.cell(0).setVPos(7,0,5.78915847869948);
	packingObject.cell(0).setVPos(8,0,5.54722426797899);
	packingObject.cell(0).setVPos(9,0,5.32116954178392);
	packingObject.cell(0).setVPos(10,0,5.12195601814524);
	packingObject.cell(0).setVPos(11,0,4.96403977604661);
	packingObject.cell(0).setVPos(12,0,4.86522797278833);
	packingObject.cell(0).setVPos(13,0,4.84241837868629);
	packingObject.cell(0).setVPos(14,0,4.90384694286031);
	packingObject.cell(0).setVPos(15,0,5.04307205176409);
	packingObject.cell(0).setVPos(16,0,5.24076557693787);
	packingObject.cell(0).setVPos(17,0,5.47315619587024);
	packingObject.cell(0).setVPos(18,0,5.719785490424);
	packingObject.cell(0).setVPos(19,0,5.9660192478911);
	packingObject.cell(0).setVPos(20,0,6.20124107011099);
	packingObject.cell(0).setVPos(21,0,6.415428749558);
	packingObject.cell(0).setVPos(22,0,6.59605565883316);


	packingObject.cell(0).setVPos(0,1,2.12868177305827);
	packingObject.cell(0).setVPos(1,1,2.36842628464961);
	packingObject.cell(0).setVPos(2,1,2.61541112729426);
	packingObject.cell(0).setVPos(3,1,2.84111427228937);
	packingObject.cell(0).setVPos(4,1,3.02003621052704);
	packingObject.cell(0).setVPos(5,1,3.13803884382293);
	packingObject.cell(0).setVPos(6,1,3.19314369894159);
	packingObject.cell(0).setVPos(7,1,3.1903807481705);
	packingObject.cell(0).setVPos(8,1,3.13597363287517);
	packingObject.cell(0).setVPos(9,1,3.03404178141157);
	packingObject.cell(0).setVPos(10,1,2.88642865022798);
	packingObject.cell(0).setVPos(11,1,2.69535256476083);
	packingObject.cell(0).setVPos(12,1,2.46809126207777);
	packingObject.cell(0).setVPos(13,1,2.22139699946827);
	packingObject.cell(0).setVPos(14,1,1.98142023179465);
	packingObject.cell(0).setVPos(15,1,1.7765110176007);
	packingObject.cell(0).setVPos(16,1,1.62711801119083);
	packingObject.cell(0).setVPos(17,1,1.54090739199095);
	packingObject.cell(0).setVPos(18,1,1.51555394120638);
	packingObject.cell(0).setVPos(19,1,1.54482535177002);
	packingObject.cell(0).setVPos(20,1,1.62333634813372);
	packingObject.cell(0).setVPos(21,1,1.7482707374639);
	packingObject.cell(0).setVPos(22,1,1.91808535139912);


	// loop over vertices, scale by a0
	double rho0 = sqrt(a0);
	for (int i=0; i<NV; i++){
		cout << "i = " << i << "; ";
		for (int d=0; d<NDIM; d++){
			packingObject.cell(0).setVPos(i,d,packingObject.cell(0).vpos(i,d)/rho0);
			cout << d << " pos = " << packingObject.cell(0).vpos(i,d) << "; ";
		}
		cout << endl;
	}
	l0 = l0/rho0;
	a0 = 1.0;

	// update cpos
	packingObject.cell(0).updateCPos();

	// loop over vertices, update vpos based on cpos
	for (int i=0; i<NV; i++){
		for (int d=0; d<NDIM; d++)
			packingObject.cell(0).setVPos(i,d,0.5*Ltmp + packingObject.cell(0).vrel(i,d));
	}

	// set asphericity values
	packingObject.cell(0).seta0(a0);
	packingObject.cell(0).setl0(l0);

	// initialize cell velocities
	packingObject.cell(0).setCVel(0,0.0);
	packingObject.cell(0).setCVel(1,0.0);

	// set force values
	packingObject.forceVals(calA0,ka,kl,gam,kb,kint,del,a);

	// set time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// print initial position to frame
	packingObject.printSystemPositions();

	// use FIRE to relax particle shape to desired force tolerance
	cout << "\t ** Running F min using FIRE on single particle to Ftol = " << Ftol << endl;
	packingObject.fireMinimizeF(Ftol, Fcheck, Kcheck);
	cout << "\t ** F min complete, Fcheck = " << Fcheck << ", Kcheck = " << Kcheck << endl;

	// print final position to frame
	packingObject.printSystemPositions();

	// end function
	return 0;
}








