/*

	Test .cpp file to test FIRE-based 
	shape relaxation

*/

// include files
#include "cellPacking2D.h"
#include "deformableParticles2D.h"

// use std name space
using namespace std;

// define PI
const double PI = 4.0*atan(1);

int main(){
	// local main variables
	int i,t,d;
	double postmp,veltmp,anew;
	double dx,dy,dnorm,dscale;
	string positionsStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/singleCells/singleCellFirePositions.dat";
	string energyStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/singleCells/singleCellFireEnergy.dat";

	// input params
	int NCELLS 	= 1;
	int NT 		= 1e6;
	int NPRINT 	= 1e1;
	double L 	= 10.0*NCELLS;


	// instantiate object
	cout << "	** Instantiating object" << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,L);

	// Basic cell variables (cells will be identical)
	int NV = 30;
	double kl,ka,gam,kb,kint,del,a;
	double asphericity,T0;

	// force and rest params
	kl 			= 1.0;
	ka 			= 1.0;
	gam 		= 0.0;
	kb 			= 10.0;
	kint 		= 1.0;
	del 		= 0.2;
	a 			= 0.0;

	// asphericity
	asphericity = 1.11;

	// temperature
	T0 = 0.0;

	// set values in setting up packing
	cout << "	** Initializing particle quantities" << endl;
	packingObject.initializeMonodisperse(NV,asphericity);
	packingObject.initializeForceConstants(kl,ka,gam,kb,kint);
	packingObject.initializeInteractionParams(del,a);

	// set cell to be in center
	packingObject.cell(0).setCPos(0,0.5*L);
	packingObject.cell(0).setCPos(1,0.5*L);
	packingObject.cell(0).regularPolygon();

	// open print objects
	packingObject.openPackingObject(positionsStr);
	packingObject.openEnergyObject(energyStr);

	/****************

		Simulation

	*****************/

	// simulation params
	double timeStepMag 			= 0.02;
	double initialPhi 			= 0.2;
	double kineticTol 			= 1e-16;
	double forceTol 			= 1e-16;

	// set time step
	packingObject.setdt(timeStepMag);

	// set initial packing fraction
	packingObject.setPackingFraction(initialPhi);

	// relax shape
	packingObject.shapeRelax(1);

	return 0;
}