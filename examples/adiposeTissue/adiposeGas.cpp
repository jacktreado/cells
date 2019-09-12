/*

	Gas of adipose tissue cells

	simplest possible, a slate of N-1 adipose cells,
	and a single deformable tumor cell that is bouncing around at fixed packing fraction
	with NVE dynamics

	TUMOR

*/


// include file
#include "cellPacking2D.h"
#include "deformableParticles2D.h"

// use std name space
using namespace std;

// define PI
const double PI = 4.0*atan(1);

int main(){
	// local main variables
	int i;
	int seed = 1;

	// output strings
	string positionsStr = "examples/adiposeTissue/adiposeGasPositions.test";
	string energyStr = "examples/adiposeTissue/adiposeGasEnergy.test";

	// input params
	int NCELLS = 17;
	int NTUMORCELLS = 1;
	int tumorNV = 12;
	int adiposeNV = 35;
	double tumorCalA = 1.17;
	double adiposeCalA = 1.05;
	
	// instantiate object
	cellPacking2D tissueObject(NCELLS,NTUMORCELLS,tumorNV,adiposeNV,tumorCalA,adiposeCalA,seed);

	// open output files
	tissueObject.openPackingObject(positionsStr);
	tissueObject.openEnergyObject(energyStr);

	// simulation setup
	int NT 				= 1e5;
	int NPRINT 			= 1e2;
	double initialPhi 	= 0.65;
	double dtMag 		= 0.03;

	tissueObject.setNT(NT);
	tissueObject.setNPRINT(NPRINT);
	tissueObject.setdt(dtMag);
	tissueObject.setPackingFraction(initialPhi);

	// perform initial overlap relief
	tissueObject.overlapRelief();

	// temperature differences
	double adiposeT, tumorT, tempRatio;
	double adiposeC, adiposel;

	tempRatio = 1.0;
	adiposeT = 1.0;
	tumorT = adiposeT*tempRatio;

	// loop over cells, give different initial velocities
	for (i=0; i<NCELLS; i++){
		if (i < NTUMORCELLS)
			tissueObject.initializeVelocities(i,tumorT);
		else
			tissueObject.initializeVelocities(i,adiposeT);
	}

	// run NVE simulation with constant force on tumor cell
	double forceScale = 0.1;
	double adiposeDamping = 0.75;
	tissueObject.tumorForce(NTUMORCELLS, forceScale, adiposeDamping);
	// tissueObject.tumorNVE();

	// return 0 to end
	return 0;
}
