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
	int seed = 1;

	// output strings
	string positionsStr = "examples/adiposeTissue/adiposeGasPositions.test";
	string energyStr = "examples/adiposeTissue/adiposeGasEnergy.test";

	// input params
	int NCELLS = 17;
	int NT = 1e7;
	int NPRINT = 1e2;
	int NTUMORCELLS = 1;
	int tumorNV = 15;
	int adiposeNV = 30;
	double tumorCalA = 1.17;
	double adiposeCalA = 1.05;
	
	// instantiate object
	cellPacking2D adiposeObject(NCELLS,NTUMORCELLS,tumorNV,adiposeNV,tumorCalA,adiposeCalA,seed);

	// setup simulation
	double adiposeT, tumorT, tempRatio;
	double adiposeC, adiposel;

	// temperature differences
	adipose


	adiposeObject.setNT(NT);
	adiposeObject.setNPRINT(NPRINT);

	// run NVE simulation
	adiposeObject.tumorNVE();



	// return 0 to end
	return 0;
}
