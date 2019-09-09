/*

	Test .cpp file to test forces and print functions
	on an individual, motile cell on a disordered background

*/

// include file
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
	ofstream vertexPrintObject;
	ofstream energyPrintObject;
	string vertexPositionsStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/singleCells/motileCellPositions.dat";
	string energyStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/singleCells/motileCellEnergy.dat";

	// open output files
	vertexPrintObject.open(vertexPositionsStr.c_str());
	if (!vertexPrintObject.is_open()){
		std::cout << "	ERROR: vertex positions output file not open, file string " << vertexPositionsStr << " is not functioning, ending." << std::endl; 
		return 1;
	}

	energyPrintObject.open(energyStr.c_str());
	if (!energyPrintObject.is_open()){
		std::cout << "	ERROR: energy output file not open, file string " << energyStr << " is not functioning, ending." << std::endl; 
		return 1;
	}


	// box variables
	double L = 40.0;

	// Basic cell variables
	int NV = 20;
	double kl,ka,gam,kb,kint,l0,a0,del,segmentMass;

	// set l0 to 1
	l0 = 1.0;

	// set interaction distance to a fraction of l0
	del = 0.2*l0;

	// set a0 accordingly (will be units of l0^2 since l0 = 1)
	double asphericity,minAsphericity;
	asphericity = 2.0;
	a0 = (NV*NV*l0*l0)/(4*PI*asphericity);

	// set unit things to 1
	kl = 1.0;
	ka = 1.0;
	kb = 2.0;
	kint = 10.0;

	// set uneeded things to 0
	gam = 0.0;
	// kb = 0.0;

	// get segment mass
	segmentMass = l0*del + (PI*del*0.5*del*0.5);

	// instantiate object
	cout << "instantiating object with a0 = " << a0 << endl;
	deformableParticles2D singleCell(NV);

	// set L
	singleCell.setL(L);

	// set initial position of cell to be the center of the box
	singleCell.setCPos(0,0.5*L);
	singleCell.setCPos(1,0.5*L);

	// set values for dimensional quantities to set scales
	singleCell.setl0(l0);
	singleCell.seta0(a0);
	singleCell.setdel(del);
	singleCell.seta(del);

	// set force constants
	singleCell.setkl(kl);
	singleCell.setka(ka);
	singleCell.setgam(gam);
	singleCell.setkb(kb);
	singleCell.setkint(kint);

	// set object to a regular polygons initially
	cout << "making a regular polygon" << endl;
	singleCell.regularPolygon(a0);




	return 0;
}