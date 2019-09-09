/*

	Test .cpp file to test forces and print functions
	on an individual cell

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
	string vertexPositionsStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/singleCells/singleCellPositions.dat";
	string energyStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/singleCells/singleCellEnergy.dat";


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
	int NV = 50;
	double kl,ka,gam,kb,kint,l0,a0,del,segmentMass;

	// set l0 to 1
	l0 = 1.0;

	// set interaction distance to a fraction of l0
	del = 0.2*l0;

	// set a0 accordingly (will be units of l0^2 since l0 = 1)
	double asphericity,minAsphericity;
	asphericity = 1.5;
	a0 = (NV*NV*l0*l0)/(4*PI*asphericity);

	// set unit things to 1
	kl = 1.0;
	ka = 1.0;
	kb = 5.0;
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

	cout << "area = " << singleCell.area() << endl;
	cout << "a0 = " << a0 << endl;
	// return 0;

	// perturb initial vertex positions
	dscale = 0.1;
	cout << "randomizing by scale = " << dscale << endl;
	for (i=0; i<NV; i++){
		// get random perturbations
		dx = drand48();
		dy = drand48();

		// normalize perturbations
		dnorm = sqrt(dx*dx + dy*dy);
		dx /= dnorm;
		dy /= dnorm;

		// rescale dx by small scale
		dx *= dscale;
		dy *= dscale;

		// perturb x direction
		singleCell.setVPos(i,0,singleCell.vpos(i,0)+dx);

		// perturb y direction
		singleCell.setVPos(i,1,singleCell.vpos(i,1)+dy);
	}

	// update cell pos
	cout << "updating center of mass position accordingly" << endl;
	singleCell.updateCPos();

	// run energy minimization
	int NT = 10000;
	int NFRAMES = 1000;
	int NPRINT = floor(NT/(NFRAMES-1));
	int printCount = 0;
	double timeScale = sqrt((segmentMass*l0*l0)/kint);
	double dt = 0.025*timeScale;

	// use damping 
	double dampingParam = 0.01*sqrt(segmentMass*kint);

	// print initial configuration
	cout << "printing INITIAL vertex positions" << endl;

	// print sim details as header
	vertexPrintObject << setw(12) << left << "START" << " " << endl;
	vertexPrintObject << setw(12) << left << "NUMCL" << setw(6) << right << 1 << endl;
	vertexPrintObject << setw(12) << left << "NUMFR" << setw(6) << right << NFRAMES+1 << endl;
	vertexPrintObject << setw(12) << left << "BOXSZ" << setw(6) << right << L << endl;

	// print cell positions
	singleCell.printVertexPositions(vertexPrintObject,0,0);

	// print simulation header
	cout << endl << endl;
	cout << "================================================" << endl << endl;
	cout << "	Beginning MD to relax shape to minimum energy" << endl;
	cout << "		* NT = " << NT << endl;
	cout << "		* NPRINT = " << NPRINT << endl;
	cout << "		* segmentMass = " << segmentMass << endl;
	cout << " 		* timeScale = " << timeScale << endl;
	cout << "		* dampingParam = " << dampingParam << endl;
	cout << "		* dt = " << dt << endl << endl;
	cout << "================================================" << endl;
	cout << endl << endl;

	// initialize velocities
	// for (i=0; i<NV; i++){
	// 	for (d=0; d<2; d++)
	// 		singleCell.setVVel(i,d,0.1*drand48());
	// }
	singleCell.setCVel(0,1.0*(l0/timeScale));

	for (t=0; t<NT; t++){
		// update vertex positions
		singleCell.verletPositionUpdate(dt);

		// update cell pos
		singleCell.updateCPos();

		// update vertex forces
		singleCell.shapeForces();

		// update vertex velocities
		singleCell.verletVelocityUpdate(dt,dampingParam);

		// output if printint time
		if (t % NPRINT == 0){
			cout << "===================" << endl << endl;
			cout << " 	t = " << t << endl << endl;
			cout << "===================" << endl;

			cout << "	* Printing vetex positions to file" << endl;
			singleCell.printVertexPositions(vertexPrintObject,0,printCount+1);

			cout << "	* Printing cell energy to file" << endl;
			singleCell.printCellEnergy(energyPrintObject,printCount+1);
			cout << endl << endl;
			printCount++;
		}
	}
	vertexPrintObject << setw(12) << left << "STERM" << " " << endl;

	cout << "test code singleCell.cpp is finished, returning 0" << endl;
	vertexPrintObject.close();
	energyPrintObject.close();
	return 0;
}
