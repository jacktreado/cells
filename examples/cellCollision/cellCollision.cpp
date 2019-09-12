/*

	Test .cpp file to test forces and print functions
	on two cells colliding

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
	string vertexPositionsStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/cellCollision/cellCollisionPositions.dat";
	string energyStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/cellCollision/cellCollisionEnergy.dat";

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
	double L = 35.0;

	// Basic cell variables (cells will be identical)
	int NV = 30;
	double kl,ka,gam,kb,kint,l0,a0,del,C,l,segmentMass;
	double asphericity;
	double polyRad;

	// set l0 to 1
	l0 = 1.0;

	// set interaction distance to a fraction of l0
	del = 0.1*l0;

	// set attraction parameter to be 0, which is in units of delta
	C = 0.0;
	l = 0.0;

	// deformability
	asphericity = 1.06;
	a0 = (NV*NV*l0*l0)/(4*PI*asphericity);

	// size
	polyRad = sqrt((2*a0)/(NV*sin(2*PI/NV)));

	// using perimeter and area force
	kl = 1.0;
	ka = 1.0;
	kb = 2.0;

	// set a strong interaction parameter
	kint = 1.0;

	// NOT using surface tension and bending force
	gam = 0.0;

	// get segment mass
	segmentMass = l0*del + (PI*del*0.5*del*0.5);

	// instantiate 2 cell objects
	deformableParticles2D cell1(NV);
	deformableParticles2D cell2(NV);

	// set box for both cells
	cell1.setL(L);
	cell2.setL(L);

	// set initial position of cells to be on either side of box,
	// on slightly different lines so the cells bounce off of one another
	// and scatter
	cell1.setCPos(0,0.3*L);
	cell1.setCPos(1,0.55*L);

	cell2.setCPos(0,0.7*L);
	cell2.setCPos(1,0.45*L);

	// set values for dimensional quantities to set scales
	cell1.setl0(l0);
	cell1.seta0(a0);
	cell1.setdel(del);
	cell1.setC(C);
	cell1.setl(l);

	cell2.setl0(l0);
	cell2.seta0(a0);
	cell2.setdel(del);
	cell2.setC(C);
	cell2.setl(l);

	// set force constants
	cell1.setkl(kl);
	cell1.setka(ka);
	cell1.setgam(gam);
	cell1.setkb(kb);
	cell1.setkint(kint);

	cell2.setkl(kl);
	cell2.setka(ka);
	cell2.setgam(gam);
	cell2.setkb(kb);
	cell2.setkint(kint);

	// initialize both objects to be regular polygons
	cout << "making both cells regular polygons" << endl;
	cell1.regularPolygon();
	cell2.regularPolygon(a0);

	cout << "perturbing perimeter slightly" << endl;
	cell1.vertexPerturbation(0.4);
	cell2.vertexPerturbation(0.4);

	cout << "initial area = " << cell1.area() << endl;
	cout << "a0 = " << a0 << endl;

	// run cell collision event
	int NT = 10000;
	int NFRAMES = 1000;
	int NPRINT = floor(NT/(NFRAMES-1));
	int printCount = 0;
	int inContact = 0;
	double timeScale = sqrt((segmentMass*l0*l0)/kint);
	double dt = 0.02*timeScale;

	// use damping 
	double dampingParam = 0.0*sqrt(segmentMass*kint);

	// print initial configuration
	cout << "printing INITIAL vertex positions" << endl;

	// print sim details as header
	vertexPrintObject << setw(12) << left << "START" << " " << endl;
	vertexPrintObject << setw(12) << left << "NUMCL" << setw(6) << right << 2 << endl;
	vertexPrintObject << setw(12) << left << "NUMFR" << setw(6) << right << NFRAMES+1 << endl;
	vertexPrintObject << setw(12) << left << "BOXSZ" << setw(6) << right << L << endl;

	// print cell positions
	cell1.printVertexPositions(vertexPrintObject,0,0);
	cell2.printVertexPositions(vertexPrintObject,1);

	// initialize velocities for collision
	double vscale = 0.2*(l0/timeScale);

	// initialize velocities for shape relaxation
	for (i=0; i<NV; i++){
		for (d=0; d<2; d++){
			cell1.setVVel(i,d,0.075*drand48());
			cell2.setVVel(i,d,0.075*drand48());
		}
	}

	// print simulation header
	cout << endl << endl;
	cout << "================================================" << endl << endl;
	cout << "	Beginning MD to relax shape to minimum energy" << endl;
	cout << "		* NT = " << NT << endl;
	cout << "		* NPRINT = " << NPRINT << endl;
	cout << "		* segmentMass = " << segmentMass << endl;
	cout << " 		* timeScale = " << timeScale << endl;
	cout << "		* dampingParam = " << dampingParam << endl;
	cout << "		* dt = " << dt << endl;
	cout << "		* collision velocity scale = " << vscale << endl << endl;
	cout << "================================================" << endl;
	cout << endl << endl;

	// loop over time to relax shape
	for (t=0; t<NT; t++){
		// update vertex positions
		cell1.verletPositionUpdate(dt);
		cell2.verletPositionUpdate(dt);

		// update cell pos
		cell1.updateCPos();
		cell2.updateCPos();

		// give velocity
		if (t == 20){
			cell1.setCVel(0,vscale);
			cell1.setCVel(1,0.0);
			cell2.setCVel(0,-vscale);
			cell2.setCVel(1,0.0);
		}

		// update vertex forces
		cell1.shapeForces();
		cell2.shapeForces();

		// calculate interaction
		inContact = cell1.segmentForce(cell2);
		
		// update vertex velocities
		cell1.verletVelocityUpdate(dt,dampingParam);
		cell2.verletVelocityUpdate(dt,dampingParam);

		// output if printint time
		if (t % NPRINT == 0){
			cout << "===================" << endl << endl;
			cout << " 	t = " << t << endl << endl;
			cout << "===================" << endl;

			cout << "	* inContact = " << inContact << endl;
			cout << "	* Printing vetex positions to file" << endl;
			cell1.printVertexPositions(vertexPrintObject,0,printCount+1);
			cell2.printVertexPositions(vertexPrintObject,1);

			cout << "	* Printing cell energy to file" << endl;
			energyPrintObject << cell1.totalKineticEnergy() + cell2.totalKineticEnergy() << " " << cell1.totalPotentialEnergy() + cell2.totalPotentialEnergy() << endl;
			cout << endl << endl;
			printCount++;
		}
	}

	vertexPrintObject << setw(12) << left << "STERM" << " " << endl;

	cout << "test code cellCollision.cpp is finished, returning 0" << endl;
	vertexPrintObject.close();
	energyPrintObject.close();

	return 0;
}