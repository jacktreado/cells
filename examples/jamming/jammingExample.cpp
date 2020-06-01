/*

	Example .cpp file to compress NCELLS to jamming
	and to compute shear modulus

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
const int NPRINT 				= 5e2;			// number of time steps between prints
const double timeStepMag 		= 0.01;			// time step in MD unit
const double phiDisk 			= 0.4;			// initial phi of SP disks
const double deltaPhi0 			= 5e-4;			// initial delta phi
const double sizeRatio 			= 1.4;			// ratio between small and large particles
const double sizeFraction		= 0.5;			// fraction of small particles

// force parameters
const double ka 			= 1.0;			// area force constant (should be = 1)
const double gam 			= 0.0;			// surface tension force constant
const double kint 			= 0.1;			// interaction energy constant
const double a 				= 0.0;			// attraction parameter 
const double del 			= 1.0;			// radius of vertices in units of l0

// tolerances
const double Ftol 			= 1e-9;		// force tolerance (for FIRE min)
const double Ptol 			= 1e-8;			// pressure tolerance

// main function
int main()
{
	// local variables

	// output files
	string posFile = "pos.test";
	string enFile = "en.test";
	string jamFile = "jam.test";
	string vdosFile = "vdos.test";

	// system details
	int NCELLS 		= 6;
	int NV			= 12;
	int seed 		= 1;
	double Ltmp 	= 1.0;

	// mechanical parameters
	double kl = 1.0;
	double kb = 0.01;
	double calA0 = 1.1;

	// instantiate main packing object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,1.0,seed);

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initializing gel at phiDisk = " << phiDisk << " using SP model" << endl;
	packingObject.initializeBidisperse(NV, phiDisk, sizeRatio, sizeFraction, del);

	// set deformability, force values
	packingObject.forceVals(calA0,ka,kl,gam,kb,kint,del,a);

	// update time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// open output files
	packingObject.openPackingObject(posFile);
	packingObject.openEnergyObject(enFile);
	packingObject.openJamObject(jamFile);
	packingObject.openStatObject(vdosFile);

	// compress to set packing fraction using FIRE, pressure relaxation
	cout << "	** jamming protocol with Ftol = " << Ftol << ", Ptol = " << Ptol << endl;
	packingObject.findJamming(deltaPhi0, Ftol, Ptol);


	cout << "	** computing VDOS, printing to " << vdosFile << endl;
	packingObject.vdos();

	return 0;

	// // compute shear modulus
	// double shearStrain = 1e-2;
	// packingObject.setShearStrain(shearStrain);
	// double G = packingObject.shearModulus();

	// // print to console
	// cout << "shear modulus near jamming = " << G << endl;

	/*
	// loop over shear strains, compute shear modulus
	int nlogpts = 50;
	double p0 = 1e-10;
	double p1 = 1e-4;
	vector<double> strains(nlogpts,0.0);
	double dp = (log10(p1) - log10(p0))/(nlogpts - 1);
	double logp, linp;
	logp = log10(p0);
	strains.at(0) = p0;
	double G1,G2,G3,G20,G30;
	double currstrain, sxy, syx;

	// print G information to file
	ofstream shearData;
	shearData.open("shear.test");

	// loop over vector, populate with points
	for (int i=1; i<nlogpts; i++){
		logp = logp + dp;
		linp = pow(10.0,logp);
		strains.at(i) = linp;
	}

	// initial shear stresses of jammed state
	G20 = packingObject.getSigmaXY();
	G30 = packingObject.getSigmaYX();

	cout << endl << endl;
	cout << "Computing shear modulus using U, sigmaXY and sigmaYX..." << endl;
	for (int pp=0; pp<nlogpts; pp++){
		// update strain
		currstrain = strains.at(pp);

		// print to console
		cout << "	* * * setting shear strain step to dg = " << currstrain << endl;

		// compute shear modulus (method 1)
		;

		// compute shear modulus (methods 2/3)
		G2 = (sxy - G20)/currstrain;
		G3 = (syx - G30)/currstrain;
		cout << "sxy = " << sxy << ", syx = " << syx << endl;

		// output different shear modulus computations
		cout << G1 << " " << G2 << " " << G3 << " " << endl;
		shearData << currstrain << " " << G1 << " " << G2 << " " << G3 << " " << endl;
	}

	shearData.close();

	// print to console, end
	cout << endl << endl;
	cout << "	** FINISHED COMPRESSING TO JAMMING, ENDING MAIN FILE" << endl;
	*/
}
























