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
const double timeStepMag 		= 0.03;			// time step in MD unit
const double T0 				= 1e-8;			// initial velocities for read-in cells

// force parameters
const double ka 				= 1.0;			// area force constant (should be = 1)
const double gam 				= 0.0;			// surface tension force constant
const double kint 				= 0.01;			// interaction energy constant
const double a 					= 0.0;			// attraction parameter 
const double del 				= 1.0;			// radius of vertices in units of l0

// tolerances
const double Ftol 				= 1e-15;		// force tolerance (for FIRE min)
const double pscale 			= 0.8;			// scale to get new pressure

// ouputs
const double dlogpTol 			= 0.1;			// log difference between pressures of different frames

// int main
int main()
{
	// input variables
	int seed, NCOMP;
	double calA0, kl, kb, pTarget, Fcheck, Pcheck, dphi;

	// inputs from command line
	string inputFile 			= "examples/jamming/dpmb_dong.pos";
	string energyFile 			= "dongEn.test";
	string jammingFile 			= "dongPos.test";
	string vdosFile 			= "dongVDOS.test";

	// simulation values
	kl = 1.0;
	kb = 0.01;
	seed = 1;
	dphi = 5e-5;

	// number of compression states after state 1 (total is NCOMP + 1)
	NCOMP = 3;

	// instantiate main packing object
	cout << "	** Reading in from file " << inputFile << endl;
	cellPacking2D packingObject(inputFile,T0,seed);

	// set deformability, force values USING DONG'S FORCE COEFFICIENTS
	double kltmp, kbtmp, l0;
	int NCELLS = packingObject.getNCELLS();
	for (int ci=0; ci<NCELLS; ci++){
		kltmp = kl/packingObject.cell(ci).getNV();

		l0 = packingObject.cell(ci).getl0();
		calA0 = packingObject.cell(ci).calA0();
		kbtmp = ((l0*l0)/(calA0*calA0))*PI*kb;

		packingObject.cell(ci).setka(1.0);
		packingObject.cell(ci).setkl(kltmp);
		packingObject.cell(ci).setkb(kbtmp);
		packingObject.cell(ci).setkint(kint);

		cout << "	** ci = " << ci << ": ka = 1, kl = " << kltmp << ", kb = " << kbtmp << endl;
	}

	// update time scale
	packingObject.vertexDPMTimeScale(timeStepMag);

	// open jamming and vdos output file
	packingObject.openJamObject(jammingFile);
	packingObject.openStatObject(vdosFile);
	packingObject.openEnergyObject(energyFile);

	// set NT and NPRINT
	packingObject.setNT(NT);
	packingObject.setNPRINT(NPRINT);

	// compress to set packing fraction using FIRE, pressure relaxation
	cout << "	** STATE 0: relaxing system force to Ftol = " << Ftol << endl;
	packingObject.fireMinimizeF(Ftol, Fcheck, Pcheck);

	// compute pressure in initial jammed state
	cout << "	** After state 0, Fcheck = " << Fcheck << endl;
	cout << "	** After state 0, Pcheck = " << Pcheck << endl << endl;

	// print initial configuration and compute VDOS
	cout << "	** computing VDOS after state 1, printing to " << vdosFile << endl << endl;
	packingObject.vdos();

	// compute degree to which size should be increased
	double phi = packingObject.packingFraction();
	double rscale = sqrt((phi + dphi)/phi);

	// compress by dphi, measure VDOS
	for (int cc=0; cc<NCOMP; cc++){
		// decrease by packing fraction
		cout << "	** Compression protocol state " << cc + 1 << "/" << NCOMP+1 << " from phi = " << phi << " to phi + dphi = " << phi + dphi << " by dphi = " << dphi << endl;
		packingObject.scaleLengths(rscale);

		// minimize forces of new state
		cout << "	** STATE " << cc + 1 << "/" << NCOMP+1 << ": relaxing system force to Ftol = " << Ftol << endl;
		packingObject.fireMinimizeF(Ftol, Fcheck, Pcheck);

		cout << "	** After state " << cc + 1 << "/" << NCOMP+1 << ", Fcheck = " << Fcheck << endl;
		cout << "	** After state " << cc + 1 << "/" << NCOMP+1 << ", Pcheck = " << Pcheck << endl << endl;
		cout << "	** computing VDOS after state " << cc + 1 << "/" << NCOMP+1 << ", printing to " << vdosFile << endl << endl;
		packingObject.vdos();
	}


	cout << "	** FINISHED COMPRESSING ABOVE JAMMING, ENDING MAIN FILE" << endl;
	return 0;	
}