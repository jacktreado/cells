#ifndef CELLPACKING2D_H
#define CELLPACKING2D_H


/*

	-- PACKING OF CELLS CLASS -- 

	* Class for systems of deformable particles, 
		holds deformable particles in packing.

	* Mostly wrapper functions for objects in cells array

	* Information on temperature, boundary changes (extension, compression, shear),

	* Note: everything in d=2, not generalizable to any d

*/

#include "deformableParticles2D.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

class cellPacking2D{
private:

	// int scalars
	int NDIM;						// spatial dimension (will always be 2)
	int NCELLS;						// number of cells
	int NT;							// number of total time steps
	int NPRINT;						// number of time steps between print steps

	// double scalars
	double seed;					// initial seed
	double dt;						// time step size
	double dt0;						// initial time step size (for dynamic time step methods)
	double T;						// system temperature
	double phi;						// system packing/volume fraction

	// boundary lengths
	std::vector<double> L;

	// virial stresses
	double sigmaXX, sigmaXY, sigmaYX, sigmaYY;
	
	// array of cells
	deformableParticles2D* cellArray;

	// contact matrix
	int* contactMatrix;				// array of contacts between particles

	// ofstream objects
	std::ofstream packingPrintObject;
	std::ofstream energyPrintObject;
	std::ofstream statPrintObject;

public:

	// Constructors and Destructors
	void defaultvars();
	cellPacking2D();
	cellPacking2D(int ncells, int nt, int nprint, double l, double s);
	cellPacking2D(int ncells, int ntumor, int tumorNV, int adiposeNV, double tumorCalA, double adiposeCalA, int s);
	cellPacking2D(std::ifstream& inputFileObject, double asphericity, double s);
	~cellPacking2D();

	// operators
	void operator=(cellPacking2D& onTheRight);	// assign one configuration to another object

	// initialize sizes
	void initializeMonodisperse(int NV, double asphericity);
	void initializeBidisperse(int NV, double sizeRatio);
	void initializeBidisperse(int NV, double sizeRatio, double asphericity);
	void initializePolydisperse(int NV, double sizeVariance, double asphericity, char distChoice);	// NEEDS WORKS

	// initialize parameters
	void initializeForceConstants(double kl, double ka, double gam, double kb, double kint);
	void initializeInteractionParams(double del, double a);

	// initialize positions
	void squareLattice();
	void hexagonalLattice();	// NEEDS WORKS

	// initialize velocities
	void initializeVelocities(double tmp0);
	void initializeVelocities(int ci, double tmp0);

	// file openers
	void openPackingObject(std::string& str){
		packingPrintObject.open(str.c_str());
		if (!packingPrintObject.is_open()) {
			std::cout << "	ERROR: packingPrintObject could not open " << str << "..." << std::endl;
			exit(1);
		}
	}

	void openEnergyObject(std::string& str){
		energyPrintObject.open(str.c_str());
		if (!energyPrintObject.is_open()) {
			std::cout << "	ERROR: energyPrintObject could not open " << str << "..." << std::endl;
			exit(1);
		}
	}

	void openStatObject(std::string& str){
		statPrintObject.open(str.c_str());
		if (!statPrintObject.is_open()) {
			std::cout << "	ERROR: statPrintObject could not open " << str << "..." << std::endl;
			exit(1);
		}
	}

	// read in initial configuration from file
	void readInFromFile(std::ifstream& inputFileObject);

	// getters
	int getNCELLS() { return NCELLS; };
	int getNT() { return NT; };
	int getNPRINT() { return NPRINT; };
	double getdt() { return dt; };
	double getT() { return T; };
	double getphi() { return phi; };

	// box len
	double getL(int d) { return L.at(d); };

	deformableParticles2D& cell(int ci);	// return cell object ci
	int nframes();							// number of frames in the simulation
	int cmindex(int ci, int cj);			// contact matrix index
	int contacts(int ci, int cj);			// contact matrix element
	int totalContacts();

	// system-wide calculations
	int totalNumberOfContacts();			// calculate total number of contacts
	double timeScale();						// MD characteristic time scale
	double packingFraction();				// calculate system packing fraction
	double shapePotentialEnergy();			// calculate potential energy due to shape
	double relaxPotentialEnergy();			// potential energy - bending energy
	double totalPotentialEnergy();			// calculate total potential energy
	double interactionPotentialEnergy();	// calculate only interaction potential energy
	double totalKineticEnergy();			// calculate total kinetic energy
	double maxForceMagnitude();				// calculate magnitude of largest force
	double maxNetForceMagnitude();			// calculate magnitude of largest net force on any cell
	double forceRMS(); 						// RMS net force
	double meanAsphericity();				// calculate system-average asphericity

	// setters
	void setNCELLS(int nc) { NCELLS = nc; };
	void setNT(int nt) { NT = nt; };
	void setNPRINT(int nprint) { NPRINT = nprint; };
	void setT(double val) { T = val; };

	// box len
	void setL(int d, double val) { L.at(d) = val; };

	// set dt: NOTE, IN NON-DIMENSIONAL FORM, so mass = rho pi r^2 = pi
	void setdt(double val) { dt0 = val*sqrt(4*atan(1)); dt = dt0; };
	void addContact(int ci, int cj);
	void deleteContact(int ci, int cj);
	void resetContacts();
	int particleContacts(int ci);
	void setPackingFraction(double val);
	void setAsphericity(double val);
	void setAsphericity(int ci, double val);
	void scaleLengths(double val);
	void rescaleVelocities(double temperature);
	int removeRattlers(int krcrs);

	/**************************

		Forces and position 
			updates

	***************************/

	void calculateForces();
	void gelationForces();
	void fverlet(int& np, double& alpha, double dampingParameter);
	void activeBrownian(double diffusionConstant);	// NEEDS WORKS



	/**************************

		Simulation Functions

	***************************/

	// NVE test function
	void cellNVE();
	void cellOverDamped();

	// looping functions
	void jammingFireRamp(double dphi, double dCalA, double asphericityTarget, double kbTarget, double phiTarget, double Ktol, double Ptol, int plotIt);
	void compressToTarget(double dphi, double phiTarget, double asphericityTarget, double Ktol, double Ptol, int plotIt, int& frameCount);

	// FIRE 2.0 relaxation functions
	void fireMinimizeP(double Ptol, double Ktol);
	void fireMinimizeF(double Ftol, double Ktol, int plotIt, int& frameCount);

	// Gelation functions
	void twoParticleContact(int NV);
	void initializeGel(int NV, double phiDisk, double sizeDispersion, double delval);
	void gelForceVals(double calA0, double kl, double ka, double gam, double kb, double kint, double del, double a);
	void qsIsoCompression(double phiTarget, double deltaPhi);
	void attractionRamp(double attractionTarget, double dAttraction);
	void gelRateExtension(double phiGel, double gelRate, double timeStepMag);
	void gelRK4();

	// Sticky SP particle functions
	void initializeStickySP(std::vector<double>& radii, double phiDisk, double sizeDispersion);
	void stickySPTriangularLattice(std::vector<double>& radii, double phiDisk);
	void stickySPGelationQS(std::vector<double>& radii, double phiGel, double dphiGel, double attractiveParam);
	void stickySPGelationRate(std::vector<double>& radii, double phiGel, double gelRate, double attractiveParam, double timeStepMag);
	void printPositionsStickySP(std::vector<double>& radii);
	void printEnergyStickySP();

	// Repulsive SP particle functions
	void fireMinimizeSP(std::vector<double>& lenscales);
	void fireMinimizeSP(std::vector<double>& radii, double attractiveParam);
	void spForces(std::vector<double>& lenscales);
	void spAttractiveForces(std::vector<double>& radii, double attractiveParam);
	void spPosVerlet();
	void spVelVerlet(std::vector<double>& lenscales);
	void spNVE(std::vector<double>& lenscales, int nt);

	// active pipeflow functions
	void initializeActiveStickySP(std::vector<double>& radii, int NV, double phiDisk, double sizeDispersion, double Lscale);
	void singleActiveCell(int NV, double phiInit, double calA0, double Dc, double Dv, double tv, double v0);
	void spActivePipeForces(std::vector<double>& radii);
	void spActivePipeWallForces(std::vector<double>& radii);
	void spActivePipeNVE(std::vector<double>& radii, double T0);
	void spActivePipeFlow(std::vector<double>& radii, double a, double v0, double Dr);

	// non-equilibrium MD functions
	void isoExtensionQS(int plotIt, int& frameCount, double phiTarget, double dphi);


	// Hopper functions
	void initializeHopperSP(std::vector<double>& radii, double w0, double w, double th, double Lmin, int NV);
	void fireMinimizeHopperSP(std::vector<double>& radii, double w0, double w, double th);
	void fireMinimizeHopperDP(double w0, double w, double th);
	void hopperForcesSP(std::vector<double>& radii, double w0, double w, double th, double g, int closed);
	void hopperForcesDP(double w0, double w, double th, double g, int closed);
	void hopperWallForcesSP(std::vector<double>& radii, double w0, double w, double th, int closed);
	void hopperWallForcesDP(double w0, double w, double th, int closed);
	void hopperSPNVE(std::vector<double>& radii, double w0, double w, double th, double T0);
	void flowHopperSP(std::vector<double>& radii, double w0, double w, double th, double g);
	void flowHopperDP(double w0, double w, double th, double g);
	double hopperPackingFraction(std::vector<double>& radii, double w0, double w, double th);
	void setHopperPackingFraction(std::vector<double>& radii, double phiNew, double w0, double w, double th);
	void hopperPosVerletSP();
	void hopperVelVerletSP(std::vector<double>& radii);
	void printHopperSP(std::vector<double>& radii, double w0, double w, double th, double g);
	void printHopperDP(double w0, double w, double th, double g);


	// tumor MD functions
	void tumorForce(int NTUMORCELLS, double forceScale, double adiposeDamping);


	// fire energy minimization
	void fireStep(int& np, double& alpha);
	void fireStep(std::vector<int>& np, std::vector<double>& alphaVec);


	// printers
	void printSystemPositions();
	void printSystemEnergy();
	void printSystemEnergy(int kmin);
	void printSystemContacts();
	void printSystemPositions(int frame);
	void printSystemEnergy(int frame, double Uval, double Kval);
	void printSystemStats();

};

#endif