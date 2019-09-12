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
	double L;						// boundary length
	double T;						// system temperature
	double phi;						// system packing/volume fraction
	
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
	void initializeInteractionParams(double del, double C, double l);

	// initialize positions
	void squareLattice();
	void hexagonalLattice();	// NEEDS WORKS

	// initialize velocities
	void initializeVelocities(double tmp0);

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
	double getL() { return L; };
	double getT() { return T; };
	double getphi() { return phi; };

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
	double meanAsphericity();				// calculate system-average asphericity

	// setters
	void setNCELLS(int nc) { NCELLS = nc; };
	void setNT(int nt) { NT = nt; };
	void setNPRINT(int nprint) { NPRINT = nprint; };
	void setL(double val) { L = val; };
	void setT(double val) { T = val; };

	void setdt(double val) { dt0 = val*timeScale(); dt = dt0; };
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

	// printers
	void printSystemPositions(int frame);
	void printSystemEnergy(int frame, double Uval, double Kval);
	void printSystemStats();




	/**************************

		Simulation Functions

	***************************/


	// looping functions
	void msFire(double dphi0, double dphiJ, double phiJGuess, double Ktol, double Ftol);
	void msDamping(double dphi, double Ktol, double Ftol, double dampingParameter);
	int potentialRelaxFire(double kineticTol, double potentialTol, int plotIt, int& frameCount);
	void jammingDamping(double dphi, double Ktol, double Utol, double dampingParameter);	// NEEDS WORKS
	void jammingFire(double dphi, double Ktol, double Utol);	// NEEDS WORKS
	void jammingFireRamp(double dphi, double dCalA, double asphericityTarget, double kbTarget, double Ktol, double Utol);

	// non-equilibrium MD functions
	void isoExtensionQS(double phiTarget, double dphi, double dampingParameter);


	// wrapper functions for simulation protocols
	void calculateForces();
	void dverlet(double dampingParameter);
	void fverlet(int& np, double& alpha, double dampingParameter);
	void activeBrownian(double diffusionConstant);	// NEEDS WORKS

	// relaxation
	void shapeRelax(int plotIt);
	void shapeRelaxRamp(double finalCalA, double initialDCalA, double kineticTol, double potentialTol);
	void overlapRelief();

	// fire energy minimization
	void fireStep(int& np, double& alpha);
};

#endif