#ifndef DEFORMABLEPARTICLES2D_H
#define DEFORMABLEPARTICLES2D_H


/*

	-- DEFORMABLE PARTICLE CLASS -- 

	Class for individual deformable particles, 
	for use as member variables in an MD class

	Holds position, velocity and force information

	note: everything in d=2, not generalizable to any d

*/

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

class deformableParticles2D{
private:
	// spatial dimension
	int NDIM;

	// number of vertices
	int NV;

	// box length cell lives in (important for PBCs)
	double L;

	// vertex positions relative to center of mass
	double* vertexPositions;

	// vertex velocity in lab frame
	double* vertexVelocity;

	// vertex accerlation in lab frame
	double* vertexAcceleration;

	// forces on vertices
	double* vertexForces;

	// cell position in lab frame
	double* cellPosition;

	// interaction potential energy
	double* interactionPotential;

	// force parameters
	double kl;			// perimeter energy scale
	double ka;			// area energy scale
	double gam;			// surface tension energy scale
	double kb;			// bending energy scale
	double kint;		// interaction energy scale

	// rest parameters
	double l0;			// rest length for edges
	double a0;			// rest area for particles
	double del;			// contact distance for two edge segments
	double C;			// 1st attraction parameter (strength)
	double l;			// 2nd attraction parameter (distance)

public:
	// constructors
	deformableParticles2D();
	deformableParticles2D(int n);
	~deformableParticles2D();

	// operators
	void operator=(deformableParticles2D& onTheRight);

	// initialization
	void initializeVertices();
	void initializeCell();
	void regularPolygon();
	void regularPolygon(double inputArea);
	void vertexPerturbation(double dscale);

	// getters (simple access)
	int getNV() { return NV; };
	double getL() { return L; };
	double getkl() { return kl; };
	double getka() { return ka; };
	double getgam() { return gam; };
	double getkb() { return kb; };
	double getkint() { return kint; };
	double getl0() { return l0; };
	double geta0() { return a0; };
	double getdel() { return del; };
	double getC() { return C; };
	double getl() {return l; }

	double vpos(int vertex, int dim);
	double vrel(int vertex, int dim);
	double vvel(int vertex, int dim);
	double vacc(int vertex, int dim);
	double vforce(int vertex, int dim);
	double cpos(int dim);
	double cvel(int dim);
	double cforce(int dim);
	double uInt(int vertex);

	// setters (simple mutation)
	void setNV(int nv) { NV = nv; };
	void setL(double val) { L = val; };
	void setkl(double val) { kl = val; };
	void setka(double val) { ka = val; };
	void setgam(double val) { gam = val; };
	void setkb(double val) { kb = val; };
	void setkint(double val) { kint = val; };
	void setl0(double val) { l0 = val; };
	void seta0(double val) { a0 = val; };
	void setdel(double val) { del = val; };
	void setC(double val) { C = val; };
	void setl(double val) { l = val; };

	void setVPos(int vertex, int dim, double val);
	void setVRel(int vertex, int dim, double val);
	void setVVel(int vertex, int dim, double val);
	void setVAcc(int vertex, int dim, double val);
	void setVForce(int vertex, int dim, double val);
	void setCPos(int dim, double val);
	void setCVel(int dim, double val);
	void setCForce(int dim, double val);
	void setUInt(int vertex, double val);
	void setAsphericity(double val);

	// update cpos
	void updateCPos();

	// scale all lengths
	void scale(double val);

	// calculations
	double area(int vertex);
	double area();									// area of cell
	double perimeter();								// perimeter of cell
	double asphericity();							// asphericity parameter
	double segmentLength(int vertex);				// distance between vertex i and i+1
	double segment(int vertex, int dim);			// vector connecting i and i+1
	double dotProduct(int v1, int v2);				// dot product between vertices v1 and v2
	double segmentDotProduct(int l1, int l2);		// dot product between segments l1 and l2

	// force functions
	void shapeForces();
	void perimeterForce();
	void segmentOverlapForce();
	void areaForce();
	void surfaceTensionForce();
	void bendForce();
	int segmentForce(deformableParticles2D& onTheRight); // return 0 or 1, depending on contact
	int radialForce(deformableParticles2D& onTheRight); 

	// energy functions
	double perimeterEnergy();
	double areaEnergy();
	double surfaceTensionEnergy();
	double bendEnergy();
	double interactionEnergy();
	double totalPotentialEnergy();
	double totalKineticEnergy();

	// overloaded operators
	void operator% (const deformableParticles2D &onTheRight); 	// interaction force calculation operator

	// integrator options
	void verletPositionUpdate(double dt);
	void verletVelocityUpdate(double dt, double dampingParam);
	
	// print functions
	void printVertexPositions(std::ofstream& vertexPrintObject, int cellID);
	void printVertexPositions(std::ofstream& vertexPrintObject, int cellID, int frame);
	void printCellEnergy(std::ofstream& energyPrintObject, int frame);
};



#endif