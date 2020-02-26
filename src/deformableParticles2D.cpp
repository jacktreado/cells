/*

	Methods file for deformableParticles2D class

	IDEAS FOR POTENTIAL SPEED UPS
		* need to check L every time calculating a segment? or can I use relative vectors?
		* in segmentForce(), determine optimal theta_c for edge feasibility
		* in segmentForce(), determine optimal buffer for cell-cell contact checking
		* add a neighbor list
		* make GPU friendly (more for sim class)? Is this possible in c++ is class uses objects of another class as member variables?
*/

// include file
#include "deformableParticles2D.h"


// namespace
using namespace std;


// constants
const double PI = 4*atan(1);



/************************

	Constructors

*************************/


// default constructor
deformableParticles2D::deformableParticles2D(){
	// dimensionality ALWAYS set to 2
	NDIM 	= 2;

	// scalar variables set to 0
	NV 		= 0;
	kl 		= 0.0;
	ka 		= 0.0;
	gam 	= 0.0;
	kb 		= 0.0;
	kint 	= 0.0;
	l0 		= 0.0;
	a0 		= 0.0;
	c0 		= cos(0.0);
	del 	= 0.0;
	a 		= 0.0;

	// pointer variables point to nullptr
	vertexPositions 		= nullptr;
	vertexVelocity	 		= nullptr;
	vertexAcceleration 		= nullptr;
	vertexForces 			= nullptr;
	cellPosition 			= nullptr;
	interactionPotential 	= nullptr;

	// periodic boundary conditions set to on
	pbc.resize(NDIM);
	for (int d=0; d<NDIM; d++)
		pbc.at(d) = 1;

	// box lengths set to 1.0
	L.resize(NDIM);
	for (int d=0; d<NDIM; d++)
		L.at(d) = 1.0;
}

// constructor to specify number of vertices
deformableParticles2D::deformableParticles2D(int n){
	// dimensionality ALWAYS set to 2
	NDIM 	= 2;

	// scalar variables set to 0
	NV 		= 0;
	kl 		= 0.0;
	ka 		= 0.0;
	gam 	= 0.0;
	kb 		= 0.0;
	kint 	= 0.0;
	l0 		= 0.0;
	a0 		= 0.0;
	del 	= 0.0;
	a 		= 0.0;

	// box lengths set to 1.0
	L.resize(NDIM);
	for (int d=0; d<NDIM; d++)
		L.at(d) = 1.0;

	// pointer variables point to nullptr
	vertexPositions 		= nullptr;
	vertexVelocity	 		= nullptr;
	vertexAcceleration 		= nullptr;
	vertexForces 			= nullptr;
	cellPosition 			= nullptr;
	interactionPotential 	= nullptr;

	// set number of vertices
	if (n <= 2){
		cout << "	ERROR: input nv = " << n << ", which needs to be at least 3. Ending." << endl;
		exit(1);
	}
	NV = n;

	// initialize vertices
	initializeVertices();

	// initialize cells
	initializeCell();
}

// destructor
deformableParticles2D::~deformableParticles2D(){
	if (vertexPositions){
		delete [] vertexPositions;
		vertexPositions = nullptr;
	}
	if (vertexVelocity){
		delete [] vertexVelocity;
		vertexVelocity = nullptr;
	}
	if (vertexAcceleration){
		delete [] vertexAcceleration;
		vertexAcceleration = nullptr;
	}
	if (vertexForces){
		delete [] vertexForces;
		vertexForces = nullptr;
	}
	if (cellPosition){
		delete [] cellPosition;
		cellPosition = nullptr;
	}
	if (interactionPotential){
		delete [] interactionPotential;
		interactionPotential = nullptr;
	}
}




/************************

	Operators

*************************/


// overloaded assignment operator (ASSUME ARRAYS IN THIS HAVE BEEN INITIALIZED!)
void deformableParticles2D::operator=(deformableParticles2D& onTheRight){
	// local variables
	int i,d;

	// dimensionality ALWAYS set to 2
	NDIM 	= 2;

	// check that array pointers do not point to null
	if (!vertexPositions){
		cout << "	ERROR: in operator=, memory for vertexPositions not yet allocated. Ending." << endl;
		exit(1);
	}
	else if (!vertexVelocity){
		cout << "	ERROR: in operator=, memory for vertexVelocity not yet allocated. Ending." << endl;
		exit(1);
	}
	else if (!vertexAcceleration){
		cout << "	ERROR: in operator=, memory for vertexAcceleration not yet allocated. Ending." << endl;
		exit(1);
	}
	else if (!vertexForces){
		cout << "	ERROR: in operator=, memory for vertexForces not yet allocated. Ending." << endl;
		exit(1);
	}
	else if (!interactionPotential){
		cout << "	ERROR: in operator=, memory for interactionPotential not yet allocated. Ending." << endl;
		exit(1);
	}
	else if (!cellPosition){
		cout << "	ERROR: in operator=, memory for cellPosition not yet allocated. Ending." << endl;
		exit(1);
	}


	// copy scalar member variables
	NV 		= onTheRight.NV;
	kl 		= onTheRight.kl;
	ka 		= onTheRight.ka;
	gam 	= onTheRight.gam;
	kb 		= onTheRight.kb;
	kint 	= onTheRight.kint;
	a0 		= onTheRight.a0;
	l0 		= onTheRight.l0;
	del 	= onTheRight.del;
	a 		= onTheRight.a;

	// set box lengths equal
	for (int d=0; d<NDIM; d++){
		pbc.at(d) = onTheRight.pbc.at(d);
		L.at(d) = onTheRight.L.at(d);
	}

	// deep copy vertex values
	for (i=0; i<NV; i++){
		for (d=0; d<NDIM; d++){
			setVPos(i,d,onTheRight.vpos(i,d));
			setVVel(i,d,onTheRight.vvel(i,d));
			setVAcc(i,d,onTheRight.vacc(i,d));
			setVForce(i,d,onTheRight.vforce(i,d));
		}		
		setUInt(i,onTheRight.uInt(i));
	}

	// deep copy com pos
	for (d=0; d<NDIM; d++)
		setCPos(d,onTheRight.cpos(d));
}


  

/************************

	Initialization

*************************/

// initialize vertex arrays
void deformableParticles2D::initializeVertices(){
	// local variables
	int i,d;

	// check if NV has been set > 0
	if (NV <= 0){
		cout << "	ERROR: in initializeVertices(), NV = " << NV << ", so not set properly. Ending." << endl;
		exit(1);
	}

	// check if memory has already been allocated
	if (vertexPositions){
		cout << "	ERROR: in initializeVertices(), memory for vertexPositions already allocated. Ending." << endl;
		exit(1);
	}
	else if (vertexVelocity){
		cout << "	ERROR: in initializeVertices(), memory for vertexVelocity already allocated. Ending." << endl;
		exit(1);
	}
	else if (vertexAcceleration){
		cout << "	ERROR: in initializeVertices(), memory for vertexAcceleration already allocated. Ending." << endl;
		exit(1);
	}
	else if (vertexForces){
		cout << "	ERROR: in initializeVertices(), memory for vertexForces already allocated. Ending." << endl;
		exit(1);
	}
	else if (interactionPotential){
		cout << "	ERROR: in initializeVertices(), memory for interactionPotential already allocated. Ending." << endl;
		exit(1);
	}


	// allocate memory
	vertexPositions = new double[NDIM*NV];
	vertexForces = new double[NDIM*NV];
	vertexVelocity = new double[NDIM*NV];
	vertexAcceleration = new double[NDIM*NV];
	interactionPotential = new double[NV];

	// set equal to 0
	for (i=0; i<NV; i++){
		setUInt(i,0.0);
		for (d=0; d<NDIM; d++){
			setVPos(i,d,0.0);
			setVVel(i,d,0.0);
			setVAcc(i,d,0.0);
			setVForce(i,d,0.0);
		}
	}
}

// initialize cell arrays
void deformableParticles2D::initializeCell(){
	// local variables
	int d;

	// check if memory has already been allocated
	if (cellPosition){
		cout << "	ERROR: in initializeCell(), memory for cellPosition already allocated. Ending." << endl;
		exit(1);
	}

	// allocate memory
	cellPosition = new double[NDIM];

	// set equal to 0
	for (d=0; d<NDIM; d++)
		setCPos(d,0.0);
}

// initialize vertex positions so cell begins as regular polygon
void deformableParticles2D::regularPolygon(){
	// local variables
	int i;
	double angleArg = 0.0;
	double polyRad;

	// check if NV has been set > 0
	if (NV <= 0){
		cout << "	ERROR: in regularPolygon(), NV = " << NV << ", so not set properly. Ending." << endl;
		exit(1);
	}
	else if (a0 < 0.1){
		cout << "	ERROR: in regularPolygon(), a0 = " << a0 << ", so too small and not set properly. Ending." << endl;
		exit(1);
	}

	// set radius of polygon
	polyRad = sqrt((2.0*a0)/(NV*sin(2.0*PI/NV)));

	// loop over vertices, set positions using rotations
	for (i=0; i<NV; i++){
		angleArg = (2.0*PI*i)/NV;
		setVRel(i,0,-polyRad*sin(angleArg));
		setVRel(i,1,polyRad*cos(angleArg));
	}
	// output
	cout << " 	-- creating regular polygon with a0 = " << a0 << ", area = " << polygonArea() << " and perimeter = " << perimeter() << ", so init calA0 = " << pow(perimeter(),2.0)/(4.0*PI*polygonArea()) << ", compare to " << NV*tan(PI/NV)/PI << endl;
}

// initialize vertex positions so cell begins as regular polygon
void deformableParticles2D::regularPolygon(double inputArea){
	// local variables
	int i;
	double angleArg = 0.0;
	double polyRad;

	// set a0 to inputArea
	a0 = inputArea;

	// check if NV has been set > 0
	if (NV <= 0){
		cout << "	ERROR: in regularPolygon(), NV = " << NV << ", so not set properly. Ending." << endl;
		exit(1);
	}
	else if (inputArea < 0.1){
		cout << "	ERROR: in regularPolygon(), inputArea = " << inputArea << ", so too small and not set properly. Ending." << endl;
		exit(1);
	}

	// set radius of polygon
	polyRad = sqrt((2*inputArea)/(NV*sin(2*PI/NV)));

	// loop over vertices, set positions using rotations
	for (i=0; i<NV; i++){
		angleArg = (2*PI*i)/NV;
		setVRel(i,0,-polyRad*sin(angleArg));
		setVRel(i,1,polyRad*cos(angleArg));
	}
}

// perturb vertex positions
void deformableParticles2D::vertexPerturbation(double dscale){
	// local variables
	int i;
	double dx,dy,dnorm;
	
	// loop over vertices, perturb
	for (i=0; i<NV; i++){
		// get random perturbations
		dx = 1.0 - 2.0*drand48();
		dy = 1.0 - 2.0*drand48();

		// normalize perturbations
		dnorm = sqrt(dx*dx + dy*dy);
		dx /= dnorm;
		dy /= dnorm;

		// rescale dx by small scale and l0
		dx *= dscale*l0;
		dy *= dscale*l0;

		// perturb x direction
		setVPos(i,0,vpos(i,0)+dx);

		// perturb y direction
		setVPos(i,1,vpos(i,1)+dy);
	}
}



/************************

	Getters

*************************/

// vertex position in lab frame (stored in vertexPositions)
double deformableParticles2D::vpos(int vertex, int dim){
	// check inputs
	if (vertex >= NV){
		cout << "	ERROR: vertex input in vpos = " << vertex << ", which is >= NV, which is NV = " << NV << ". Ending." << endl;
		exit(1);
	}
	else if (dim >= NDIM){
		cout << "	ERROR: dim input in vpos = " << dim << ", which is >= NDIM. Ending." << endl;
		exit(1);
	}

	// local variable
	int index = NDIM*vertex + dim;

	// return value
	return vertexPositions[index];
}


// vertex position relative to center of mass 
double deformableParticles2D::vrel(int vertex, int dim){
	// check inputs
	if (vertex >= NV){
		cout << "	ERROR: vertex input in vpos = " << vertex << ", which is >= NV. Ending." << endl;
		exit(1);
	}
	else if (dim >= NDIM){
		cout << "	ERROR: dim input in vpos = " << dim << ", which is >= NDIM. Ending." << endl;
		exit(1);
	}

	// local variable
	int index = NDIM*vertex + dim;
	double relDist;

	// get relative distance with MIC
	relDist = distance(vertexPositions[index],cellPosition[dim],dim);

	// return value
	return relDist;
}


// vertex velocity in lab frame
double deformableParticles2D::vvel(int vertex, int dim){
	// check input
	if (vertex >= NV){
		cout << "	ERROR: vertex input in vvel = " << vertex << ", which is >= NV. Ending." << endl;
		exit(1);
	}
	else if (dim >= NDIM){
		cout << "	ERROR: dim input in vvel = " << dim << ", which is >= NDIM. Ending." << endl;
		exit(1);
	}

	// local variable
	int index = NDIM*vertex + dim;

	// return value
	return vertexVelocity[index];
}


// vertex acceleration in lab frame
double deformableParticles2D::vacc(int vertex, int dim){
	// check input
	if (vertex >= NV){
		cout << "	ERROR: vertex input in vacc = " << vertex << ", which is >= NV. Ending." << endl;
		exit(1);
	}
	else if (dim >= NDIM){
		cout << "	ERROR: dim input in vacc = " << dim << ", which is >= NDIM. Ending." << endl;
		exit(1);
	}

	// local variable
	int index = NDIM*vertex + dim;

	// return value
	return vertexAcceleration[index];
}


// net force in lab frame
double deformableParticles2D::vforce(int vertex, int dim){
	// check inputs
	if (vertex >= NV){
		cout << "	ERROR: vertex input in vforce = " << vertex << ", which is >= NV. Ending." << endl;
		exit(1);
	}
	else if (dim >= NDIM){
		cout << "	ERROR: dim input in vforce = " << dim << ", which is >= NDIM. Ending." << endl;
		exit(1);
	}

	// local variable
	int index = NDIM*vertex + dim;

	// return value
	return vertexForces[index];
}



// cell center of mass position
double deformableParticles2D::cpos(int dim){
	// check input
	if (dim >= NDIM){
		cout << "	ERROR: dim input in cpos = " << dim << ", which is >= NDIM. Ending." << endl;
		exit(1);
	}

	// return value
	return cellPosition[dim];
}


// cell center of mass velocity
double deformableParticles2D::cvel(int dim){
	// check input
	if (dim >= NDIM){
		cout << "	ERROR: dim input in cvel = " << dim << ", which is >= NDIM. Ending." << endl;
		exit(1);
	}

	// local variables
	int i;
	double velSum = 0.0;

	// loop over vertices
	for (i=0; i<NV; i++)
		velSum += vvel(i,dim);

	return velSum/NV;
}


// net force on the cell center of mass
double deformableParticles2D::cforce(int dim){
	// check input
	if (dim >= NDIM){
		cout << "	ERROR: dim input in cforce = " << dim << ", which is >= NDIM. Ending." << endl;
		exit(1);
	}

	// local variables
	int i;
	double forceSum = 0.0;

	// loop over vertices
	for (i=0; i<NV; i++)
		forceSum += vforce(i,dim);

	// return value
	return forceSum;
}


// interaction potential energy per vertex (energy on segment starting at vertex i)
double deformableParticles2D::uInt(int vertex){
	// check inputs
	if (vertex >= NV){
		cout << "	ERROR: vertex input in uInt = " << vertex << ", which is >= NV. Ending." << endl;
		exit(1);
	}

	// return value
	return interactionPotential[vertex];
}


// distance between points
double deformableParticles2D::distance(double p2, double p1, int d){
	// local variables
	double dp;

	// calculate distance between points
	dp = p2 - p1;

	// if pbc, check minimum image distance
	if (pbc.at(d) == 1)
		dp -= L.at(d)*round(dp/L.at(d));

	// return value
	return dp;
}



// distance between two cells: points from this to cell onTheRight
double deformableParticles2D::cellDistance(deformableParticles2D& onTheRight, int d){
	// local variables
	double distComp;

	// calculate normal distance component between centers
	distComp = onTheRight.cpos(d) - cpos(d);

	// if pbcs are on, check minimum image distance
	if (pbc.at(d) == 1)
		distComp -= L.at(d)*round(distComp/L.at(d));

	// return distance component
	return distComp;
}


/************************

	Setters

*************************/


// set position of vertex in lab frame
void deformableParticles2D::setVPos(int vertex, int dim, double val){
	// check inputs
	if (vertex >= NV){
		cout << "	ERROR: vertex input in setVPos = " << vertex << ", which is >= NV. Ending." << endl;
		exit(1);
	}
	else if (dim >= NDIM){
		cout << "	ERROR: dim input in setVPos = " << dim << ", which is >= NDIM. Ending." << endl;
		exit(1);
	}

	// local variable
	int index = NDIM*vertex + dim;

	// set value
	vertexPositions[index] = val;
}


// set position of vertex in cell frame
void deformableParticles2D::setVRel(int vertex, int dim, double val){
	// check inputs
	if (vertex >= NV){
		cout << "	ERROR: vertex input in setVRel = " << vertex << ", which is >= NV. Ending." << endl;
		exit(1);
	}
	else if (dim >= NDIM){
		cout << "	ERROR: dim input in setVRel = " << dim << ", which is >= NDIM. Ending." << endl;
		exit(1);
	}

	// local variable
	int index = NDIM*vertex + dim;

	// set value
	vertexPositions[index] = val + cpos(dim);
}


// set velocity of vertices
void deformableParticles2D::setVVel(int vertex, int dim, double val){
	// check input
	if (vertex >= NV){
		cout << "	ERROR: vertex input in setVVel = " << vertex << ", which is >= NV. Ending." << endl;
		exit(1);
	}
	else if (dim >= NDIM){
		cout << "	ERROR: dim input in setVVel = " << dim << ", which is >= NDIM. Ending." << endl;
		exit(1);
	}

	// local variable
	int index = NDIM*vertex + dim;

	// set value
	vertexVelocity[index] = val;
}


// set acceleration on cell vertex
void deformableParticles2D::setVAcc(int vertex, int dim, double val){
	// check input
	if (vertex >= NV){
		cout << "	ERROR: vertex input in setVAcc = " << vertex << ", which is >= NV. Ending." << endl;
		exit(1);
	}
	else if (dim >= NDIM){
		cout << "	ERROR: dim input in setVAcc = " << dim << ", which is >= NDIM. Ending." << endl;
		exit(1);
	}

	// local variable
	int index = NDIM*vertex + dim;

	// set value
	vertexAcceleration[index] = val;
}


// set force on cell vertex
void deformableParticles2D::setVForce(int vertex, int dim, double val){
	// check inputs
	if (vertex >= NV){
		cout << "	ERROR: vertex input in setVForce = " << vertex << ", which is >= NV. Ending." << endl;
		exit(1);
	}
	else if (dim >= NDIM){
		cout << "	ERROR: dim input in setVForce = " << dim << ", which is >= NDIM. Ending." << endl;
		exit(1);
	}

	// local variable
	int index = NDIM*vertex + dim;

	// set value
	vertexForces[index] = val;
}


// set position of cell center of mass
void deformableParticles2D::setCPos(int dim, double val){
	// check input
	if (dim >= NDIM){
		cout << "	ERROR: dim input in setCPos = " << dim << ", which is >= NDIM. Ending." << endl;
		exit(1);
	}

	// set value
	cellPosition[dim] = val;
}


// impose a global velocity on the entirety of the cell by dividing velocity between each vertex
void deformableParticles2D::setCVel(int dim, double val){
	// local variables
	int i;

	// loop over vertices, divide velocity val between each vertex
	for (i=0; i<NV; i++)
		setVVel(i,dim,val);

}


// impose a net force on the entirety of the cell by dividing force between each vertex
void deformableParticles2D::setCForce(int dim, double val){
	// local variables
	int i;

	// loop over vertices, divide force val between each vertex
	for (i=0; i<NV; i++)
		setVForce(i,dim,val/NV);
}


// set interaction potential energy
void deformableParticles2D::setUInt(int vertex, double val){
	// check inputs
	if (vertex >= NV){
		cout << "	ERROR: vertex input in setUInt = " << vertex << ", which is >= NV. Ending." << endl;
		exit(1);
	}

	// set value
	interactionPotential[vertex] = val;
}


// set asphericity by changing a0
void deformableParticles2D::setAsphericity(double val){
	if (val < 1.0){
		cout << "	ERROR: trying to set asphericity to be < 1, ending." << endl;
		exit(1);
	}

	// set area based on asphericity
	a0 = pow(NV*l0,2.0)/(4.0*PI*val);
}

// set asphericity by changing l0
void deformableParticles2D::setAsphericityConstA(double val){
	if (val < 1.0){
		cout << "	ERROR: trying to set asphericity to be < 1, ending." << endl;
		exit(1);
	}

	l0 = (1.0/NV)*sqrt(4*PI*a0*val);
}


// update cpos based on vpos
void deformableParticles2D::updateCPos(){
	// local variables
	int i,d;
	double cposx = 0.0;
	double cposy = 0.0;

	// loop over vertices
	for (i=0; i<NV; i++){
		cposx += vrel(i,0) + cpos(0);
		cposy += vrel(i,1) + cpos(1);
	}

	// divide by NV
	cposx /= NV;
	cposy /= NV;

	// check PBCs
	if (pbc.at(0) == 1){
		if (cposx > L.at(0))
			cposx -= L.at(0);
		else if (cposx < 0)
			cposx += L.at(0);
	}
	
	if (pbc.at(1) == 1){
		if (cposy > L.at(1))
			cposy -= L.at(1);
		else if (cposy < 0)
			cposy += L.at(1);
	}

	// divide cpos by NV to get centroid
	setCPos(0,cposx);
	setCPos(1,cposy);
}

// scale all lengths by factor
void deformableParticles2D::scale(double val){
	// local variables
	int i,d;

	// loop over vertex positions, scale
	for (i=0; i<NV; i++){
		for (d=0; d<NDIM; d++)
			setVRel(i,d,vrel(i,d)*val);
	}

	// scale force-related dimensional parameters
	l0 *= val;
	a0 *= pow(val,NDIM);
}



/************************

	Calculations

*************************/


// calculate polygon area
double deformableParticles2D::polygonArea(){
	// local variables
	int i, kstart, km1, kp1, kcurr, knext;
	bool cwpos;
	double maxy = 0;
	double totalArea = 0.0;

	// find index of maximum y value
	for (i=0; i<NV; i++){
		if (vrel(i,1) > maxy){
			maxy = vrel(i,1);
			kstart = i;
		}
	}

	// get indices left and right of kstart
	km1 	= (kstart + NV - 1) % NV;
	kp1 	= (kstart + 1) % NV;

	// determine which way is clockwise
	if (vrel(kp1,0) > vrel(km1,0))
		cwpos = true;
	else
		cwpos = false;

	// initialize indices for counting
	kcurr = kstart;

	// loop over vertices
	for (i=0; i<NV; i++){
		// determine index of next vertex
		if (cwpos)
			knext = (kcurr + 1) % NV;
		else
			knext = (kcurr + NV - 1) % NV;

		// add to Dong's area formula
		totalArea += vrel(knext,0)*vrel(kcurr,1) - vrel(kcurr,0)*vrel(knext,1);

		// update k
		kcurr = knext;
	}

	// halve area
	totalArea *= 0.5;

	// check area sign
	if (totalArea < 0)
		totalArea *= -1;

	// return area
	return totalArea;
}

// calculate cell area
// NOTE: OLD VERSION USED TRIANGULAR AREA FROM area(i) FUNCTION
double deformableParticles2D::area(){
	// local variables
	int vi, vim1;
	double li, lim1, ci;
	double uxi, uyi, uxim1, uyim1;
	double thetai, cvxi;
	double avi, avtot;
	double av = (PI*pow(0.5*del*l0,2))/(2.0*PI);

	// loop over vertices, add correct area fraction
	avtot = 0.0;
	for (vi=0; vi<NV; vi++){
		// get previous vertex
		vim1 = (vi - 1 + NV) % NV;

		// check cosine
		ci = segmentCosine(vi);

		// if parallel segments, just take half vert area, else calc arc area
		if (abs(ci - 1.0) < 1e-8)
			avi = PI*av;
		else{
			// get arc fraction for given vertex
			thetai = acos(ci);

			// get segment lengths
			li = segmentLength(vi);
			lim1 = segmentLength(vim1);

			// get components of unit vectors
			uxi = segment(vi,0)/li;
			uyi = segment(vi,1)/li;
			uxim1 = segment(vim1,0)/lim1;
			uyim1 = segment(vim1,1)/lim1;

			// check convexity
			cvxi = uxim1*uyi - uxi*uyim1;

			// calculate fraction for vertex i
			if (cvxi > 0.0)
				avi = av*(PI + thetai);
			else
				avi = av*(PI - thetai);

			// add to area
			avtot += avi;
		}

		// check for error
		if (avi != avi){
			cout << "	ERROR: vertex area calculated is a nan, vi = " << vi << ", thetai = " << thetai << ". " << endl;
			cout << "	** uxi = " << uxi << ", uyi = " << uyi << "; uxim1 = " << uyim1 << ", uyim1 = " << uyim1 << endl;
			cout << " 	** segmentCosine(vi) = " << segmentCosine(vi) << endl;
			cout << "  	** Ending." << endl;
			exit(1);
		}
	}

	// return area
	return polygonArea() + avtot;
}

// calculate cell perimeter
double deformableParticles2D::perimeter(){
	// local variables
	int i;
	double totalPerimeter = 0.0;

	// loop over segments, add segment length to perimeter
	for (i=0; i<NV; i++)
		totalPerimeter += segmentLength(i);

	// return perimeter
	return totalPerimeter;
}


// calculate instantaneous asphericity
double deformableParticles2D::asphericity(){
	return pow(perimeter(),2)/(4*PI*polygonArea());
}


// calculate preferred asphericity
double deformableParticles2D::calA0(){
	return pow(NV*l0,2.0)/(4.0*PI*a0);
}


// calculate segment length between vertex+1 and vertex
double deformableParticles2D::segmentLength(int vertex){
	// check inputs
	if (vertex >= NV){
		cout << "	ERROR: vertex input in segmentLength = " << vertex << ", which is >= NV. Ending." << endl;
		exit(1);
	}

	// local variables
	int d;
	double segLength = 0.0;

	// get segment length
	for (d=0; d<NDIM; d++)
		segLength += pow(segment(vertex,d),2);

	// take square root
	segLength = sqrt(segLength);

	// return segment length
	return segLength;
}


// get vectorial component of segment between vertex+1 and vertex
double deformableParticles2D::segment(int vertex, int dim){
	// check inputs
	if (vertex >= NV){
		cout << "	ERROR: vertex input in segment = " << vertex << ", which is >= NV. Ending." << endl;
		exit(1);
	}
	else if (dim >= NDIM){
		cout << "	ERROR: vertex input in segment = " << dim << ", which is >= NDIM. Ending." << endl;
		exit(1);
	}

	// local variables
	int ip1;
	double seg;

	// wrap vertex labels
	ip1 = (vertex+1) % NV;

	// check minimum image
	seg = distance(vpos(ip1,dim),vpos(vertex,dim),dim);

	// return dim component of segment vector
	return seg;
}


// get dot product between vertex v1 and v2
double deformableParticles2D::dotProduct(int v1, int v2){
	// local variables
	int d;
	double val = 0.0;

	// loop over dim
	for (d=0; d<NDIM; d++)
		val += vrel(v1,d)*vrel(v2,d);

	// return value
	return val;
}


// get dot product between segments l1 and l2
double deformableParticles2D::segmentDotProduct(int l1, int l2){
	// local variables
	int l1p1,l2p1;
	double val;

	// wrap vertex labels
	l1p1 = (l1+1) % NV;
	l2p1 = (l2+1) % NV;

	// calc by vertex dot products
	val = dotProduct(l1p1,l2p1) - dotProduct(l1,l2p1) - dotProduct(l2,l1p1) + dotProduct(l1,l2);

	// return value
	return val;
}

double deformableParticles2D::segmentCosine(int vi){
	// local variables
	int vim1, d;
	double dp, li, lim1, ui, uim1;

	// wrap label of vi - 1
	vim1 = (vi-1+NV) % NV;

	// get segment length
	li = segmentLength(vi);
	lim1 = segmentLength(vim1);

	// calculate dot product of unit vectors
	dp = 0.0;
	for (d=0; d<NDIM; d++){
		ui = segment(vi,d)/li;
		uim1 = segment(vim1,d)/lim1;
		dp += ui*uim1;
	}

	// return normalized dot product
	return dp;
}



/************************

	Shape Forces

*************************/

void deformableParticles2D::shapeForces(){
	// calculate perimeter force on vertex i
	if (kl > 0)
		perimeterForce();

	// calculate perimeter force on vertex i
	if (ka > 0)
		areaForce();

	// calculate perimeter force on vertex i
	if (gam > 0)
		surfaceTensionForce();

	// calculate perimeter force on vertex i
	if (kb > 0)
		bendForce();
}


// RECALCULATE SHAPE FORCES FROM BALANCED ENERGY!!!
void deformableParticles2D::balancedShapeForces(){
	// local variables                                                                                                                                                                                                                 
	int d,i,im1,im2,ip1;
	double ftmp,fxTmp,fyTmp;
	double K,Kl,Ka,Kb;
	double aStrain,lStrainI,lStrainIm1;
	double Cim1,Ci,sign_i,sign_im1;

	// bending energy parameters
	double lim2,lim1,li,lip1; 					// segment lengths
	double cim1,ci,cip1;						// angle cosines
	double uim2,uim1,ui,uip1;					// segment unit vector elements
	double cStrainIp1,cStrainI,cStrainIm1;		// cosine strains
	double normC0;								// cosine strain renormalization

	// total area
	double totalArea = polygonArea();
	// c0 = 1.0;

	// loop over vertices, calculate each force that is active
	for (i=0; i<NV; i++){
		// wrap vertices
		im2 = (i-2+NV) % NV;
		im1 = (i-1+NV) % NV;
		ip1 = (i+1) % NV;

		// calculate segment lengths
		lim1 = segmentLength(im1);
		li = segmentLength(i);

		// calculate perimeter force
		if (kl > 0){
			// get constants                                                                                                                                                                                                           
	        Cim1 = 1.0/lim1;
	        Ci = 1.0/li;

	        // get segment strains
	        lStrainI = (li/l0) - 1.0;
	        lStrainIm1 = (lim1/l0) - 1.0;

	        // get effective segment length spring constant
	        Kl = kl*NV*l0;

	        // loop over dimensions, add to force                                                                                                                                                                                      
	        for (d=0; d<NDIM; d++){
	            ftmp = Ci*lStrainI*segment(i,d) - Cim1*lStrainIm1*segment(im1,d);
	            ftmp *= Kl;
	            setVForce(i,d,vforce(i,d)+ftmp);
	        }
		}

		// calculate area force
		if (ka > 0){
			// effective area spring constant
			Ka = (ka*a0/(l0*l0));

			// calculate area strain
			aStrain = (totalArea/a0) - 1.0;

			// calculate force term in each direction (based on calc from notes)
			fxTmp = -Ka*aStrain*0.5*(vrel(ip1,1) - vrel(im1,1));
			fyTmp = -Ka*aStrain*0.5*(vrel(im1,0) - vrel(ip1,0));

			// add to force on vertices
			setVForce(i,0,vforce(i,0)+fxTmp);
			setVForce(i,1,vforce(i,1)+fyTmp);
		}

		// calculate bending force
		if (kb > 0){
			// define segment lengths
			lim2 = segmentLength(im2);
			lip1 = segmentLength(ip1);

			// define cosines
			cim1 	= segmentCosine(im1);
			ci 		= segmentCosine(i);
			cip1 	= segmentCosine(ip1);

			// cosine normalization factor
			normC0 = 1 + c0;

			// effective bending spring
			Kb = (kb*NV*l0*l0)/normC0;

			// define cosine strains
			cStrainIm1 = (cim1 - c0)/normC0;
			cStrainI = (ci - c0)/normC0;
			cStrainIp1 = (cip1 - c0)/normC0;

			// add to force in each direction
			for(d=0; d<NDIM; d++){
				// define unit vector segment elements
				uim2 = segment(im2,d)/lim2;
				uim1 = segment(im1,d)/lim1;
				ui = segment(i,d)/li;
				uip1 = segment(ip1,d)/lip1;

				// get force element
				ftmp = cStrainIm1*(uim2 - cim1*uim1)/lim1;
				ftmp += cStrainIp1*(cip1*ui - uip1)/li;
				ftmp += cStrainI*(((ui - ci*uim1)/lim1) + ((ci*ui - uim1)/li));
				ftmp *= -Kb;
				
				// add to vectorial force
				setVForce(i,d,vforce(i,d)+ftmp);
			}
		}
	}

	// calculate perimeter force on vertex i
	if (gam > 0)
		surfaceTensionForce();
}

void deformableParticles2D::perimeterForce(){
	// local variables                                                                                                                                                                                                                 
	int d,i,im1;
	double Cim1,Ci,ftmp;

	// loop over vertices, calculate forces                                                                                                                                                                                            
    for (i=0; i<NV; i++){
        // wrap vertices                                                                                                                                                                                                           
        im1 = (i-1+NV) % NV;

        // get constants                                                                                                                                                                                                           
        Cim1 = 1.0-(l0/segmentLength(im1));
        Ci = 1.0-(l0/segmentLength(i));

        // loop over dimensions, add to force                                                                                                                                                                                      
        for (d=0; d<NDIM; d++){
            ftmp = Ci*segment(i,d) - Cim1*segment(im1,d);
            ftmp *= kl;
            setVForce(i,d,vforce(i,d)+ftmp);
        }
    }
}

void deformableParticles2D::areaForce(){
	// local variables
	int i,ip1,im1,d;
	double sign_i, sign_im1;
	double totalArea = polygonArea();
	double fxTmp, fyTmp; // scalar force


	// loop over vertices
	for (i=0; i<NV; i++){
		// determine ip1,im1 vertices
		im1 = (i-1+NV) % NV;
		ip1 = (i+1) % NV;

		// get sign of grad terms
		sign_i = vrel(i,0)*vrel(ip1,1) - vrel(ip1,0)*vrel(i,1);
		if (sign_i < 0)
			sign_i = -1.0;
		else
			sign_i = 1.0;

		sign_im1 = vrel(im1,0)*vrel(i,1) - vrel(i,0)*vrel(im1,1);
		if (sign_im1 < 0)
			sign_im1 = -1.0;
		else
			sign_im1 = 1.0;

		// calculate force term in each direction (based on calc from notes)
		fxTmp = -ka*0.5*(totalArea-a0)*(sign_i*vrel(ip1,1) - sign_im1*vrel(im1,1));
		fyTmp = -ka*0.5*(totalArea-a0)*(sign_im1*vrel(im1,0) - sign_i*vrel(ip1,0));

		// add to force on vertices
		setVForce(i,0,vforce(i,0)+fxTmp);
		setVForce(i,1,vforce(i,1)+fyTmp);
	}
}

void deformableParticles2D::surfaceTensionForce(){
	// local variables
	int i,im1,d;
	double iUnitVectorComponent,im1UnitVectorComponent,ftmp;

	// loop over vertices
	for (i=0; i<NV; i++){
		// make sure loop is cyclic
		im1 = (i-1+NV) % NV;

		// add to force
		for (d=0; d<NDIM; d++){
			// get unit vector components
			im1UnitVectorComponent = segment(im1,d)/segmentLength(im1);
			iUnitVectorComponent = segment(i,d)/segmentLength(i);
			ftmp = -gam*(iUnitVectorComponent - im1UnitVectorComponent);

			// add to force
			setVForce(i,d,vforce(i,d)+ftmp);
		}
	} 
}

void deformableParticles2D::bendForce(){
	// local variables
	int i,ip1,im1,im2,d;				// indices
	double lim2,lim1,li,lip1; 			// segment lengths
	double Kim2, Kim1, Ki, Kip1;		// big K constants: scale segment unit vectors
	double kim2_im1,kim1,kim1_i,ki;		// little K constants: define big K constants
	double cim2,cim1,ci;				// angle cosines
	double ftmp;						// scalar component of force

	// loop over vertices
	for (i=0; i<NV; i++){
		// wrap vertex labels
		im2 = (i-2+NV) % NV;
		im1 = (i-1+NV) % NV;
		ip1 = (i+1) % NV;
		
		// define segment lengths
		lim2 = segmentLength(im2);
		lim1 = segmentLength(im1);
		li = segmentLength(i);
		lip1 = segmentLength(ip1);

		// define cosines
		cim2 = segmentDotProduct(im2,im1)/(lim2*lim1);
		cim1 = segmentDotProduct(im1,i)/(lim1*li);
		ci = segmentDotProduct(i,ip1)/(li*lip1);

		// define little k constants
		kim2_im1 = (cim2-c0)/lim1;
		kim1 = (cim1-c0)/lim1;
		kim1_i = (cim1-c0)/li;
		ki = (ci-c0)/li;

		// define big K constants
		Kim2 = kim2_im1;
		Kim1 = Kim2*cim2 + kim1*cim1 + kim1_i;
		Ki = ki*ci + kim1 + kim1_i*cim1;
		Kip1 = ki;

		// add to force in each direction
		for(d=0; d<NDIM; d++){
			// get ftmp
			ftmp = Kip1*segment(ip1,d)/lip1 - Ki*segment(i,d)/li;
			ftmp += Kim1*segment(im1,d)/lim1 - Kim2*segment(im2,d)/lim2;
			ftmp *= kb;
			
			// add to vectorial force
			setVForce(i,d,vforce(i,d)+ftmp);
		}
	}
}



/************************

	Interaction Forces

	* Using L as a private member variable works, 
		but doing simulations with fixed boundary will be trickier

*************************/


// force between segments between vertices
// 		* interacting parts are circulolines, need to calculate dmin between two line segments
// 		* dmin calculation strategy outlined in deformableParticleForces.pdf
// 		* also updates interaction potential
// 		* need to pass in box length to properly check image distances
// 		* !! NOTE POSSIBLE BUG: if L in this and L' in onTheRight are NOT equal, 
// 			then PBCs are meaningless and errors will 
int deformableParticles2D::segmentForce(deformableParticles2D &onTheRight){
	// return variable
	int inContact = 0;

	// local variables
	int i,j,d;

	// -------------------------
	// 
	// 	   Distance cutoff
	//
	// -------------------------

	// section variables
	double centerDistance = 0.0; 			// center-to-center distance
	vector<double> deltaMuNu(NDIM,0.0);		// vector to store center-to-center distance vector
	double muREff,nuREff,buffer;			// variables to determine effect cell contact distance
	double distTmp = 0.0;					// temporary distance, for mimimum image convention (MIC)


	// calculate connecting vector deltaMuNu
	for (d=0; d<NDIM; d++){
		// get distance
		distTmp = cellDistance(onTheRight,d);

		// save distance in vector
		deltaMuNu.at(d) = distTmp;

		// calculate vector norm
		centerDistance += pow(distTmp,2);
	}
	centerDistance = sqrt(centerDistance);

	// get effect radii
	muREff = sqrt(area()/PI);
	nuREff = sqrt(onTheRight.area()/PI);
	buffer = 0.1*perimeter();

	// if not close enough, return 0
	if ((muREff + nuREff + 0.5*(del*l0 + onTheRight.del*onTheRight.l0) + buffer) < centerDistance)
		return 0;

	// -------------------------
	// 
	// 	   Feasible edges
	//
	// -------------------------

	// section variables
	vector<int> possibleMuEdges;			// set of feasible edges on cell mu (this)
	vector<int> possibleNuEdges;			// set of feasible edges on cell nu (onTheRight)
	double vertexCenterDotProduct;			// dot product between vertex location and center-to-center vector
	double dpTol = 0.2;						// cutoff for dot product between vertex location and center-to-center vector

	// initialize force vectors here
	vector< vector<double> > muForce;
	vector< vector<double> > nuForce;
	vector<double> tmpForceVec(NDIM,0.0);


	// loop over vertices i on mu (this) and check dot products with deltaMuNu to determine feasibility
	for (i=0; i<NV; i++){
		// initialize dot product between r_mu,i and deltaMuNu to be 0
		vertexCenterDotProduct = 0.0;

		// calculate dot product 
		for (d=0; d<NDIM; d++)
			vertexCenterDotProduct += vrel(i,d)*deltaMuNu.at(d);

		// check sign (+ means feasible contact)
		if (vertexCenterDotProduct > dpTol)
			possibleMuEdges.push_back(i);

		// push back entries to vector for force calculation
		muForce.push_back(tmpForceVec);
	}

	// loop over vertices j on nu (onTheRight) and check dot products with deltaMuNu to determine feasibility
	for (j=0; j<onTheRight.getNV(); j++){
		// initialize dot product between r_nu,j and deltaMuNu to be 0
		vertexCenterDotProduct = 0.0;

		// calculate dot product 
		for (d=0; d<NDIM; d++)
			vertexCenterDotProduct += onTheRight.vrel(j,d)*deltaMuNu.at(d);

		// check sign (- means feasible contact)
		if (vertexCenterDotProduct < -dpTol)
			possibleNuEdges.push_back(j);

		// push back entries to vector for force calculation
		nuForce.push_back(tmpForceVec);
	}

	// ----------------------------------------
	// 
	// 	  Feasible contacts / Force calculation
	//
	// ----------------------------------------

	cout << "	*** ERROR: segment force under construction at this time, need to incorporate better contact";
	cout << "				checking. Ending program." << endl;
	exit(1);

	// add forces to mu vertices
	for (i=0; i<NV; i++){
		for (d=0; d<NDIM; d++){
			setVForce(i,d,vforce(i,d) + muForce.at(i).at(d));
		}
	}

	// add forces to nu vertices
	for (j=0; j<onTheRight.getNV(); j++){
		for (d=0; d<NDIM; d++){
			onTheRight.setVForce(j,d,onTheRight.vforce(j,d) + nuForce.at(j).at(d));
		}
	}

	// return if or if not in contact
	return inContact;
}


// force between vertices
// 		* interacting parts are disks of diameter delta = l_0
// 		* also updates interaction potential

int deformableParticles2D::vertexForce(deformableParticles2D &onTheRight, vector<double>& fij, vector<double>& rij){
	// return variable
	int inContact = 0;

	// local variables
	int i,j,d,dd;

	// -------------------------
	// 
	// 	   Distance cutoff
	//
	// -------------------------

	// section variables
	double centerDistance = 0.0; 			// center-to-center distance
	vector<double> deltaMuNu(NDIM,0.0);		// vector to store center-to-center distance vector
	double muREff,nuREff,buffer;			// variables to determine effect cell contact distance
	double distTmp = 0.0;					// temporary distance, for mimimum image convention (MIC)

	// calculate connecting vector deltaMuNu
	for (d=0; d<NDIM; d++){
		// get distance
		distTmp = cellDistance(onTheRight,d);

		// save distance in vector
		rij.at(d) = -distTmp;
		deltaMuNu.at(d) = distTmp;

		// calculate vector norm
		centerDistance += pow(distTmp,2);
	}
	centerDistance = sqrt(centerDistance);

	// get effect radii
	muREff = sqrt(area()/PI);
	nuREff = sqrt(onTheRight.area()/PI);
	buffer = 0.1*perimeter();

	// if not close enough, return 0
	if ((muREff + nuREff + 0.5*(del*l0 + onTheRight.del*onTheRight.l0) + buffer) < centerDistance)
		return 0;


	// -------------------------
	// 
	// 	   Vertex forces
	//
	// -------------------------

	double forceScale = kint;				// force scale
	double energyScale = forceScale;		// energy scale
	double distScale = 0.0;					// distance scale
	double p1 = 1.0 + a;					// edge of interaction zone, units of delta
	double ftmp = 0.0;						// temporary force variable
	double uTmp = 0.0;						// temporary energy variable
	double vertexDist = 0.0;				// distance variable
	vector<double> vertexVec(NDIM,0.0);		// vector to hold vectorial distance quantity
	double contactDistance = 0.0;			// contact distance variable

	// loop over vertex pairs, check for contact
	for (i=0; i<NV; i++){
		for (j=0; j<onTheRight.NV; j++){

			// get distance between vertices i and j
			vertexDist = 0.0;
			for (d=0; d<NDIM; d++){
				// get distance to nearest image
				distTmp = distance(onTheRight.vpos(j,d),vpos(i,d),d);

				// add to vertex distance
				vertexVec.at(d) = distTmp;

				// add to scalar distance
				vertexDist += distTmp*distTmp;
			}

			// get contact distance
			contactDistance = 0.5*(del*l0 + onTheRight.del*onTheRight.l0);

			// get vertex distance
			vertexDist = sqrt(vertexDist);

			// check overlap distances
			if (vertexDist < contactDistance*p1){
				// set inContact to 1 for return
				inContact = 1;

				// define scaled distance (x = distance/contact distance)
				distScale = vertexDist/contactDistance;

				// update force and energy scales
				forceScale = kint;
				energyScale = kint * contactDistance;
				// energyScale = forceScale;

				// IF in zone to use repulsive force (and, if a > 0, bottom of attractive well)
				if (vertexDist < contactDistance){

					// add to interaction potential
					uTmp = 0.5 * energyScale * pow(1 - distScale,2) - (energyScale*a*a)/6.0;
		
					setUInt(i,uInt(i) + 0.5*uTmp);
					onTheRight.setUInt(j,onTheRight.uInt(j) + 0.5*uTmp);

					// add to vectorial forces
					for (d=0; d<NDIM; d++){
						// get force value
						ftmp = -forceScale * (1 - distScale) * vertexVec.at(d) / vertexDist;

						// add to force on i
						setVForce(i,d,vforce(i,d) + ftmp);

						// subtract off complement from force on j
						onTheRight.setVForce(j,d,onTheRight.vforce(j,d) - ftmp);

						// calculate contribution from vertex-vertex interaction
						// for (dd=0; dd<NDIM; dd++)
							// virialStress.at(NDIM*d + dd) += ftmp*vertexVec.at(dd);
						fij.at(d) += ftmp;
					}
				}

				// IF a > 0, in attractive well
				else if (vertexDist >= contactDistance && vertexDist < contactDistance*p1 && a > 0.0){

					// add to interaction potential
					uTmp =  (-energyScale/(6.0*a)) * pow(distScale - 1.0,2) * (2*(distScale - 1.0) - 3*a) - (energyScale*a*a)/6.0;
		
					setUInt(i,uInt(i) + 0.5*uTmp);
					onTheRight.setUInt(j,onTheRight.uInt(j) + 0.5*uTmp);

					// add to vectorial forces
					for (d=0; d<NDIM; d++){
						// get force value
						ftmp = -(forceScale/a) * (distScale - 1.0) * (distScale - p1) * vertexVec.at(d) / vertexDist;

						// add to force on i
						setVForce(i,d,vforce(i,d) + ftmp);

						// subtract off complement from force on j
						onTheRight.setVForce(j,d,onTheRight.vforce(j,d) - ftmp);

						// calculate contribution from vertex-vertex interaction
						fij.at(d) += ftmp;
					}
				}
			}
		}
	}

	// return if in contact or not
	return inContact;
}

int deformableParticles2D::vertexForce(deformableParticles2D &onTheRight, vector<double>& fij, vector<double>& rij, double aij){
	// return variable
	int inContact = 0;

	// local variables
	int i,j,d,dd;

	// -------------------------
	// 
	// 	   Distance cutoff
	//
	// -------------------------

	// section variables
	double centerDistance = 0.0; 			// center-to-center distance
	vector<double> deltaMuNu(NDIM,0.0);		// vector to store center-to-center distance vector
	double muREff,nuREff,buffer;			// variables to determine effect cell contact distance
	double distTmp = 0.0;					// temporary distance, for mimimum image convention (MIC)

	// calculate connecting vector deltaMuNu
	for (d=0; d<NDIM; d++){
		// get distance
		distTmp = cellDistance(onTheRight,d);

		// save distance in vector
		rij.at(d) = -distTmp;
		deltaMuNu.at(d) = distTmp;

		// calculate vector norm
		centerDistance += pow(distTmp,2);
	}
	centerDistance = sqrt(centerDistance);

	// get effect radii
	muREff = sqrt(area()/PI);
	nuREff = sqrt(onTheRight.area()/PI);
	buffer = 0.1*perimeter();

	// if not close enough, return 0
	if ((muREff + nuREff + 0.5*(del*l0 + onTheRight.del*onTheRight.l0) + buffer) < centerDistance)
		return 0;


	// -------------------------
	// 
	// 	   Vertex forces
	//
	// -------------------------

	double forceScale = kint;					// force scale
	double energyScale = forceScale;			// energy scale
	double distScale = 0.0;						// distance scale
	double l1 = aij;							// magnitude of interaction
	double l2 = a;								// interaction range
	double p1 = 1.0 + aij;						// point where attraction is max (in units of delta)
	double p2 = 1.0 + a;						// point where attractive zone ends (in units of delta)
	double ftmp = 0.0;							// temporary force variable
	double uTmp = 0.0;							// temporary energy variable
	double vertexDist = 0.0;					// distance variable
	vector<double> vertexVec(NDIM,0.0);			// vector to hold vectorial distance quantity
	double contactDistance = 0.0;				// contact distance variable

	// loop over vertex pairs, check for contact
	for (i=0; i<NV; i++){
		for (j=0; j<onTheRight.NV; j++){

			// get distance between vertices i and j
			vertexDist = 0.0;
			for (d=0; d<NDIM; d++){
				// get distance to nearest image
				distTmp = distance(onTheRight.vpos(j,d),vpos(i,d),d);

				// add to vertex distance
				vertexVec.at(d) = distTmp;

				// add to scalar distance
				vertexDist += distTmp*distTmp;
			}

			// get contact distance
			contactDistance = 0.5*(del*l0 + onTheRight.del*onTheRight.l0);

			// get vertex distance
			vertexDist = sqrt(vertexDist);

			// check overlap distances
			if (vertexDist < contactDistance*p2){
				// set inContact to 1 for return
				inContact = 1;

				// define scaled distance (x = distance/contact distance)
				distScale = vertexDist/contactDistance;

				// update force and energy scales
				forceScale = kint;
				energyScale = kint * contactDistance;

				// IF in zone to use repulsive force (and, if aij > 0, bottom of attractive well)
				if (vertexDist < contactDistance*p1){

					// add to interaction potential
					uTmp = 0.5 * energyScale * (pow(1 - distScale,2) - l1*l2);
		
					setUInt(i,uInt(i) + 0.5*uTmp);
					onTheRight.setUInt(j,onTheRight.uInt(j) + 0.5*uTmp);

					// add to vectorial forces
					for (d=0; d<NDIM; d++){
						// get force value
						ftmp = -forceScale * (1 - distScale) * vertexVec.at(d) / vertexDist;

						// add to force on i
						setVForce(i,d,vforce(i,d) + ftmp);

						// subtract off complement from force on j
						onTheRight.setVForce(j,d,onTheRight.vforce(j,d) - ftmp);

						// calculate force contribution from vertex-vertex interaction
						fij.at(d) += ftmp;
					}
				}

				// IF aij > 0, in attractive well
				else if (vertexDist >= contactDistance*p1){

					// add to interaction potential
					uTmp =  -(0.5*energyScale*l1/(l2 - l1)) * pow(distScale - 1 - l2,2);
		
					setUInt(i,uInt(i) + 0.5*uTmp);
					onTheRight.setUInt(j,onTheRight.uInt(j) + 0.5*uTmp);

					// add to vectorial forces
					for (d=0; d<NDIM; d++){
						// get force value
						ftmp = -(forceScale*l1/(l2 - l1)) * (distScale - 1.0 - l2) * vertexVec.at(d) / vertexDist;

						// add to force on i
						setVForce(i,d,vforce(i,d) + ftmp);

						// subtract off complement from force on j
						onTheRight.setVForce(j,d,onTheRight.vforce(j,d) - ftmp);

						// calculate contribution from vertex-vertex interaction
						fij.at(d) += ftmp;
					}
				}
			}
		}
	}

	// return if in contact or not
	return inContact;
}


// check how many vertex-vertex attractive contacts exist between two cells
int deformableParticles2D::pwAttractiveContacts(deformableParticles2D &onTheRight){
	// return variable
	int numAttractiveContacts = 0;

	// local variables
	int i,j,d,dd;

	// -------------------------
	// 
	// 	   Distance cutoff
	//
	// -------------------------

	// section variables
	double centerDistance = 0.0; 			// center-to-center distance
	vector<double> deltaMuNu(NDIM,0.0);		// vector to store center-to-center distance vector
	double muREff,nuREff,buffer;			// variables to determine effect cell contact distance
	double distTmp = 0.0;					// temporary distance, for mimimum image convention (MIC)

	// calculate connecting vector deltaMuNu
	for (d=0; d<NDIM; d++){
		// get distance
		distTmp = cellDistance(onTheRight,d);

		// save distance in vector
		deltaMuNu.at(d) = distTmp;

		// calculate vector norm
		centerDistance += pow(distTmp,2);
	}
	centerDistance = sqrt(centerDistance);

	// get effect radii
	muREff = sqrt(area()/PI);
	nuREff = sqrt(onTheRight.area()/PI);
	buffer = 0.1*perimeter();

	// if not close enough, return 0
	if ((muREff + nuREff + 0.5*(del*l0 + onTheRight.del*onTheRight.l0) + buffer) < centerDistance)
		return 0;


	// -------------------------
	// 
	// 	   Vertex forces
	//
	// -------------------------

	double distScale = 0.0;						// distance scale
	double l2 = a;								// interaction range
	double p2 = 1.0 + a;						// point where attractive zone ends (in units of delta)
	double vertexDist = 0.0;					// distance variable
	vector<double> vertexVec(NDIM,0.0);			// vector to hold vectorial distance quantity
	double contactDistance = 0.0;				// contact distance variable

	// loop over vertex pairs, check for contact
	for (i=0; i<NV; i++){
		for (j=0; j<onTheRight.NV; j++){

			// get distance between vertices i and j
			vertexDist = 0.0;
			for (d=0; d<NDIM; d++){
				// get distance to nearest image
				distTmp = distance(onTheRight.vpos(j,d),vpos(i,d),d);

				// add to vertex distance
				vertexVec.at(d) = distTmp;

				// add to scalar distance
				vertexDist += distTmp*distTmp;
			}

			// get contact distance
			contactDistance = 0.5*(del*l0 + onTheRight.del*onTheRight.l0);

			// get vertex distance
			vertexDist = sqrt(vertexDist);

			// check overlap distance, IF in attractive zone, increment number of attractive contacts
			if (vertexDist >= contactDistance && vertexDist < contactDistance*p2 && l2 > 0.0)
				numAttractiveContacts++;
		}
	}

	// return if in contact or not
	return numAttractiveContacts;
}


// harmonic spring force between centers of overlapping particles
int deformableParticles2D::radialForce(deformableParticles2D &onTheRight, double bscale){
	// local variables
	int i, j, d, maxI, maxJ;
	double viDotProduct, vDotProduct, ftmp, maxDotProduct, uTmp;
	double c1, c2;
	int inContact = 0;

	// calculate distance
	double centerDistance = 0.0; 			// center-to-center distance
	vector<double> deltaMuNu(NDIM,0.0);		// vector to store center-to-center distance vector
	double muREff,nuREff,buffer;			// variables to determine effect cell contact distance
	double distTmp = 0.0;					// temporary distance, for mimimum image convention (MIC)
	double contactDistance = 0.0;			// distance based on vertices closest to pair

	// buffer
	buffer = bscale*sqrt(area()/PI);

	// calculate connecting vector deltaMuNu
	for (d=0; d<NDIM; d++){
		// get distance
		distTmp = cellDistance(onTheRight,d);

		// save distance in vector
		deltaMuNu.at(d) = distTmp;

		// calculate vector norm
		centerDistance += pow(distTmp,2);
	}
	centerDistance = sqrt(centerDistance);

	// get contact distance by assuming regular polygons
	c1 = sqrt((2.0*a0)/(NV*sin(2*PI/NV)));
	c2 = sqrt((2.0*onTheRight.geta0())/(onTheRight.getNV()*sin(2*PI/onTheRight.getNV())));
	// c1 = sqrt(area()/PI);
	// c2 = sqrt(onTheRight.area()/PI);
	contactDistance = c1 + c2 + buffer;

	cout << "contactDistance = " << contactDistance << ", centerDistance = " << centerDistance << endl;

	// check distance between centers
	if (centerDistance < contactDistance){
		inContact = 1;
		for (d=0; d<NDIM; d++){
			ftmp = -(kint/contactDistance)*(1 - (centerDistance/contactDistance))*deltaMuNu.at(d)/centerDistance;

			// distribute force to vertices
			for (i=0; i<NV; i++)
				setVForce(i,d,vforce(i,d) + ftmp/NV);
			for (j=0; j<onTheRight.getNV(); j++)
				onTheRight.setVForce(j,d,onTheRight.vforce(j,d) - ftmp/onTheRight.getNV());
		}

		uTmp =  0.5 * kint * pow((1 - (centerDistance/contactDistance)),2);

		for (i=0; i<NV; i++)
			setUInt(i,uInt(i) + uTmp/NV);
		for (j=0; j<onTheRight.getNV(); j++)
			onTheRight.setUInt(j,onTheRight.uInt(j) + uTmp/onTheRight.getNV());
	}

	return inContact;
}

/************************

	Energies

*************************/

double deformableParticles2D::perimeterEnergy(){
	// local variables
	int i;
	double val = 0.0;
	double El;


	if (kl > 0){
		// effective energy scale
		El = 0.5*kl*NV*l0*l0;

		// sum squared perimeter strains
		for (i=0; i<NV; i++)
			val += pow((segmentLength(i)/l0) - 1.0,2.0);

		// multiply by energy scale
		val *= El;

		// return energy value
		return val;
	}
	else
		return 0.0;
	

	/*
	NOTE: ABOVE VERSION FOR BALANCED FORCE VERSION!

	THIS VERSION FOR ORGINAL MODEL, WITH DIMENSIONAL STRAINS
	if (kl > 0){
		// energy based on perimeter deviations
		for (i=0; i<NV; i++)
			val += 0.5*kl*pow(segmentLength(i) - l0,2.0);

		// return value
		return val;
	}
	else
		return 0.0;
	*/
}

double deformableParticles2D::areaEnergy(){
	// local variables
	double val = 0.0;
	double Ea, aStrain;

	if (ka > 0){
		// effective area energy scale 
		Ea = 0.5*ka*pow(a0/l0,2.0);

		// dimensionless strain
		aStrain = (polygonArea()/a0) - 1.0;

		// calculate and return energy
		val = Ea*aStrain*aStrain;
		return val;
	}
	else
		return 0.0;

	/*
	NOTE: ABOVE VERSION FOR BALANCED FORCE VERSION!

	THIS VERSION FOR ORGINAL MODEL, WITH DIMENSIONAL STRAINS

	if (ka > 0){
		// calculate energy
		val = 0.5*ka*pow(area()-a0,2.0);

		// return value
		return val;
	}
	else
		return 0.0;
	*/
}

double deformableParticles2D::surfaceTensionEnergy(){
	// local variables
	int i;
	double val = 0.0;

	if (gam > 0){
		// loop over vertices to calculate energy
		for (i=0; i<NV; i++)
			val += gam*segmentLength(i);

		// return value
		return val;
	}
	else
		return 0.0;
}

double deformableParticles2D::bendEnergy(){
	// local variables
	int i;
	double val = 0.0;
	double Eb, cStrainI, normC0;

	if (kb > 0){
		// bending energy scale
		Eb = 0.5*kb*NV*l0*l0;

		// cosine normalization
		normC0 = 1 + c0;

		// loop over vertices to calculate energy
		for (i=0; i<NV; i++){
			cStrainI = (segmentCosine(i) - c0)/normC0;
			val += cStrainI*cStrainI;
		}
		val *= Eb;

		// return value
		return val;
	}
	else
		return 0.0;

	/*
	NOTE: ABOVE VERSION FOR BALANCED FORCE VERSION!

	THIS VERSION FOR ORGINAL MODEL, WITH DIMENSIONAL STRAINS

	// check if energy is activated
	if (kb > 0){
		// loop over vertices to calculate energy
		for (i=0; i<NV; i++)
			val += 0.5*kb*pow(segmentCosine(i) - c0,2.0);

		// return value
		return val;
	}
	else
		return 0.0;
	*/
}

double deformableParticles2D::interactionEnergy(){
	// local variables
	int i;
	double val = 0.0;

	if (kint > 0){
		// loop over vertices to calculate energy
		for (i=0; i<NV; i++)
			val += uInt(i);

		// return value
		return val;
	}
	else
		return 0.0;
}

double deformableParticles2D::totalPotentialEnergy(){
	// local variables
	int i;
	double totalUInt = 0.0;

	// loop over vertices, add to totalUInt
	for (i=0; i<NV; i++)
		totalUInt += uInt(i);


	// return total potential energy
	return perimeterEnergy() + areaEnergy() + surfaceTensionEnergy() + bendEnergy() + totalUInt;
}

double deformableParticles2D::totalKineticEnergy(){
	// local variables
	int i,d;
	double val = 0.0;

	// vertex mass
	double segmentMass = PI*pow(0.5*del*l0,2.0);

	// loop over vertices, get kinetic energy
	for (i=0; i<NV; i++){
		for (d=0; d<NDIM; d++)
			val += 0.5*segmentMass*vvel(i,d)*vvel(i,d);
	}

	// return value
	return val;
}





/************************

	MD Integrators

*************************/


void deformableParticles2D::verletPositionUpdate(double dt){
	// local variables
	int i,d;
	double postmp;

	// update vertex positions
	for (i=0; i<NV; i++){
		for (d=0; d<NDIM; d++){
			// update positions using velocity-Verlet with PBCs
			postmp = vpos(i,d) + dt*vvel(i,d) + 0.5*vacc(i,d)*dt*dt;
			setVPos(i,d,postmp);

			// set forces to 0
			setVForce(i,d,0.0);
		}

		// set uint to 0
		setUInt(i,0.0);
	}
}


void deformableParticles2D::verletVelocityUpdate(double dt){
	// local variables
	int i,d;
	double veltmp,anew,segmentMass,b;
	double ftmp, dampNum, dampDenom, dampUpdate;

	// get segment mass
	segmentMass = PI*pow(0.5*del*l0,2);

	// update vertex velocities
	for (i=0; i<NV; i++){
		for (d=0; d<NDIM; d++){
			// get current velocities	
			veltmp = vvel(i,d);

			// get the new acceleration from forces with damping
			anew = vforce(i,d)/segmentMass;

			// update velocity
			veltmp += 0.5*dt*(anew + vacc(i,d));

			// set new velocity and acceleration
			setVVel(i,d,veltmp);
			setVAcc(i,d,anew);
		}
	}
}

void deformableParticles2D::verletVelocityUpdate(double dt, double dampingParam){
    // local variables                                                                                                                                                                                                                 
    int i,d;
    double veltmp,anew,segmentMass,b;
    double ftmp, dampNum, dampDenom, dampUpdate;

    // get segment mass                                                                                                                                                                                                                
    segmentMass = PI*pow(0.5*del*l0,2);

    // scale damping                                                                                                                                                                                                                   
    b = dampingParam*sqrt(segmentMass);

    // update vertex velocities                                                                                                                                                                                                        
    for (i=0; i<NV; i++){
		for (d=0; d<NDIM; d++){
            // get current velocities                                                                                                                                                                                          
            veltmp = vvel(i,d);

            // calculate force damping update                                                                                                                                                                                  
            ftmp = vforce(i,d);
            dampNum = b*(veltmp - 0.5*vacc(i,d)*dt);
            dampDenom = 1.0 + 0.5*b*dt;
            dampUpdate = (ftmp - dampNum)/dampDenom;

			// update force                                                                                                                                                                                                    
            setVForce(i,d,dampUpdate);

	        // get the new acceleration from forces with damping                                                                                                                                                               
            anew = vforce(i,d)/segmentMass;

			// update velocity                                                                                                                                                                                                 
			veltmp += 0.5*dt*(anew + vacc(i,d));

            // set new velocity and acceleration                                                                                                                                                                               
			setVVel(i,d,veltmp);
            setVAcc(i,d,anew);
		}
    }
}




/************************

	Print Functions

*************************/

void deformableParticles2D::printVertexPositions(ofstream& vertexPrintObject, int cellID, int frame){
	// print variables
	int wID = 6; 	// width of ID
	int wNAME = 12;	// width of descriptor

	// print number of vertices as header
	vertexPrintObject << setw(wNAME) << left << "NEWFR" << " " << endl;
	vertexPrintObject << setw(wNAME) << left << "FRAME" << setw(wID) << right << frame << endl;
	printVertexPositions(vertexPrintObject,cellID);
}

void deformableParticles2D::printVertexPositions(ofstream& vertexPrintObject, int cellID){
	// local variables
	int i,d;

	// print variables
	int wID = 6; 	// width of ID
	int wNAME = 12;	// width of descriptor
	int wNUM = 30;	// width of number
	int p = 16;		// variable precision

	// print number of vertices as header
	vertexPrintObject << setw(wNAME) << left << "NVERT" << setw(wID) << right << NV << endl;

	// print cell info
	vertexPrintObject << setw(wNAME) << left << "CELLP";
	vertexPrintObject << setw(wID) << right << cellID;
	for (d=0; d<NDIM; d++)
		vertexPrintObject << setw(wNUM) << setprecision(p) << right << cpos(d);
	vertexPrintObject << setw(wNUM) << setprecision(p) << l0;
	vertexPrintObject << setw(wNUM) << setprecision(p) << a0;
	vertexPrintObject << setw(wNUM) << setprecision(p) << del;
	vertexPrintObject << setw(wNUM) << setprecision(p) << polygonArea();
	vertexPrintObject << setw(wNUM) << setprecision(p) << area();
	vertexPrintObject << setw(wNUM) << setprecision(p) << perimeter();
	vertexPrintObject << endl;

	// print vertex column header
	vertexPrintObject << setw(wNAME) << left << "CINFO" << setw(wID) << right << "vID";
	vertexPrintObject << setw(wNUM) << right << "X pos";
	vertexPrintObject << setw(wNUM) << right << "Y pos";
	vertexPrintObject << setw(wNUM) << right << "X vel";
	vertexPrintObject << setw(wNUM) << right << "Y vel";
	vertexPrintObject << setw(wNUM) << right << "X frc";
	vertexPrintObject << setw(wNUM) << right << "Y frc";
	vertexPrintObject << endl;

	// loop over vertices, print
	for (i=0; i<NV; i++){
		vertexPrintObject << setw(wNAME) << left << "VERTP";
		vertexPrintObject << setw(wID) << right << i;
		for (d=0; d<NDIM; d++){
			if (vpos(i,d) - cpos(d) > 0.5*L.at(d))
				vertexPrintObject << setw(wNUM) << setprecision(p) << right << vpos(i,d) - L.at(d);
			else if (cpos(d) - vpos(i,d) > 0.5*L.at(d))
				vertexPrintObject << setw(wNUM) << setprecision(p) << right << vpos(i,d) + L.at(d);
			else
				vertexPrintObject << setw(wNUM) << setprecision(p) << right << vpos(i,d);
		}
		for (d=0; d<NDIM; d++)
			vertexPrintObject << setw(wNUM) << setprecision(p) << right << vvel(i,d);
		for (d=0; d<NDIM; d++)
			vertexPrintObject << setw(wNUM) << setprecision(p) << right << vforce(i,d);
		vertexPrintObject << endl;
	}
}

void deformableParticles2D::printCellEnergy(ofstream& energyPrintObject, int frame){
	// local variables
	double uPerimeter, uArea, uSurfaceTension, uBend, uInteraction, uTotal, KTotal;

	// get energies
	uPerimeter = perimeterEnergy();
	uArea = areaEnergy();
	uSurfaceTension = surfaceTensionEnergy();
	uBend = bendEnergy();
	uInteraction = interactionEnergy();
	uTotal = uPerimeter + uArea + uSurfaceTension + uBend + uInteraction;
	KTotal = totalKineticEnergy();

	// print
	energyPrintObject << setw(6) << right << frame;
	energyPrintObject << setw(30) << setprecision(16) << right << uPerimeter;
	energyPrintObject << setw(30) << setprecision(16) << right << uArea;
	energyPrintObject << setw(30) << setprecision(16) << right << uSurfaceTension;
	energyPrintObject << setw(30) << setprecision(16) << right << uBend;
	energyPrintObject << setw(30) << setprecision(16) << right << uTotal;
	energyPrintObject << setw(30) << setprecision(16) << right << KTotal;
	energyPrintObject << endl;
}




