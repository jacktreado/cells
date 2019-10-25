/*

	Methods file for cellPacking2D class

*/


// include file
#include "deformableParticles2D.h"
#include "cellPacking2D.h"

// namespace
using namespace std;

// constants
const double PI = 4*atan(1);


/************************

	Constructors &
	Destructor

*************************/

// set up variables
void cellPacking2D::defaultvars(){
	// dimensionality ALWAYS set to 2
	NDIM 	= 2;

	// scalar variables set to 0
	NCELLS 		= 0;
	NT 			= 0;
	NPRINT 		= 0;
	dt 			= 0.0;
	dt0 		= 0.0;
	L 			= 0.0;
	T 			= -1.0;
	phi 		= -1.0;
	sigmaXX 	= 0.0;
	sigmaXY 	= 0.0;
	sigmaYX 	= 0.0;
	sigmaYY 	= 0.0;

	// pointer variables point to nullptr
	cellArray 				= nullptr;
	contactMatrix	 		= nullptr;

	// seed random numbers randomly
	srand48(seed);
}



// default constructor
cellPacking2D::cellPacking2D(){
	defaultvars();
}



// overloaded constructor with system information as arguments
cellPacking2D::cellPacking2D(int ncells, int nt, int nprint, double l, double s){
	// local variables
	int i,NC;

	// set initial seed
	seed = s;

	// first use starting variabls
	defaultvars();

	// set member variables based on inputs
	NCELLS 		= ncells;
	NT 			= nt;
	NPRINT 		= nprint;
	L 			= l;

	// test for error in inputs
	if (NCELLS <= 0){
		cout << "	ERROR: in overloaded operator, NCELLS <= 0 so cannot initialize arrays, ending code here." << endl;
		exit(1);
	}
	else if (NT <= 0){
		cout << "	ERROR: in overloaded operator, NT <= 0 so cannot run simulations, ending code here." << endl;
		exit(1);
	}
	else if (NPRINT <= 0){
		cout << "	ERROR: in overloaded operator, NPRINT <= 0 so cannot run simulations, ending code here." << endl;
		exit(1);
	}
	else if (L <= 0.0){
		cout << "	ERROR: in overloaded operator, L <= 0.0 so cannot run simulations, ending code here." << endl;
		exit(1);
	}

	// test that memory has not yet been initialized
	if (contactMatrix){
		cout << "	ERROR: in overloaded operator, contactMatrix ptr already initialized, ending code here." << endl;
		exit(1);
	}

	// initialize cell array
	cellArray = new deformableParticles2D[NCELLS];

	// initialize contact matrix
	NC = NCELLS*(NCELLS-1)/2;
	contactMatrix = new int[NC];

	// intialize contact matrix to 0
	for (i=0; i<NC; i++)
		contactMatrix[i] = 0;
}

// overloaded constructor for tumorous tissue simulation
cellPacking2D::cellPacking2D(int ncells, int ntumor, int tumorNV, int adiposeNV, double tumorCalA, double adiposeCalA, int s){
	// local variables
	int i, ci, NC;

	// set initial seed
	seed = s;

	// first use starting variabls
	defaultvars();

	// set member variables based on inputs
	NCELLS 		= ncells;
	L 			= 10.0*NCELLS;

	// test for error in inputs
	if (NCELLS <= 0){
		cout << "	ERROR: in adipose tissue constructor, NCELLS <= 0 so cannot initialize arrays, ending code here." << endl;
		exit(1);
	}
	else if (ntumor <= 0){
		cout << "	ERROR: in adipose tissue constructor, ntumor <= 0 so cannot initialize arrays, ending code here." << endl;
		exit(1);
	}
	else if (tumorNV <= 2){
		cout << "	ERROR: in adipose tissue constructor, tumorNV = " << tumorNV << " which is too small for cells, ending." << endl;
	}
	else if (adiposeNV <= 2){
		cout << "	ERROR: in adipose tissue constructor, adiposeNV = " << adiposeNV << " which is too small for cells, ending." << endl;
	}

	// test that memory has not yet been initialized
	if (contactMatrix){
		cout << "	ERROR: in overloaded operator, contactMatrix ptr already initialized, ending code here." << endl;
		exit(1);
	}

	// initialize cell array
	cellArray = new deformableParticles2D[NCELLS];

	// initialize contact matrix
	NC = NCELLS*(NCELLS-1)/2;
	contactMatrix = new int[NC];

	// intialize contact matrix to 0
	for (i=0; i<NC; i++)
		contactMatrix[i] = 0;

	// initialize cells as ntumor tumor cells, rest adipose
	for (ci=0; ci<NCELLS; ci++){

		// set box length for each cell
		cell(ci).setL(L);

		// tumor cells
		if (ci < ntumor){
			// initialize arrays in tumor cells
			cell(ci).setNV(tumorNV);
			cell(ci).initializeVertices();
			cell(ci).initializeCell();

			// set tumor asphericity
			cell(ci).setAsphericity(tumorCalA);
		}

		// adipose cells
		else{
			// initialize arrays in adipose cells
			cell(ci).setNV(adiposeNV);
			cell(ci).initializeVertices();
			cell(ci).initializeCell();

			// set adipose asphericity
			cell(ci).setAsphericity(adiposeCalA);
		}
	}

	// force constants for adipose simulation
	double kl,ka,gam,kb,kint,del,a;

	// set force constant values
	kl 		= 2.0;
	ka 		= 2.0;
	gam 	= 0.0;
	kb 		= 0.0;
	kint 	= 1.0;
	del 	= 1.0;		// USING VERTEX FORCES, SO DEL MUST BE SET TO L0!
	a 		= 0.0;

	// set values
	initializeForceConstants(kl,ka,gam,kb,kint);
	initializeInteractionParams(del,a);

	// turn off bending energy for tumor cells
	for (ci=0; ci<ntumor; ci++)
		cell(ci).setkb(0.0);

	// initialize cell positions on a square lattice
	squareLattice();
}

// overloaded constructor with file stream object as argument
cellPacking2D::cellPacking2D(ifstream& inputFileObject, double asphericity, double s){
	// set initial seed
	seed = s;

	// set variables to default
	defaultvars();

	// local variables
	int NC, ci, vi, nv, d;
	double l0tmp, l0, a0, posTmp;

	// read in simulation information
	inputFileObject >> NCELLS;
	inputFileObject >> L;
	inputFileObject >> l0tmp;

	cout << "NCELLS = " << NCELLS << endl;
	cout << "L = " << L << endl;

	// reset length scales to be in units of l0
	L /= l0tmp;

	// test for error in inputs
	if (NCELLS <= 0){
		cout << "	ERROR: in overloaded operator, NCELLS <= 0 so cannot initialize arrays, ending code here." << endl;
		exit(1);
	}
	else if (L <= 0.0){
		cout << "	ERROR: in overloaded operator, L <= 0.0 so cannot run simulations, ending code here." << endl;
		exit(1);
	}

	// test that memory has not yet been initialized
	if (contactMatrix){
		cout << "	ERROR: in overloaded operator, contactMatrix ptr already initialized, ending code here." << endl;
		exit(1);
	}

	// initialize cell array
	cellArray = new deformableParticles2D[NCELLS];

	// loop over cells
	for (ci=0; ci<NCELLS; ci++){
		// number of vertices
		inputFileObject >> nv;

		// initialize cell objects
		cell(ci).setL(L);
		cell(ci).setNV(nv);
		cell(ci).initializeVertices();
		cell(ci).initializeCell();

		// set cell com position
		for (d=0; d<NDIM; d++){
			inputFileObject >> posTmp;
			cell(ci).setCPos(d,posTmp/l0tmp);
		}

		// set vertex positions
		for (vi=0; vi<nv; vi++){
			for (d=0; d<NDIM; d++){
				inputFileObject >> posTmp;
				cell(ci).setVPos(vi,d,posTmp/l0tmp);
			}
		}

		// determine new l0 for given area
		l0 = (1.0/nv)*sqrt(4*PI*cell(ci).area()*asphericity);
		a0 = (nv*nv*l0*l0)/(4*PI*asphericity);

		// set a0 to enforce asphericity in cell i
		cell(ci).seta0(a0);
		cell(ci).setl0(l0);
	}

	// initialize contact matrix
	NC = NCELLS*(NCELLS-1)/2;
	contactMatrix = new int[NC];

	// intialize contact matrix to 0
	for (ci=0; ci<NC; ci++)
		contactMatrix[ci] = 0;
}


// destructor
cellPacking2D::~cellPacking2D(){
	if (cellArray){
		delete [] cellArray;
		cellArray = nullptr;
	}
	if (contactMatrix){
		delete [] contactMatrix;
		contactMatrix = nullptr;
	}

	// close file objects
	packingPrintObject.close();
	energyPrintObject.close();
	statPrintObject.close();
}



/************************

	Operators

*************************/

// assignment operator (ASSUME EVERYTHING HAS BEEN INITIALIZED)
void cellPacking2D::operator=(cellPacking2D& onTheRight){
	// local variables
	int ci,cj;

	// update packing fraction of onTheRight
	onTheRight.phi = onTheRight.packingFraction();

	// copy scalar variables
	NCELLS 	= onTheRight.NCELLS;
	NT 		= onTheRight.NT;
	NPRINT 	= onTheRight.NPRINT;

	seed 	= onTheRight.seed;
	dt 		= onTheRight.dt;
	dt0 	= onTheRight.dt0;
	L 		= onTheRight.L;
	T 		= onTheRight.T;
	phi 	= onTheRight.phi;

	sigmaXX = onTheRight.sigmaXX;
	sigmaXY = onTheRight.sigmaXY;
	sigmaYX = onTheRight.sigmaYX;
	sigmaYY = onTheRight.sigmaYY;

	// test that memory has not yet been initialized
	if (cellArray){
		cout << "	ERROR: in overloaded operator, cellArray ptr already initialized, ending code here." << endl;
		exit(1);
	}
	if (contactMatrix){
		cout << "	ERROR: in overloaded operator, contactMatrix ptr already initialized, ending code here." << endl;
		exit(1);
	}

	// initialize cell array
	cellArray = new deformableParticles2D[NCELLS];
	contactMatrix = new int[NCELLS*(NCELLS-1)/2];

	// deep copy cell objects and contact matrix
	for (ci=0; ci<NCELLS; ci++){
		// copy cell objects (using overloaded operator in deformableParticle2D class)
		cell(ci) = onTheRight.cell(ci);

		// copy elements from contact matrix
		for (cj=ci+1; cj<NCELLS; cj++){
			if (onTheRight.contacts(ci,cj))
				addContact(ci,cj);
		}
	}
}





/************************

	Initialization

*************************/


// monodisperse system
void cellPacking2D::initializeMonodisperse(int NV, double asphericity){
	// local variables
	int i;
	double a0;

	// check that NCELLS has been initialized
	if (NCELLS <= 0){
		cout << "	ERROR: in initializeMonodisperse(), NCELLS <= 0 so cannot initialize arrays, ending code here." << endl;
		exit(1);
	}
	if (NV <= 2){
		cout << "	ERROR: in initializeMonodisperse(), NV = " << NV << " which is too small for cells, ending." << endl;
		exit(1);
	}

	// check that cell array memory has been allocated
	if (!cellArray){
		cout << "	ERROR: in initializeMonodisperse(),cellArray memory has not been allocated, ending." << endl;
		exit(1);
	}

	// setup each cell
	for (i=0; i<NCELLS; i++){
		// set box length for each cell
		cell(i).setL(L);

		// need to initialize stuff in each cell object
		cell(i).setNV(NV);
		cell(i).initializeVertices();
		cell(i).initializeCell();

		// calculate a0 based on NV and asphericity
		a0 = (NV*NV)/(4.0*PI*asphericity);

		// set a0 to enforce asphericity in cell i
		cell(i).seta0(a0);
	}
}

// bidisperse system initialized to regular polygons
void cellPacking2D::initializeBidisperse(int NV, double sizeRatio){
	// local variables
	int i,nv,nvLarge,smallIndex;
	double a0,asphericity;

	// check that NCELLS has been initialized
	if (NCELLS <= 0){
		cout << "	ERROR: in initializeMonodisperse(), NCELLS <= 0 so cannot initialize arrays, ending code here." << endl;
		exit(1);
	}
	if (NV <= 2){
		cout << "	ERROR: in initializeMonodisperse(), NV = " << NV << " which is too small for cells, ending." << endl;
	}

	// check that cell array memory has been allocated
	if (!cellArray){
		cout << "	ERROR: in initializeBidisperse(),cellArray memory has not been allocated, ending." << endl;
		exit(1);
	}

	// calculate smaller nv to enforce bidispersity
	nvLarge = round(sizeRatio*NV);
	smallIndex = NCELLS/2;

	// setup each cell
	for (i=0; i<NCELLS; i++){
		// set box length for each cell
		cell(i).setL(L);

		// determine number of vertices
		if (i >= smallIndex)
			nv = nvLarge;
		else
			nv = NV;

		// need to initialize stuff in each cell object
		cell(i).setNV(nv);
		cell(i).initializeVertices();
		cell(i).initializeCell();

		// calculate a0 based on fact that they are regular polygons
		asphericity = nv*tan(PI/nv)/PI;
		cell(i).setAsphericity(asphericity);
	}
}


// bidisperse system
void cellPacking2D::initializeBidisperse(int NV, double sizeRatio, double asphericity){
	// local variables
	int i,nv,nvLarge,smallIndex;
	double a0;

	// check that NCELLS has been initialized
	if (NCELLS <= 0){
		cout << "	ERROR: in initializeMonodisperse(), NCELLS <= 0 so cannot initialize arrays, ending code here." << endl;
		exit(1);
	}
	if (NV <= 2){
		cout << "	ERROR: in initializeMonodisperse(), NV = " << NV << " which is too small for cells, ending." << endl;
	}

	// check that cell array memory has been allocated
	if (!cellArray){
		cout << "	ERROR: in initializeBidisperse(),cellArray memory has not been allocated, ending." << endl;
		exit(1);
	}

	// calculate smaller nv to enforce bidispersity
	nvLarge = round(sizeRatio*NV);
	smallIndex = NCELLS/2;

	// setup each cell
	for (i=0; i<NCELLS; i++){
		// set box length for each cell
		cell(i).setL(L);

		// determine number of vertices
		if (i >= smallIndex)
			nv = nvLarge;
		else
			nv = NV;

		// need to initialize stuff in each cell object
		cell(i).setNV(nv);
		cell(i).initializeVertices();
		cell(i).initializeCell();
		cell(i).setAsphericity(asphericity);
	}
}


// set force constants for each cell in system
void cellPacking2D::initializeForceConstants(double kl, double ka, double gam, double kb, double kint){
	// local variables
	int i;

	// check that NCELLS has been initialized
	if (NCELLS <= 0){
		cout << "	ERROR: in initializeForceConstants(), NCELLS <= 0 so cannot initialize arrays, ending code here." << endl;
		exit(1);
	}

	// loop over cells, set inputs
	for (i=0; i<NCELLS; i++){
		cell(i).setkl(kl);
		cell(i).setka(ka);
		cell(i).setgam(gam);
		cell(i).setkb(kb);
		cell(i).setkint(kint);
	}
}


// set interaction parameters
void cellPacking2D::initializeInteractionParams(double del, double a){
	// local variables
	int i;

	// check that NCELLS has been initialized
	if (NCELLS <= 0){
		cout << "	ERROR: in initializeInteractionParams(), NCELLS <= 0 so cannot initialize arrays, ending code here." << endl;
		exit(1);
	}

	// loop over cells
	for (i=0; i<NCELLS; i++){
		cell(i).setdel(del*cell(i).getl0());
		cell(i).seta(a);
	}
}


// set initial positions to populate a square lattice
void cellPacking2D::squareLattice(){
	// local variables
	int ci, xIndex, yIndex;
	int gridPoints = round(1.5*NCELLS);

	double xpos, ypos;
	double buffer = 0.05*L;
	double spacing = (L - 2.0*buffer)/(gridPoints-1);

	// loop over cells, give random initial cell positions, initialize as regular polygons
	for (ci=0; ci<NCELLS; ci++){
		// get random lattice indices
		xIndex = ceil(drand48()*(gridPoints-1));
		yIndex = ceil(drand48()*(gridPoints-1));

		// get locations
		xpos = buffer + xIndex*spacing;
		ypos = buffer + yIndex*spacing;

		// set positions
		cell(ci).setCPos(0,xpos);
		cell(ci).setCPos(1,ypos);

		// initialize as regular polygon
		cell(ci).regularPolygon();

		// perturb vertices a lil bit
		cell(ci).vertexPerturbation(0.1);

		// output cell asphericities
		cout << "initializing cell " << ci << " on square lattice, initial asphericity = " << cell(ci).asphericity();
		cout << " and calA0 = " << ((double)cell(ci).getNV()*cell(ci).getNV()*cell(ci).getl0()*cell(ci).getl0())/(4.0*PI*cell(ci).geta0()) << endl;
	}

	// calculate packing fraction
	phi = packingFraction();

	// print statement
	cout << "Particles initialized on square lattice with initial packing fraction phi = " << phi << endl;
}


// set positions as if particles were bidisperse spheres
void cellPacking2D::bidisperseDisks(double sizeRatio, double phiT){
	// local variables
	int it = 0;
	int itMax = 1e5;
	int ci;
	double xpos, ypos, T0;
	double phiTol = 1e-8;
	double rSmall, rLarge;
	double phiInit = 0.1;
	double phi = phiInit;
	double dphi = 0.01;
	double dphiTmp = dphi;

	// set random initial positions, cells as regular polygon
	for (ci=0; ci<NCELLS; ci++){
		// get random location in box
		xpos = L*drand48();
		ypos = L*drand48();

		// set as initial position of com
		cell(ci).setCPos(0,xpos);
		cell(ci).setCPos(1,ypos);

		// set particles as regular polygons
		cell(ci).regularPolygon();
	}

	// initialize velocities
	T0 = 1e-8;
	cout << "\n\nIn initial bidisperse disk packing, setting initial T = " << T0 << endl;
	initializeVelocities(T0);

	// loop over packing fractions
	while (phi < phiT){

		// relax overlaps
		cout << "Relieving overlaps of particles as disks" << endl;
		overlapRelief(phi);

		// update phi appropriately
		if (phi + dphi < phiT)
			phi += dphi;
		else{
			dphiTmp = dphi;
			while(phi + dphiTmp > phiT && dphiTmp > 1e-8)
				dphiTmp *= 0.9;
			phi += dphiTmp;
		}
		cout << "updating phi to phi + dphi = " << phi << endl;
	}
} 


// set initial COM velocities based on temperature
void cellPacking2D::initializeVelocities(double tmp0){
	// local variables
	int ci,d;
	double rv, vscale, ek;
	vector<double> vmean(NDIM,0.0);

	// set initial temperature
	T = tmp0;

	// initialize random velocities
	for (ci=0; ci<NCELLS; ci++){

		// add to velocities and mean
		for (d=0; d<NDIM; d++){
			// get random direction
			rv = drand48();

			cell(ci).setCVel(d,rv);
			vmean.at(d) += rv;
		}

	}

	// get mean
	for (d=0; d<NDIM; d++)
		vmean.at(d) /= NCELLS;

	// subtract of mean, calc EK
	ek = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		for (d=0; d<NDIM; d++){
			// subtract of com motion
			cell(ci).setCVel(d,cell(ci).cvel(d) - vmean.at(d));

			// calc ek
			ek += 0.5*pow(cell(ci).cvel(d),2);
		}
	}

	// get vscale
	vscale = sqrt(T/(cell(0).area()*ek));
	for (ci=0; ci<NCELLS; ci++){
    	for (d=0; d<NDIM; d++)
        	cell(ci).setCVel(d,cell(ci).cvel(d)*vscale);
    }
}


void cellPacking2D::initializeVelocities(int ci, double tmp0){
	// local variables
	int d;
	double rv, vscale, ek;

	// create random velocities
	for (d=0; d<NDIM; d++){
		rv = drand48();
		cell(ci).setCVel(d,rv);
	}

	// calc kinetic energy
	ek = 0.0;
	for (d=0; d<NDIM; d++)
		ek += 0.5*pow(cell(ci).cvel(d),2);

	// get vscale
	vscale = sqrt(tmp0/(cell(ci).area()*ek));
	for (d=0; d<NDIM; d++)
    	cell(ci).setCVel(d,cell(ci).cvel(d)*vscale);
}







/************************

	Getters

*************************/

// return object of cell ci
deformableParticles2D& cellPacking2D::cell(int ci){
	// check input is OK
	if (ci >= NCELLS || ci < 0){
		cout << "	ERROR: in cell(), ci = " << ci << " which is out of bounds for cellArray, ending code here." << endl;
		exit(1);
	}

	// return object from array of cells
	return cellArray[ci];
}


// number of frames in simulation, given NT and NPRINT
int cellPacking2D::nframes(){
	return ceil(NT/NPRINT);
}



// index of contact matrix
int cellPacking2D::cmindex(int ci, int cj){
	// check input is OK
	if (ci >= NCELLS || ci < 0){
		cout << "	ERROR: in cmindex(), ci = " << ci << " which is out of bounds for contactMatrix, ending code here." << endl;
		exit(1);
	}
	if (cj >= NCELLS || cj < 0){
		cout << "	ERROR: in cmindex(), cj = " << cj << " which is out of bounds for contactMatrix, ending code here." << endl;
		exit(1);
	}

	if (ci > cj)
		return NCELLS*cj + ci - (cj+1)*(cj+2)/2;
	else
		return NCELLS*ci + cj - (ci+1)*(ci+2)/2;
}


// entry in contact matrix
int cellPacking2D::contacts(int ci, int cj){
	// local variables
	int index = cmindex(ci,cj);
	int totalEntries = NCELLS*(NCELLS-1)/2;

	// return entry in contactMatrix array
	return contactMatrix[cmindex(ci,cj)];
}





/************************

	Calculations

*************************/


// total number of particle-particle contacts
int cellPacking2D::totalNumberOfContacts(){
	// local variables
	int ci,cj;
	int val = 0;

	// loop over particle contacts
	for (ci=0; ci<NCELLS; ci++){
		for (cj=ci+1; cj<NCELLS; cj++)
			val += contacts(ci,cj);
	}

	// return number of contacts
	return val;
}


// characteristic time scale in system, based on vertex spring fluctuations
double cellPacking2D::timeScale(){
	// local variables
	double charMass, charLength, charEnergy, val;

	// characteristic mass
	// charMass = cell(0).getdel()*(cell(0).getl0() + PI*0.5*cell(0).getdel());
	charMass = 0.25*PI*cell(0).getdel()*cell(0).getdel(); 		// needs to be this for mass if vertex forces only
	charLength = cell(0).getdel();
	charEnergy = cell(0).getkint();

	// time scale value
	val = sqrt(charMass*charLength*charLength/charEnergy);

	return val;
}


// Calculate the packing fraction
double cellPacking2D::packingFraction(){
	// local variables
	int ci, vi;
	double vrad, segangle;
	double val = 0.0;

	// loop over cells, packing fraction is : triangular area + 0.5*delta*perimeter area + area of circular corners
	for (ci=0; ci<NCELLS; ci++){
		// val += cell(ci).area() + (0.5*cell(ci).getdel())*cell(ci).perimeter() + PI*0.25*cell(ci).getdel()*cell(ci).getdel();
		val += cell(ci).area();		// needs to be this if vertex forces only

		// add contribution from vertices
		for (vi=0; vi<cell(ci).getNV(); vi++){
			// vertex radius
			vrad = 0.5*cell(ci).getdel()*cell(ci).getl0();

			// angle between successive segments
			segangle = acos(cell(ci).segmentCosine(vi));

			// add extra area
			val += 0.5*pow(vrad,2)*(PI+segangle);
		}
	}

	// divide by box area
	val /= pow(L,NDIM);

	// return value
	return val;
}


// Calculate shape potential energy in the system
double cellPacking2D::shapePotentialEnergy(){
	// local variables
	int ci;
	double val = 0.0;

	// loop over cells, add potential energy
	for (ci=0; ci<NCELLS; ci++)
		val += (cell(ci).perimeterEnergy() + cell(ci).areaEnergy());

	// return value
	return val;
}

double cellPacking2D::relaxPotentialEnergy(){
	// local variables
	int ci;
	double val = 0.0;

	// loop over cells, add potential energy
	for (ci=0; ci<NCELLS; ci++)
		val += (cell(ci).perimeterEnergy() + cell(ci).areaEnergy() + cell(ci).interactionEnergy());

	// return value
	return val;
}

// Calculate total potential energy in system
double cellPacking2D::totalPotentialEnergy(){
	// local variables
	int ci;
	double val = 0.0;

	// loop over cells, add potential energy
	for (ci=0; ci<NCELLS; ci++)
		val += cell(ci).totalPotentialEnergy();

	// return value
	return val;
}


// Calculate interaction potential between cells
double cellPacking2D::interactionPotentialEnergy(){
	// local variables	
	int ci,vi;
	double val = 0.0;

	// loop over cells, sum over interaction potential only
	for (ci=0; ci<NCELLS; ci++){
		for (vi=0; vi<cell(ci).getNV(); vi++)
			val += cell(ci).uInt(vi);
	}

	// return value
	return val;
}


// Calculate total kinetic energy in system
double cellPacking2D::totalKineticEnergy(){
	// local variables
	int ci;
	double val = 0.0;

	// loop over cells, add kinetic energy
	for (ci=0; ci<NCELLS; ci++)
		val += cell(ci).totalKineticEnergy();

	// return value
	return val;
}


// Calculate magnitude of largest force
double cellPacking2D::maxForceMagnitude(){
	// local variables
	int ci,vi,d,nv;
	double ftmp, maxForce;

	// set initial maxForce
	maxForce = 0.0;

	// loop over all forces on all vertices
	for (ci=0; ci<NCELLS; ci++){
		// number of vertices on cell ci
		nv = cell(ci).getNV();

		// loop over vertices
		for (vi=0; vi<nv; vi++){

			// check force
			for (d=0; d<NDIM; d++){
				ftmp = cell(ci).vforce(vi,d)*cell(ci).vforce(vi,d);
				if (ftmp > maxForce)
					maxForce = ftmp;
			}
		}
	}

	// return max force
	return sqrt(maxForce);
}

double cellPacking2D::maxNetForceMagnitude(){
	// local variables
	int ci,vi,d,nv;
	double maxForce, testForce;
	vector<double> ftmp(NDIM,0.0);

	// set initial maxForce
	maxForce = 0.0;

	// loop over all forces on all vertices
	for (ci=0; ci<NCELLS; ci++){
		// number of vertices on cell ci
		nv = cell(ci).getNV();

		// reset force
		for (d=0; d<NDIM; d++)
			ftmp.at(d) = 0.0;

		// loop over vertices
		for (vi=0; vi<nv; vi++){

			// calc force
			for (d=0; d<NDIM; d++)
				ftmp.at(d) += cell(ci).vforce(vi,d);
		}

		// get initial test force
		testForce = 0.0;
		for (d=0; d<NDIM; d++)
			testForce += ftmp.at(d)*ftmp.at(d); 

		// check again max
		if (testForce > maxForce)
			maxForce = testForce;
	}

	// return max force
	return sqrt(maxForce);
}

double cellPacking2D::forceRMS(){
	int ci, vi, d;
	int NVTOTAL = 0;
	double frms = 0.0;

	// loop over forces, calc total force norm
	for (ci=0; ci<NCELLS; ci++){
		NVTOTAL += cell(ci).getNV();
		for (vi=0; vi<cell(ci).getNV(); vi++){
			for (d=0; d<NDIM; d++)
				frms += pow(cell(ci).vforce(vi,d),2);
		}
	}

	// get force scale
	frms = sqrt(frms/(NDIM*NVTOTAL));

	// return
	return frms;
}

// calculate mean asphericity of all particles
double cellPacking2D::meanAsphericity(){
	// local variables
	int ci;
	double val = 0.0;

	// loop over cells, calc mean asphericity
	for (ci=0; ci<NCELLS; ci++)
		val += cell(ci).asphericity();
	val /= NCELLS;

	return val;
}


/************************

	Setters

*************************/


// set value in contact matrix to 1
void cellPacking2D::addContact(int ci, int cj){
	contactMatrix[cmindex(ci,cj)] = 1;
}


// set value in contact matrix to 0
void cellPacking2D::deleteContact(int ci, int cj){
	contactMatrix[cmindex(ci,cj)] = 0;
}


// delete all contacts to reset matrix to 0
void cellPacking2D::resetContacts(){
	for (int ci=0; ci<NCELLS; ci++){
		for (int cj=ci+1; cj<NCELLS; cj++)
			deleteContact(ci,cj);
	}
}

int cellPacking2D::particleContacts(int ci){
	// local variables
	int cj, nc = 0;

	// loop over other cells, count contacts
	for (cj=0; cj<NCELLS; cj++)
		if (cj != ci)
			nc += contacts(ci,cj);

	// return
	return nc;
}


// set packing fraction to desired value
void cellPacking2D::setPackingFraction(double val){
	// local variables
	int i;
	double scaleFactor;

	// update packing fraction
	// phi = packingFraction();

	// calculate val to scale lengths with
	scaleFactor	= pow(val/phi,1.0/NDIM);

	// scale all lengths by scale factor
	scaleLengths(scaleFactor);

	// update new phi
	phi = packingFraction();
}


// scale all length scales in the system by scaleFactor
void cellPacking2D::scaleLengths(double scaleFactor){
	// local variables
	int i;

	// // scale time based on dim analysis of MD time scale
	// dt *= pow(scaleFactor,0.5*NDIM);
	// dt0 *= pow(scaleFactor,0.5*NDIM);

	// loop over cells, use scale to change lengths
	for (i=0; i<NCELLS; i++)
		cell(i).scale(scaleFactor);
}

void cellPacking2D::setAsphericity(double val){
	// local variables
	int ci;

	// set all cells to specified asphericity values
	for (ci=0; ci<NCELLS; ci++)
		cell(ci).setAsphericityConstA(val);
}

// change asphericity on cell ci
void cellPacking2D::setAsphericity(int ci, double val){
	cell(ci).setAsphericityConstA(val);
}


// rescale velocities according to set temperature
void cellPacking2D::rescaleVelocities(double temperature){
	// local variables
	int ci, vi, d;
	double currentTemp, vscale;

	// get current temperature
	currentTemp = totalKineticEnergy();

	// get vscale
	vscale = sqrt(temperature/currentTemp);

	// rescale
	for (ci=0; ci<NCELLS; ci++){
		for (vi=0; vi<cell(ci).getNV(); vi++){
			for (d=0; d<NDIM; d++)
				cell(ci).setVVel(vi,d,vscale*cell(ci).vvel(vi,d));
		}
	}
}


// remove rattlers
int cellPacking2D::removeRattlers(int krcrs){
	int ci, cj, r, nr, nm;

	// monitor recursion depth
	krcrs++;

	// number of rattlers
	nr = 0;

	// number of "marginal" rattlers to be removed
	nm = 0;

	// loop over rows, eliminate contacts to rattlers
	for (ci=0; ci<NCELLS; ci++) {
		// get number of contacts
		r = particleContacts(ci);

		// remove from network if r <= DOF, delete contacts
		if (r <= NDIM) {
			// increment # of rattlers
			nr++;

			// if in contact, remove contacts
			if (r > 0) {
				nm++;

				for (cj=0; cj<NCELLS; cj++) {
					// delete contact between ci and cj
					if (ci != cj)
						deleteContact(ci,cj);
				}
			}
		}
	}

	if (krcrs > 100) {
		cout << "max recursive depth reached, be wary of output" << endl;
		return -1;
	}

	// recursively check once rattlers are removed
	if (nm == 0)
		return nr;
	else
		return removeRattlers(krcrs);
}



/**************************

	Forces and position 
		updates

***************************/


// calculate all forces, both shape and pairwise
void cellPacking2D::calculateForces(){
	// local variables
	int ci,cj,d,dd,inContact;

	// vector to store pairwise forces
	// vector<double> virialStress(NDIM*NDIM,0.0);
	vector<double> fij(NDIM,0.0);
	vector<double> rij(NDIM,0.0);

	// reset virial stresses to 0
	sigmaXX = 0.0;
	sigmaXY = 0.0;
	sigmaYX = 0.0;
	sigmaYY = 0.0;

	// reset forces
	for (ci=0; ci<NCELLS; ci++){
		for (d=0; d<NDIM; d++)
			cell(ci).setCForce(d,0.0);
	}

	// loop over cells and cell pairs, calculate shape and interaction forces
	for (ci=0; ci<NCELLS; ci++){
		// loop over pairs, add info to contact matrix
		for (cj=ci+1; cj<NCELLS; cj++){
			// reset virial stresses
			for (d=0; d<NDIM; d++){
				fij.at(d) = 0.0;
				rij.at(d) = 0.0;
			}

			// calculate forces
			inContact = cell(ci).vertexForce(cell(cj),fij,rij);
			if (inContact == 1)
				addContact(ci,cj);

			// DEBUG: OUTPUT PAIRWISE FORCES
			// cout << "ci = " << ci << ", cj = " << cj << ": " << fij.at(0) << " " << fij.at(1) << endl;

			// compute virial stresses
			sigmaXX += fij.at(0)*rij.at(0);
			sigmaXY += fij.at(1)*rij.at(0);
			sigmaYX += fij.at(0)*rij.at(1);
			sigmaYY += fij.at(1)*rij.at(1);
		}

		// forces on vertices due to shape
		cell(ci).shapeForces();
	}

	// normalize virial stresses by the area
	sigmaXX /= pow(L,2);
	sigmaXY /= pow(L,2);
	sigmaYX /= pow(L,2);
	sigmaYY /= pow(L,2);
}

// fire energy minimization
void cellPacking2D::fireStep(int& np, double& alpha){
	// HARD CODE IN FIRE PARAMETERS
	const double alpha0 	= 0.25;
	const double finc 		= 1.1;
	const double fdec 		= 0.5;
	const double falpha 	= 0.99;
	const double dtmax 		= 10*dt0;
	const int NMIN 			= 5;

	// local variables
	int ci,vi,d;
	double P,vstarnrm,fstarnrm,vtmp,ftmp;

	// calculate P and norms
	P = 0.0;
	vstarnrm = 0.0;
	fstarnrm = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		for (vi=0; vi<cell(ci).getNV(); vi++){
			for (d=0; d<NDIM; d++){
				// get tmp variables
				ftmp = cell(ci).vforce(vi,d);
				vtmp = cell(ci).vvel(vi,d);

				// calculate based on all vertices on all cells
				P += ftmp*vtmp;
				vstarnrm += vtmp*vtmp;
				fstarnrm += ftmp*ftmp;
			}
		}
	}

	// get norms
	vstarnrm = sqrt(vstarnrm);
	fstarnrm = sqrt(fstarnrm);

	// update velocities if forces are acting
	if (fstarnrm > 0){
		for (ci=0; ci<NCELLS; ci++){
			for (vi=0; vi<cell(ci).getNV(); vi++){
				for (d=0; d<NDIM; d++){
					vtmp = (1 - alpha)*cell(ci).vvel(vi,d) + alpha*(cell(ci).vforce(vi,d)/fstarnrm)*vstarnrm;
					cell(ci).setVVel(vi,d,vtmp);
				}
			}
		}
	}

	// save current kinetic energy to temperature variable
	T = totalKineticEnergy();

	// update P and alpha
	if (P > 0 && np > NMIN){

		// increase dt
		if (dt * finc < dtmax)
			dt *= finc;
		else
			dt = dtmax;

		// decrease alpha
		alpha *= falpha;

		np++;
	}
	else if (P <= 0) {
		// decrease time step, but only to limit
		dt *= fdec;

		// set global velocity vector to zero
		for (ci=0; ci<NCELLS; ci++){
			for (vi=0; vi<cell(ci).getNV(); vi++){
				for (d=0; d<NDIM; d++){
					cell(ci).setVVel(vi,d,0.0);
				}
			}
		}

		// set alpha -> alpha0
		alpha = alpha0;

		// set np -> 0
		np = 0;
	}
	else if (P > 0 && np <= NMIN) 
		np++;
}

// fire energy minimization with vector of alpha values, individual P values for each cell
void cellPacking2D::fireStep(vector<int>& np, vector<double>& alpha){
	// HARD CODE IN FIRE PARAMETERS
	const double alpha0 	= 0.1;
	const double finc 		= 1.1;
	const double fdec 		= 0.1;
	const double falpha 	= 0.99;
	const double dtmax 		= 10*dt0;
	const int NMIN 			= 5;

	// local variables
	int ci,vi,d;
	double vstarnrm,fstarnrm,vtmp,ftmp;

	// updated temperature
	T = 0;

	// P values
	vector<double> P(NCELLS,0.0);

	// loop over cells
	for (ci=0; ci<NCELLS; ci++){

		// reset norms
		vstarnrm = 0.0;
		fstarnrm = 0.0;

		// loop over vertices
		for (vi=0; vi<cell(ci).getNV(); vi++){
			for (d=0; d<NDIM; d++){
				// get tmp variables
				ftmp = cell(ci).vforce(vi,d);
				vtmp = cell(ci).vvel(vi,d);

				// calculate based on all vertices on all cells
				P.at(ci) += ftmp*vtmp;
				vstarnrm += vtmp*vtmp;
				fstarnrm += ftmp*ftmp;
			}
		}

		// get norms
		vstarnrm = sqrt(vstarnrm);
		fstarnrm = sqrt(fstarnrm);

		// update velocities if forces are acting
		if (fstarnrm > 0){
			for (vi=0; vi<cell(ci).getNV(); vi++){
				for (d=0; d<NDIM; d++){
					vtmp = (1 - alpha.at(ci))*cell(ci).vvel(vi,d) + alpha.at(ci)*(cell(ci).vforce(vi,d)/fstarnrm)*vstarnrm;
					cell(ci).setVVel(vi,d,vtmp);
				}
			}
		}

		// update temperature, cannot use function because we are looping over cells
		for (vi=0; vi<cell(ci).getNV(); vi++){
			for (d=0; d<NDIM; d++){
				T += 0.5*cell(ci).vvel(vi,d)*cell(ci).vvel(vi,d)*0.25*PI*cell(ci).getdel()*cell(ci).getdel();
			}
		}

		// update P and alpha values
		if (P.at(ci) > 0){

			if (np.at(ci) > NMIN){
				// increase dt
				if (dt * finc < dtmax)
					dt *= finc;
				else
					dt = dtmax;

				// decrease alpha
				alpha.at(ci) *= falpha;

				// increase np counter
				np.at(ci)++;
			}
			else
				np.at(ci)++;

		}
		else {
			// decrease time step
			dt *= fdec;

			// set global velocity vector to zero
			for (vi=0; vi<cell(ci).getNV(); vi++){
				for (d=0; d<NDIM; d++)
					cell(ci).setVVel(vi,d,0.0);
			}

			// set alpha -> alpha0
			alpha.at(ci) = alpha0;

			// set np -> 0
			np.at(ci) = 0;
		}
	}
}

// perform single step of a verlet integration with FIRE minimization
void cellPacking2D::fverlet(int& np, double& alpha, double dampingParameter){
	// local variables
	int i;

	// perform FIRE step
	fireStep(np,alpha);

	// update positions
	for (i=0; i<NCELLS; i++){
		cell(i).verletPositionUpdate(dt);
		cell(i).updateCPos();
	}

	// reset contacts before force calculation
	resetContacts();

	// calculate forces
	calculateForces();

	// update velocities
	for (i=0; i<NCELLS; i++)
		cell(i).verletVelocityUpdate(dt,dampingParameter);
}



/**************************

	Simulation Functions

***************************/


// NVE dynamics
void cellPacking2D::cellNVE(){
	// local variables
	int t, ci;
	double U,K;

	// check that NT, NPRINT have been assigned values
	if (NT <= 0){
		cout << "	* ERROR: tried to run tumorNVE() function without proper NT, ending." << endl;
		exit(1);
	}
	else if (NPRINT <= 0){
		cout << "	* ERROR: tried to run tumorNVE() function without proper NPRINT, ending." << endl;
		exit(1);
	}

	// run NVE for allotted time
	for (t=0; t<NT; t++){
		// calculate energies
		U = totalPotentialEnergy();
		K = totalKineticEnergy();

		// print data first to get the initial condition
		if (t % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 	CELL NVE, t = " << t << endl << endl;
			cout << "===================================================" << endl;

			if (packingPrintObject.is_open()){
				cout << "	* Printing vetex positions to file" << endl;
				printSystemPositions();
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing cell energy to file" << endl;
				printSystemEnergy();
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* U 		= " << U << endl;
			cout << "	* K 		= " << K << endl;
			cout << "	* E 		= " << U + K << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << "	* nc 		= " << totalNumberOfContacts() << endl;
			cout << endl << endl;
		}

		// use velocity verlet to advance time

		// update positions
		for (ci=0; ci<NCELLS; ci++){
			cell(ci).verletPositionUpdate(dt);
			cell(ci).updateCPos();
		}

		// reset contacts before force calculation
		resetContacts();

		// calculate forces
		calculateForces();

		// update velocities
		for (ci=0; ci<NCELLS; ci++)
			cell(ci).verletVelocityUpdate(dt);
	}

	// print something to conclude
	cout << "t = NT, which = " << t << ", so ending NVE run!" << endl;
}


// test overdamped dynamics
void cellPacking2D::cellOverDamped(){
	// local variables
	int t, ci, vi, d;
	double U, K;

	// check that NT, NPRINT have been assigned values
	if (NT <= 0){
		cout << "	* ERROR: tried to run tumorNVE() function without proper NT, ending." << endl;
		exit(1);
	}
	else if (NPRINT <= 0){
		cout << "	* ERROR: tried to run tumorNVE() function without proper NPRINT, ending." << endl;
		exit(1);
	}

	// set time step
	dt = 0.1;

	// run overdamped using RK4 for allotted time
	for (t=0; t<NT; t++){
		// calculate energies
		U = totalPotentialEnergy();
		K = totalKineticEnergy();

		// print data first to get the initial condition
		if (t % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 	CELL overdamped dynamics with RK4, t = " << t << endl << endl;
			cout << "===================================================" << endl;

			if (packingPrintObject.is_open()){
				cout << "	* Printing vetex positions to file" << endl;
				printSystemPositions();
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing cell energy to file" << endl;
				printSystemEnergy();
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* U 		= " << U << endl;
			cout << "	* K 		= " << K << endl;
			cout << "	* E 		= " << U + K << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << "	* nc 		= " << totalNumberOfContacts() << endl;
			cout << endl << endl;
		}

		// reset forces before position update
		for (ci=0; ci<NCELLS; ci++){
			// loop over vertices
			for (vi=0; vi<cell(ci).getNV(); vi++){
				// loop over dimensions, update positions and reset forces for next time
				for (d=0; d<NDIM; d++){
					// reset forces
					cell(ci).setVForce(vi,d,0.0);

					// reset interaction energy
					cell(ci).setUInt(vi,0.0);
				}
			}
		}

		// reset contacts before force calculation
		resetContacts();

		// use RK4 to update positions
		gelRK4();
	}

	// print something to conclude
	cout << "t = NT, which = " << t << ", so ending NVE run!" << endl;
}

// LOOPING FUNCTIONS

// prepare jammed packing by ramping asphericity
void cellPacking2D::jammingFireRamp(double dphi, double dCalA, double asphericityTarget, double kbTarget, double phiTarget, double Ktol, double Ptol, int plotIt){
	// local variables
	int i, ci, nr, kr, isjammed, k, kmax, asphericityLow, kbLow;
	double Knew, Pvirial, energyScale, dkb;

	// update packing fraction
	phi = packingFraction();

	// initialize structure to unjammed
	isjammed = 0;
	nr = NCELLS;

	// initialize current calA
	asphericityLow = 1;
	kbLow = 1;

	// set k and kmax
	k = 0;
	kmax = 1e6;

	// bending energy increment
	dkb = 1e-3;

	// initialize energies
	Knew = totalKineticEnergy();
	energyScale = cell(0).getkint();

	// loop over by alternating compression steps and 
	// shape change steps
	while (isjammed == 0 && k < kmax){
		// output to console
		cout << "===================================================" << endl << endl << endl;
		cout << " 	Jamming by isotropic compression" << endl << endl;
		cout << "===================================================" << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* phi 			= " << phi << endl;
		cout << "	* dCalA 		= " << dCalA << endl;
		cout << "	* calA0 (0) 		= " << cell(0).calA0()<< endl;
		cout << "	* calA0 (end) 		= " << cell(NCELLS-1).calA0() << endl;
		cout << "	* low calA 		= " << asphericityLow << endl;
		cout << "	* target calA0 		= " << asphericityTarget << endl;
		cout << "	* # of rattlers = " << nr << endl;
		cout << "	* isjammed = " << isjammed << endl;
		cout << endl << endl;

		// update calA and bending energies if below target
		if (asphericityLow == 1){
			// reset calA check
			asphericityLow = 0;

			// update relevant cal A values
			for (ci=0; ci < NCELLS; ci++){
				if (cell(ci).calA0() < asphericityTarget){
					asphericityLow = 1;
					setAsphericity(ci,cell(ci).calA0()+dCalA);

					// randomize vertices
					cell(ci).vertexPerturbation(0.1);
				}
			}
		}

		if (kbLow == 1){
			// reset kb check
			kbLow = 0;

			// update kb
			for (ci=0; ci < NCELLS; ci++){
				if (cell(ci).getkb() < kbTarget){
					// indicate kb still low
					kbLow = 1;
					
					// increment kb on cell ci
					cell(ci).setkb(cell(ci).getkb() + dkb);
				}
			}
		}

		// relax shapes (energies calculated in relax function)
		// potentialRelaxFire(Ktol,Utol,plotIt,k);
		fireMinimizeP(Ptol, Ktol);
		Knew = totalKineticEnergy();
		Pvirial = 0.5*(sigmaXX + sigmaYY);
		energyScale = cell(0).getkint();

		// remove rattlers
		kr = 0;
		nr = removeRattlers(kr);

		// check for target packing fraction
		phi = packingFraction();
		if (phi > phiTarget){
			cout << "	** phi > phiTarget, so ending!" << endl;
			cout << "	** Final phi = " << phi << endl;
			break;
		}
		else
			setPackingFraction(phi + dphi);

		// check whether or not system has jammed
		if (abs(Pvirial) > 2*Ptol*energyScale && Knew < NCELLS*cell(0).getNV()*energyScale*Ktol)
			isjammed = 1;

		if (isjammed){

			// increase cal A for all cells below target calA0
			if (asphericityLow == 1){

				// reset jamming check
				isjammed = 0;

				// loop over cells until target asphericity found
				while (asphericityLow == 1){
					// reset asphericityLow
					asphericityLow = 0;

					// increase iterator
					k++;

					// update asphericities
					for (ci=0; ci < NCELLS; ci++){
						if (cell(ci).calA0() < asphericityTarget){
							asphericityLow = 1;
							setAsphericity(ci,cell(ci).calA0()+dCalA);
						}
					}

					// relax energy
					fireMinimizeP(Ptol, Ktol);
				}
			}
			else{
				cout << "Mechanically-stable state found!" << endl;
				cout << "Note that this is just MS, probably overjammed, need root search to access true jamming point" << endl;
				cout << "Final phi = " << phi << endl;
				cout << "Writing final configuration to file." << endl;
				printSystemPositions(k++);
				packingPrintObject << setw(12) << left << "STERM" << " " << endl;
				break;
			}
		}
	}

	if (k == kmax){
		cout << "	ERROR: particles could not jam in allotted time, ending." << endl;
		exit(1);
	}
	else{
		cout << " Jamming a success! End asphericity is = " << meanAsphericity() << endl;
		cout << " final k = " << k << endl;
	}
}



// compress packing to target phi quasistatically
void cellPacking2D::compressToTarget(double dphi, double phiTarget, double asphericityTarget, double Ktol, double Ptol, int plotIt, int& frameCount){
	// local variables
	int i, ci, k, kmax, NPHISTEPS;
	double dCalA, dphiTmp;

	// update packing fraction
	phi = packingFraction();

	// determine calA steps to have asphericity increase before phi target is met
	NPHISTEPS = floor((phiTarget - phi)/dphi) - 1;
	dCalA = (asphericityTarget - cell(NCELLS-1).calA0())/NPHISTEPS;

	// set k and kmax
	frameCount = 0;
	k = 0;
	kmax = 1e5;

	// loop over by alternating compression steps and 
	// shape change steps
	while (phi <  phiTarget && k < kmax){
		// output to console
		cout << "===================================================" << endl << endl << endl;
		cout << " 	Jamming by isotropic compression" << endl << endl;
		cout << "===================================================" << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* phi 			= " << phi << endl;
		cout << "	* dCalA 		= " << dCalA << endl;
		cout << "	* calA0 (0) 		= " << cell(0).calA0()<< endl;
		cout << "	* calA0 (end) 		= " << cell(NCELLS-1).calA0() << endl;
		cout << "	* target calA0 		= " << asphericityTarget << endl;
		cout << endl << endl;

		// relax shapes (energies calculated in relax function)
		fireMinimizeP(Ptol, Ktol);
		k = frameCount;

		// increase asphericity
		for (ci=0; ci<NCELLS; ci++){
			// if below asphericity target, increase
			if (cell(ci).calA0() < asphericityTarget)
				setAsphericity(ci,cell(ci).calA0()+dCalA);
		}

		// increase packing fraction
		phi = packingFraction();
		if (phi + dphi < phiTarget)
			setPackingFraction(phi + dphi);
		else{
			// use dummy phi step to take smaller incremental steps
			dphiTmp = dphi;

			// loop until step is below target
			while (phi + dphiTmp > phiTarget && dphiTmp > 1e-6)
				dphiTmp *= 0.5;

			// increase phi
			setPackingFraction(phi + dphiTmp);
		}

		if (phi > phiTarget){
			cout << "	** phi > phiTarget, so ending!" << endl;
			cout << "	** Final phi = " << phi << endl;
			break;
		}
		else
			setPackingFraction(phi + dphi);
	}

	if (k == kmax){
		cout << "	ERROR: particles could not jam in allotted time, ending." << endl;
		exit(1);
	}
	else{
		cout << " Compression a success! End asphericity is = " << meanAsphericity() << endl;
		cout << " final k = " << k << endl;
	}
}


// FIRE 2.0 energy minimzation with backstepping if P < 0
void cellPacking2D::fireMinimizeP(double Ptol, double Ktol){
	// HARD CODE IN FIRE PARAMETERS
	const double alpha0 	= 0.25;
	const double finc 		= 1.05;
	const double fdec 		= 0.5;
	const double falpha 	= 0.99;
	const double dtmax 		= 10*dt0;
	const double dtmin 		= 0.05*dt0;
	const int NMIN 			= 20;
	const int NNEGMAX 		= 2000;
	const int NDELAY 		= 1000;
	int npPos				= 0;
	int npNeg 				= 0;
	int npPMIN				= 0;
	double alpha 			= alpha0;
	double alphat 			= alpha;
	double t 				= 0.0;
	double forceScale 		= cell(0).getkint();
	double energyScale		= forceScale;
	double P 				= 0;

	// local variables
	int ci,vi,d,k,kmax;
	double vstarnrm,fstarnrm,vtmp,ftmp;
	double Knew, Pvirial, Kcheck, Pcheck;

	// variable to test for potential energy minimization
	bool converged = false;

	// reset time step
	dt = dt0;

	// initialize forces
	resetContacts();
	calculateForces();

	// initialize virial pressure from pressure from last time
	Pvirial = 0.5*(sigmaXX + sigmaYY);

	// initialize energy and force tracking, pressure
	Knew = totalKineticEnergy();

	// scale P and K for convergence checking
	Pcheck = Pvirial/(energyScale*NCELLS*cell(0).getNV());
	Kcheck = Knew/(energyScale*NCELLS*cell(0).getNV());

	// reset velocities to 0
	for (ci=0; ci<NCELLS; ci++){
		for (vi=0; vi<cell(ci).getNV(); vi++){
			for (d=0; d<NDIM; d++)
				cell(ci).setVVel(vi,d,0.0);
		}
	}

	// iterate through MD time until system converged
	kmax = 5e5;
	for (k=0; k<kmax; k++){

		// output some information to console
		if (k % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 	FIRE MINIMIZATION, k = " << k << endl << endl;
			cout << "===================================================" << endl;			
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << Kcheck << endl;
			cout << "	* Pvirial 	= " << Pcheck << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << "	* alpha 	= " << alpha << endl;
			cout << "	* alphat 	= " << alphat << endl;
			cout << "	* P 		= " << P << endl;
			cout << endl << endl;
		}

		// Step 1. calculate P and norms
		P = 0.0;
		vstarnrm = 0.0;
		fstarnrm = 0.0;
		for (ci=0; ci<NCELLS; ci++){
			for (vi=0; vi<cell(ci).getNV(); vi++){
				for (d=0; d<NDIM; d++){
					// get tmp variables
					ftmp = cell(ci).vforce(vi,d);
					vtmp = cell(ci).vvel(vi,d);

					// calculate based on all vertices on all cells
					P += ftmp*vtmp;
					vstarnrm += vtmp*vtmp;
					fstarnrm += ftmp*ftmp;
				}
			}
		}

		// get norms
		vstarnrm = sqrt(vstarnrm);
		fstarnrm = sqrt(fstarnrm);


		// Step 2. Adjust simulation based on net motion of system
		if (P > 0){
			// increment pos counter
			npPos++;

			// reset neg counter
			npNeg = 0;

			// update alpha_t for next time
			alphat = alpha;

			// alter sim if enough positive steps taken
			if (npPos > NMIN){
				// change time step
				if (dt*finc < dtmax)
					dt *= finc;
				else
					dt = dtmax;

				// decrease alpha
				alpha *= falpha;
			}
		}
		else{
			// reset pos counter
			npPos = 0;

			// rest neg counter
			npNeg++;

			// check for stuck sim
			if (npNeg > NNEGMAX)
				break;

			// decrease time step if past initial delay
			if (k > NMIN){
				// decrease time step 
				if (dt*fdec > dtmin)
					dt *= fdec;
				else
					dt = dtmin;

				// change alpha
				alpha = alpha0;
				alphat = alpha;
			}

			// take half step backwards
			for (ci=0; ci<NCELLS; ci++){
				for (vi=0; vi<cell(ci).getNV(); vi++){
					for (d=0; d<NDIM; d++)
						cell(ci).setVPos(vi,d,cell(ci).vpos(vi,d) - 0.5*dt*cell(ci).vvel(vi,d));
				}
			}

			// reset velocities to 0
			for (ci=0; ci<NCELLS; ci++){
				for (vi=0; vi<cell(ci).getNV(); vi++){
					for (d=0; d<NDIM; d++)
						cell(ci).setVVel(vi,d,0.0);
				}
			}
		}

		// update velocities if forces are acting
		if (fstarnrm > 0){
			for (ci=0; ci<NCELLS; ci++){
				for (vi=0; vi<cell(ci).getNV(); vi++){
					for (d=0; d<NDIM; d++){
						vtmp = (1 - alphat)*cell(ci).vvel(vi,d) + alphat*(cell(ci).vforce(vi,d)/fstarnrm)*vstarnrm;
						cell(ci).setVVel(vi,d,vtmp);
					}
				}
			}
		}

		// do verlet update
		for (ci=0; ci<NCELLS; ci++){
			cell(ci).verletPositionUpdate(dt);
			cell(ci).updateCPos();
		}

		// reset contacts before force calculation
		resetContacts();

		// calculate forces
		calculateForces();

		// update velocities
		for (ci=0; ci<NCELLS; ci++)
			cell(ci).verletVelocityUpdate(dt);

		// update t
		t += dt;

		// track energy and forces
		Knew = totalKineticEnergy();
		Pvirial = 0.5*(sigmaXX + sigmaYY);

		// scale P and K for convergence checking
		Pcheck = Pvirial/(energyScale*NCELLS*cell(0).getNV());
		Kcheck = Knew/(energyScale*NCELLS*cell(0).getNV());

		// update if Pvirial under tol
		if (abs(Pcheck) < Ptol)
			npPMIN++;
		else
			npPMIN = 0;

		// check for convergence
		converged = (abs(Pcheck) < Ptol && npPMIN > NMIN);
		converged = (converged || (abs(Pcheck) > Ptol && Kcheck < Ktol));

		if (converged){
			cout << "	** FIRE has converged!" << endl;
			cout << "	** Kcheck = " << Kcheck << endl;
			cout << "	** Pcheck = " << Pcheck << endl;
			cout << "	** k = " << k << ", t = " << t << endl;
			cout << "	** Breaking out of FIRE protocol." << endl;

			// print config and energy
			if (packingPrintObject.is_open()){
				cout << "	* Printing vetex positions to file" << endl;
				printSystemPositions();
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing cell energy to file" << endl;
				printSystemEnergy();
			}

			if (statPrintObject.is_open()){
				cout << "	* Printing cell contacts to file" << endl;
				printSystemContacts();
			}

			break;
		}
	}

	// reset dt to be original value before ending function
	dt = dt0;

	// if no convergence, just stop
	if (k == kmax)
		cout << "	** FIRE not converged in kmax = " << kmax << " force evaluations" << endl;
}


void cellPacking2D::fireMinimizeF(double Ftol, double Ktol, int plotIt, int& frameCount){
	// HARD CODE IN FIRE PARAMETERS
	const double alpha0 	= 0.25;
	const double finc 		= 1.1;
	const double fdec 		= 0.5;
	const double falpha 	= 0.99;
	const double dtmax 		= 10*dt0;
	const double dtmin 		= 0.02*dt0;
	const int NMIN 			= 20;
	const int NNEGMAX 		= 2000;
	const int NDELAY 		= 1000;
	int npPos				= 0;
	int npNeg 				= 0;
	int npPMIN				= 0;
	double alpha 			= alpha0;
	double alphat 			= alpha;
	double t 				= 0.0;
	double energyScale 		= cell(0).getkint();
	double forceScale		= energyScale/cell(0).getl0();

	// local variables
	int ci,vi,d,k,kmax;
	double P,vstarnrm,fstarnrm,vtmp,ftmp;
	double Knew, fRMS;

	// variable to test for potential energy minimization
	bool converged = false;

	// reset time step
	dt = dt0;

	// initialize root-mean-square force from last time
	fRMS = forceRMS();;

	// initialize forces
	resetContacts();
	calculateForces();

	// initialize energy and force tracking, pressure
	Knew = totalKineticEnergy();

	// reset velocities to 0
	for (ci=0; ci<NCELLS; ci++){
		for (vi=0; vi<cell(ci).getNV(); vi++){
			for (d=0; d<NDIM; d++)
				cell(ci).setVVel(vi,d,0.0);
		}
	}

	// iterate through MD time until system converged
	kmax = 5e5;
	for (k=0; k<kmax; k++){

		// output some information to console
		if (k % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 	FIRE MINIMIZATION, k = " << k << ", frame = " << frameCount << endl << endl;
			cout << "===================================================" << endl;

			if (plotIt == 1){
				if (packingPrintObject.is_open()){
					cout << "	* Printing vetex positions to file" << endl;
					printSystemPositions(frameCount);
				}
				
				if (energyPrintObject.is_open()){
					cout << "	* Printing cell energy to file" << endl;
					printSystemEnergy(frameCount,fRMS/forceScale,Knew/energyScale);
				}
				frameCount++;
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << Knew/(NCELLS*cell(0).getNV()*energyScale*Ktol) << endl;
			cout << "	* fRMS 		= " << fRMS/(Ftol*forceScale) << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << "	* alpha 	= " << alpha << endl;
			cout << "	* P 		= " << P << endl;
			cout << endl << endl;
		}

		// Step 1. calculate P and norms
		P = 0.0;
		vstarnrm = 0.0;
		fstarnrm = 0.0;
		for (ci=0; ci<NCELLS; ci++){
			for (vi=0; vi<cell(ci).getNV(); vi++){
				for (d=0; d<NDIM; d++){
					// get tmp variables
					ftmp = cell(ci).vforce(vi,d);
					vtmp = cell(ci).vvel(vi,d);

					// calculate based on all vertices on all cells
					P += ftmp*vtmp;
					vstarnrm += vtmp*vtmp;
					fstarnrm += ftmp*ftmp;
				}
			}
		}

		// get norms
		vstarnrm = sqrt(vstarnrm);
		fstarnrm = sqrt(fstarnrm);


		// Step 2. Adjust simulation based on net motion of system
		if (P > 0){
			// increment pos counter
			npPos++;

			// reset neg counter
			npNeg = 0;

			// update alpha_t for next time
			alphat = alpha;

			// alter sim if enough positive steps taken
			if (npPos > NMIN){
				// change time step
				if (dt*finc < dtmax)
					dt *= finc;
				else
					dt = dtmax;

				// decrease alpha
				alpha *= falpha;
			}
		}
		else{
			// reset pos counter
			npPos = 0;

			// rest neg counter
			npNeg++;

			// check for stuck sim
			if (npNeg > NNEGMAX)
				break;

			// decrease time step if past initial delay
			if (k > NMIN){
				// decrease time step 
				if (dt*fdec > dtmin)
					dt *= fdec;
				else
					dt = dtmin;

				// change alpha
				alpha = alpha0;
				alphat = alpha;
			}

			// take half step backwards
			for (ci=0; ci<NCELLS; ci++){
				for (vi=0; vi<cell(ci).getNV(); vi++){
					for (d=0; d<NDIM; d++)
						cell(ci).setVPos(vi,d,cell(ci).vpos(vi,d) - 0.5*dt*cell(ci).vvel(vi,d));
				}
			}

			// reset velocities to 0
			for (ci=0; ci<NCELLS; ci++){
				for (vi=0; vi<cell(ci).getNV(); vi++){
					for (d=0; d<NDIM; d++)
						cell(ci).setVVel(vi,d,0.0);
				}
			}
		}

		// update velocities if forces are acting
		if (fstarnrm > 0){
			for (ci=0; ci<NCELLS; ci++){
				for (vi=0; vi<cell(ci).getNV(); vi++){
					for (d=0; d<NDIM; d++){
						vtmp = (1 - alphat)*cell(ci).vvel(vi,d) + alphat*(cell(ci).vforce(vi,d)/fstarnrm)*vstarnrm;
						cell(ci).setVVel(vi,d,vtmp);
					}
				}
			}
		}

		// do verlet update
		for (ci=0; ci<NCELLS; ci++){
			cell(ci).verletPositionUpdate(dt);
			cell(ci).updateCPos();
		}

		// reset contacts before force calculation
		resetContacts();

		// calculate forces
		calculateForces();

		// update velocities
		for (ci=0; ci<NCELLS; ci++)
			cell(ci).verletVelocityUpdate(dt,0.0);

		// update t
		t += dt;

		// track energy and forces
		Knew = totalKineticEnergy();
		fRMS = forceRMS();
		energyScale = cell(0).getkint();
		forceScale = energyScale/cell(0).getdel();

		// update if root-mean-square force under tol
		if (abs(fRMS) < Ftol*energyScale)
			npPMIN++;
		else
			npPMIN = 0;

		// check for convergence
		converged = (abs(fRMS) < Ftol*forceScale && npPMIN > NMIN);
		converged = (converged || (abs(fRMS) > 2*Ftol*forceScale && Knew < NCELLS*cell(0).getNV()*energyScale*Ktol));

		if (converged){
			cout << "	** FIRE has converged!" << endl;
			cout << "	** k = " << k << ", t = " << t << endl;
			cout << "	** rms F = " << fRMS << endl;
			cout << "	** Breaking out of FIRE protocol." << endl;
			break;
		}
	}

	// reset dt to be original value before ending function
	dt = dt0;

	// if no convergence, just stop
	if (k == kmax)
		cout << "	** FIRE not converged in kmax = " << kmax << " force evaluations" << endl;
}


// GELATION FUNCTIONS


// function to test isotropic extension protocol
void cellPacking2D::twoParticleContact(int NV){
	// local variables
	int ci, vi, d, nvtmp;
	double rtmp, l0tmp, a0tmp, p0tmp;

	// box length from 3 particle diameters (each has unit radius)
	L = 6.0;

	// number of vertices on both particles
	nvtmp = NV;

	// generate cells and sizes
	for (ci=0; ci<NCELLS; ci++){
		// boundary information ( SET PBCS TO 1 for plant cells )
		cell(ci).setL(L);
		cell(ci).setpbc(0,1);
		cell(ci).setpbc(1,1);

		// set number of vertices
		cell(ci).setNV(nvtmp);

		// array information
		cell(ci).initializeVertices();
		cell(ci).initializeCell();

		// unit radius
		rtmp = 1.0;

		// calculate a0 and l0 based on fact that they are regular polygons
		a0tmp = 0.5*nvtmp*pow(rtmp,2.0)*sin(2.0*PI/nvtmp);
		l0tmp = 2.0*rtmp*sin(PI/nvtmp);
		p0tmp = l0tmp*nvtmp;

		// set preferred area and length 
		cell(ci).seta0(a0tmp);
		cell(ci).setl0(l0tmp);
		cell(ci).setdel(1.0);
	}

	// initialize particle positions
	cout << "		-- Ininitializing cell positions" << endl;
	// set as initial position of com
	cell(0).setCPos(0,0.5*L - 1.0 - 0.5*cell(0).getl0());
	cell(0).setCPos(1,0.5*L);

	cell(1).setCPos(0,0.5*L + 1.0 + 0.5*cell(0).getl0());
	cell(1).setCPos(1,0.5*L);

	// initialize vertices as a regular polygon
	for (ci=0; ci<NCELLS; ci++)
		cell(ci).regularPolygon();

	// initial time scales ( = sqrt(m*sigma/f_0) = sqrt(PI)))
	cout << "		-- Ininitializing time scale" << endl;

	// packing fraction
	phi = packingFraction();

	// print config and energy
	if (packingPrintObject.is_open()){
		cout << "	* Printing vetex positions to file" << endl;
		printSystemPositions();
	}
	
	if (energyPrintObject.is_open()){
		cout << "	* Printing cell energy to file" << endl;
		printSystemEnergy();
	}

	if (statPrintObject.is_open()){
		cout << "	* Printing cell contacts to file" << endl;
		printSystemContacts();
	}
}

// initialize plant cell particles as disks at input packing fraction
void cellPacking2D::initializeGel(int NV, double phiDisk, double sizeDispersion, double delval){
	// local variables
	int ci, vi, d, nvtmp;
	double calA;
	double xpos, ypos;
	double xmin, xmax, ymin, ymax;
	double r1, r2, g1, radsum;
	double rtmp, a0tmp, l0tmp, p0tmp, calA0tmp;
	double meanArea = 0.0;
	double lc = 0.0;

	// seed random number generator
	srand48(23562457*seed);

	// minimum number of vertices
	const int nvmin = 12;

	// output to console
	cout << "		-- In gelation initialization, initializing plant cells and relaxing initial overlaps as repulsive SP particles" << endl;

	// initialize length scales as gaussian random variables (will becomes area square roots)
	radsum = 0.0;
	vector<double> lenscales(NCELLS,0.0);
	for (ci=0; ci<NCELLS; ci++){
		// generate random numbers
		r1 = drand48();
		r2 = drand48();

		// calculate gaussian random variable using Box-Muller transform
		g1 = sqrt(-2.0*log(r1))*cos(2*PI*r2);

		// get radius
		lenscales.at(ci) = g1*sizeDispersion + 1.0;

		// add to lenscales sum for boundary size
		radsum += lenscales.at(ci)*lenscales.at(ci);
	}

	// determine box length from particle sizes and input packing fraction
	L = sqrt(PI*radsum/phiDisk);

	// set phi to input
	phi = phiDisk;

	// reseed rng
	srand48(56835698*seed);

	// initialize cell information
	cout << "		-- Ininitializing plant cell objects" << endl;
	for (ci=0; ci<NCELLS; ci++){
		// boundary information ( SET PBCS TO 1 for plant cells )
		cell(ci).setL(L);
		cell(ci).setpbc(0,1);
		cell(ci).setpbc(1,1);

		// number of vertices ( SIGMA SETS # OF VERTS )
		nvtmp = ceil(lenscales.at(ci)*NV);
		if (nvtmp > nvmin)
 			cell(ci).setNV(nvtmp);
		else
			cell(ci).setNV(nvmin);

		// array information
		cell(ci).initializeVertices();
		cell(ci).initializeCell();

		// initial length of polygon side
		l0tmp = 2.0*lenscales.at(ci)*sin(PI/nvtmp);

		// use rtmp slightly smaller than lenscale, so no overlaps at end
		rtmp = lenscales.at(ci) - 0.25*delval*l0tmp;

		// calculate a0 and l0 based on fact that they are regular polygons
		a0tmp = 0.5*nvtmp*pow(rtmp,2.0)*sin(2.0*PI/nvtmp);
		l0tmp = 2.0*rtmp*sin(PI/nvtmp);
		p0tmp = l0tmp*nvtmp;

		// set preferred area and length 
		cell(ci).seta0(a0tmp);
		cell(ci).setl0(l0tmp);
		cell(ci).setdel(delval);
	}

	// initialize particle positions
	cout << "		-- Ininitializing cell positions" << endl;
	for (ci=0; ci<NCELLS; ci++){
		// set min and max values of positions
		xmin = 0;
		xmax = L;
		ymin = 0;
		ymax = L;

		// get random location in box
		xpos = (xmax-xmin)*drand48() + xmin;
		ypos = (ymax-ymin)*drand48() + ymin;

		// set as initial position of com
		cell(ci).setCPos(0,xpos);
		cell(ci).setCPos(1,ypos);

		// initialize vertices as a regular polygon
		cell(ci).regularPolygon();
	}

	// initial time scales ( = sqrt(m*sigma/f_0) = sqrt(PI)))
	cout << "		-- Ininitializing time scale" << endl;
	dt = 0.1*sqrt(PI);
	dt0 = dt;

	// use FIRE in PBC box to relax overlaps
	cout << "		-- Using FIRE to relax overlaps..." << endl;
	fireMinimizeSP(lenscales);
}

// set forces values for plant cells
void cellPacking2D::gelForceVals(double calA0, double kl, double ka, double gam, double kb, double kint, double del, double a){
	// local variables
	int ci;

	// loop over cells, add values
	for (ci=0; ci<NCELLS; ci++){
		// set calA0 by altering a0
		cell(ci).setAsphericity(calA0);

		// set shape force scales
		cell(ci).setkl(kl);
		cell(ci).setka(ka);
		cell(ci).setgam(gam);
		cell(ci).setkb(kb);

		// set interaction force scales
		cell(ci).setkint(kint);
		cell(ci).setdel(del);
		cell(ci).seta(a);
	}
}


// compress isostatically
void cellPacking2D::qsIsoCompression(double phiTarget, double deltaPhi){
	// local variables
	double phi, phi0, phiNew, dphi;
	int NSTEPS, k;

	// tolerances
	const double Ktol = 1e-15;
	const double Ptol = 1e-6;

	// relax shapes (energies calculated in relax function)
	cout << "	** IN qsIsoCompression, performing initial relaxation" << endl;
	fireMinimizeP(Ptol, Ktol);

	// get initial packing fraction
	phi = packingFraction();
	phi0 = phi;

	// determine number of steps to target
	NSTEPS = ceil((phiTarget - phi0)/deltaPhi) + 1;

	// update new dphi to make steps even
	dphi = (phiTarget - phi)/NSTEPS;

	// iterator
	k = 0;

	// loop until phi is the correct value
	while (k < NSTEPS){
		// update iterator
		k++;

		// output to console
		cout << "===================================================" << endl << endl << endl;
		cout << " 	quasistatic isotropic compression to target phi = " << phiTarget << endl << endl;
		cout << "===================================================" << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* NSTEPS 		= " << NSTEPS << endl;		
		cout << "	* dphi 			= " << dphi << endl << endl;
		cout << "	AFTER LAST MINIMIZATION:" << endl;
		cout << "	* phi 			= " << phi << endl;
		cout << endl << endl;

		// increase packing fraction to new phi value
		phiNew = phi0 + k*dphi;
		setPackingFraction(phiNew);

		// relax shapes (energies calculated in relax function)
		fireMinimizeP(Ptol, Ktol);

		// calculate phi after minimization
		phi = packingFraction();
	}
}

// increase attraction quasi-statically, relax after each increment
void cellPacking2D::attractionRamp(double attractionTarget, double dAttraction){
	// local variables
	int k, ci, NSTEPS;
	double da, currAttraction;

	// tolerances
	const double Ktol = 1e-12;
	const double Ptol = 1e-6;

	// relax shapes (energies calculated in relax function)
	cout << "	** IN attractionRamp, performing initial relaxation" << endl;
	fireMinimizeP(Ptol, Ktol);

	// get current attraction (assuming all attraction the same)
	currAttraction = cell(0).geta();

	// determine step for changing asphericity
	NSTEPS = ceil((attractionTarget - currAttraction)/dAttraction) + 1;

	// update new dphi to make steps even
	da = (attractionTarget - currAttraction)/NSTEPS;

	// iterator
	k = 0;

	// loop until attraction is the correct value
	while (k < NSTEPS){
		// update iterator
		k++;

		// output to console
		cout << "===================================================" << endl << endl << endl;
		cout << " 	Ramping attraction to aTarget = " << attractionTarget << endl << endl;
		cout << "===================================================" << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* NSTEPS 		= " << NSTEPS << endl;
		cout << "	* current a 	= " << currAttraction << endl;
		cout << "	* next a 		= " << currAttraction + da << endl;
		cout << endl << endl;

		// update attraction in all of the cells
		for (ci=0; ci<NCELLS; ci++)
			cell(ci).seta(currAttraction + da);

		// update current attraction
		currAttraction += da;

		// relax potential energy
		fireMinimizeP(Ptol, Ktol);
	}
}

// decrease packing fraction at rate
void cellPacking2D::gelRateExtension(double phiGel, double gelRate, double timeStepMag){
	// local variables
	int ci, vi, d, k;
	double t = 0.0;
	double phi0, phitmp, veltmp;
	cellPacking2D checkObj;

	// damping
	const double b = 0.5;

	// update time step
	dt = timeStepMag*sqrt(PI);
	dt0 = dt;

	// get initial packing fraction
	phi0 = packingFraction();
	phitmp = phi0;

	// iterator
	k = 0;

	// loop time until packing fraction below threshold
	while (phitmp > phiGel){
		// update phi
		phi = phitmp;
		phitmp = phi0/(gelRate*t + 1.0);
		setPackingFraction(phitmp);

		// update iterator
		k++;

		// print to console/file
		if (k % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 		gel rate sim, k = " << k << ", t = " << t << endl << endl;
			cout << "===================================================" << endl;

			// print config and energy
			if (packingPrintObject.is_open()){
				cout << "	* Printing vetex positions to file" << endl;
				printSystemPositions();
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing cell energy to file" << endl;
				printSystemEnergy();
			}

			if (statPrintObject.is_open()){
				cout << "	* Printing cell contacts to file" << endl;
				printSystemContacts();
			}
			
			cout << "	* phitmp 	= " << phitmp << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* K 		= " << totalKineticEnergy() << endl;
			cout << "	* U 		= " << totalPotentialEnergy() << endl;
			cout << "	* nc 		= " << totalNumberOfContacts() << endl;
			cout << endl << endl;
		}

		// // DEBUG: PRINT OUT EVERY TIME FOR k > 78000
		// if (k > 78000){
		// 	cout << "Printing everything, nan coming up soon!" << endl;

		// 	// print config and energy
		// 	if (packingPrintObject.is_open()){
		// 		cout << "	* Printing vetex positions to file" << endl;
		// 		printSystemPositions();
		// 	}
			
		// 	if (energyPrintObject.is_open()){
		// 		cout << "	* Printing cell energy to file" << endl;
		// 		printSystemEnergy();
		// 	}

		// 	if (k > 79000){
		// 		cout << "should have found a nan by now, ending" << endl;
		// 		exit(1);
		// 	}

		// }

		// do verlet update
		for (ci=0; ci<NCELLS; ci++){
			cell(ci).verletPositionUpdate(dt);
			cell(ci).updateCPos();
		}

		// reset contacts before force calculation
		resetContacts();

		// calculate forces
		calculateForces();

		// update velocities
		for (ci=0; ci<NCELLS; ci++)
			cell(ci).verletVelocityUpdate(dt,b);

		/*
		// reset forces and uint for RK4 update below
		for (ci=0; ci<NCELLS; ci++){
			// loop over vertices
			for (vi=0; vi<cell(ci).getNV(); vi++){
				// loop over dimensions, update positions and reset forces for next time
				for (d=0; d<NDIM; d++){
					// get velocities (= forces in overdamped regime)
					// veltmp = cell(ci).vforce(vi,d);

					// // update positions (EULER STEP)
					// cell(ci).setVPos(vi,d,cell(ci).vpos(vi,d) + dt*veltmp);

					// update velocities
					cell(ci).setVVel(vi,d,veltmp);

					// reset forces
					cell(ci).setVForce(vi,d,0.0);

					// reset interaction energy
					cell(ci).setUInt(vi,0.0);
				}
			}
		}

		// reset contacts before force calculation
		resetContacts();

		// use RK4 to update positions
		// calculateForces();
		gelRK4();
		*/

		// increment
		t += dt;
	}
}


void cellPacking2D::gelRK4(){
	// local variables
	int ci, vi, d;
	double veltmp;
	cellPacking2D k2, k3, k4;

	// Update positions based on current forces
	calculateForces();

	// update positions for k2
	k2 = *this;
	for (ci=0; ci<k2.NCELLS; ci++){
		for (vi=0; vi<k2.cell(ci).getNV(); vi++){
			for (d=0; d<k2.NDIM; d++)
				k2.cell(ci).setVPos(vi,d,cell(ci).vpos(vi,d) + 0.5*dt*cell(ci).vforce(vi,d));
		}
	}

	// calculate forces due to new positions in k2
	k2.calculateForces();

	// update positions for k3
	k3 = k2;
	for (ci=0; ci<k3.NCELLS; ci++){
		for (vi=0; vi<k3.cell(ci).getNV(); vi++){
			for (d=0; d<k3.NDIM; d++)
				k3.cell(ci).setVPos(vi,d,cell(ci).vpos(vi,d) + 0.5*dt*k2.cell(ci).vforce(vi,d));
		}
	}

	// calculate forces due to new positions in k3
	k3.calculateForces();

	// update positions for k4
	k4 = k3;
	for (ci=0; ci<k4.NCELLS; ci++){
		for (vi=0; vi<k4.cell(ci).getNV(); vi++){
			for (d=0; d<k4.NDIM; d++)
				k4.cell(ci).setVPos(vi,d,cell(ci).vpos(vi,d) + dt*k3.cell(ci).vforce(vi,d));
		}
	}

	// update positions in this based on RK4
	for (ci=0; ci<k2.NCELLS; ci++){
		for (vi=0; vi<k2.cell(ci).getNV(); vi++){
			for (d=0; d<k2.NDIM; d++){
				veltmp = (1.0/6.0)*(cell(ci).vforce(vi,d) + 2.0*k2.cell(ci).vforce(vi,d) + 2.0*k3.cell(ci).vforce(vi,d) + k4.cell(ci).vforce(vi,d));
				cell(ci).setVPos(vi,d,cell(ci).vpos(vi,d) + dt*veltmp);
				cell(ci).setVVel(vi,d,veltmp);
			}
		}
	}
}




// FUNCTIONS FOR INTIIAL SP RELAXATION STAGE
void cellPacking2D::fireMinimizeSP(vector<double>& lenscales){
	// HARD CODE IN FIRE PARAMETERS
	const double alpha0 	= 0.25;
	const double finc 		= 1.01;
	const double fdec 		= 0.5;
	const double falpha 	= 0.99;
	const double dtmax 		= 10*dt0;
	const double dtmin 		= 0.02*dt0;
	const int NMIN 			= 20;
	const int NNEGMAX 		= 2000;
	const int NDELAY 		= 1000;
	int npPos				= 0;
	int npNeg 				= 0;
	int npPMIN				= 0;
	int closed 				= 1;
	double alpha 			= alpha0;
	double alphat 			= alpha;
	double t 				= 0.0;
	double Ptol 			= 1e-8;
	double Ktol 			= 1e-12;
	bool converged 			= false;

	// local variables
	int ci,vi,d,itr,itrMax;
	double P,vstarnrm,fstarnrm,vtmp,ftmp,ptmp;
	double Knew, Pvirial;
	double Kcheck, Pcheck;

	// reset time step
	dt = dt0;

	// reset velocities to 0
	for (ci=0; ci<NCELLS; ci++){
		for (d=0; d<NDIM; d++){
			cell(ci).setCVel(d,0.0);
			cell(ci).setCForce(d,0.0);
		}
	}

	// initialize forces (neglect damping forces, only interactions)
	resetContacts();
	spForces(lenscales);

	// initialize virial pressure from pressure from last time
	Pvirial = 0.5*(sigmaXX + sigmaYY);

	// update kinetic energy based on com velocity
	Knew = 0.0;
	for (ci=0; ci<NCELLS; ci++)
		Knew += 0.5*(PI*pow(lenscales.at(ci),2))*sqrt(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2));

	// iterate through MD time until system converged
	itrMax = 5e5;
	for (itr=0; itr<itrMax; itr++){

		// output some information to console
		if (itr % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 	FIRE MINIMIZATION, itr = " << itr << endl << endl;
			cout << "===================================================" << endl;

			// print if object has been opened already
			if (packingPrintObject.is_open()){
				cout << "	* Printing SP center positions to file" << endl;
				printSystemPositions();
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing SP energy to file" << endl;
				printSystemEnergy();
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << Knew << endl;
			cout << "	* Pvirial 	= " << Pvirial << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << "	* alpha 	= " << alpha << endl;
			cout << "	* P 		= " << P << endl;
			cout << endl << endl;
		}

		// Step 1. calculate P and norms
		P = 0.0;
		vstarnrm = 0.0;
		fstarnrm = 0.0;
		for (ci=0; ci<NCELLS; ci++){
			for (d=0; d<NDIM; d++){
				// get tmp variables
				ftmp = cell(ci).cforce(d);
				vtmp = cell(ci).cvel(d);

				// calculate based on all vertices on all cells
				P += ftmp*vtmp;
				vstarnrm += vtmp*vtmp;
				fstarnrm += ftmp*ftmp;
			}
		}

		// get norms
		vstarnrm = sqrt(vstarnrm);
		fstarnrm = sqrt(fstarnrm);


		// Step 2. Adjust simulation based on net motion of system
		if (P > 0){
			// increment pos counter
			npPos++;

			// reset neg counter
			npNeg = 0;

			// update alpha_t for next time
			alphat = alpha;

			// alter sim if enough positive steps taken
			if (npPos > NMIN){
				// change time step
				if (dt*finc < dtmax)
					dt *= finc;
				else
					dt = dtmax;

				// decrease alpha
				alpha *= falpha;
			}
		}
		else{
			// reset pos counter
			npPos = 0;

			// rest neg counter
			npNeg++;

			// check for stuck sim
			if (npNeg > NNEGMAX)
				break;

			// decrease time step if past initial delay
			if (itr > NMIN){
				// decrease time step 
				if (dt*fdec > dtmin)
					dt *= fdec;
				else
					dt = dtmin;

				// change alpha
				alpha = alpha0;
				alphat = alpha;
			}

			// take half step backwards
			for (ci=0; ci<NCELLS; ci++){
				for (d=0; d<NDIM; d++)
					cell(ci).setCPos(d,cell(ci).cpos(d) - 0.5*dt*cell(ci).cvel(d));
			}

			// reset velocities to 0
			for (ci=0; ci<NCELLS; ci++){
				for (d=0; d<NDIM; d++)
					cell(ci).setCVel(d,0.0);
			}
		}

		// update velocities if forces are acting
		if (fstarnrm > 0){
			for (ci=0; ci<NCELLS; ci++){
				for (d=0; d<NDIM; d++){
					vtmp = (1 - alphat)*cell(ci).cvel(d) + alphat*(cell(ci).cforce(d)/fstarnrm)*vstarnrm;
					cell(ci).setCVel(d,vtmp);
				}
			}
		}

		// verlet position update
		spPosVerlet();

		// reset contacts before force calculation
		resetContacts();

		// calculate forces between disks (with door closed)
		spForces(lenscales);

		// verlet velocity update
		spVelVerlet(lenscales);

		// update t
		t += dt;

		// update virial pressure
		Pvirial = 0.5*(sigmaXX + sigmaYY);

		// update kinetic energy based on com velocity
		Knew = 0.0;
		for (ci=0; ci<NCELLS; ci++)
			Knew += 0.5*(PI*pow(lenscales.at(ci),2))*sqrt(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2));

		// update if Pvirial under tol
		if (abs(Pvirial) < Ptol)
			npPMIN++;
		else
			npPMIN = 0;

		// calc check variables
		Kcheck = Knew/(NDIM*NCELLS);
		Pcheck = Pvirial/(NDIM*NCELLS);

		// check for convergence
		converged = (abs(Pvirial) < Pvirial && npPMIN > NMIN);
		converged = (converged || (abs(Pvirial) > 2*Ptol && Kcheck < Ktol));

		if (converged){
			cout << "	** FIRE has converged!" << endl;
			cout << "	** Kcheck = " << Kcheck << endl;
			cout << "	** Pcheck = " << Pcheck << endl;
			cout << "	** itr = " << itr << ", t = " << t << endl;
			cout << "	** Breaking out of FIRE protocol." << endl;
			break;
		}
	}

	// reset dt to be original value before ending function
	dt = dt0;

	// if no convergence, just stop
	if (itr == itrMax)
		cout << "	** FIRE not converged in itrMax = " << itr << " force evaluations" << endl;
}

void cellPacking2D::spForces(vector<double>& lenscales){
	// local variables
	int ci, cj, vi, d;
	double contactDistance = 0.0; 
	double centerDistance = 0.0; 
	double overlap = 0.0;
	double uv = 0.0;
	double ftmp, utmp;
	vector<double> distanceVec(NDIM,0.0);

	// reset virial stresses to 0
	sigmaXX = 0.0;
	sigmaXY = 0.0;
	sigmaYX = 0.0;
	sigmaYY = 0.0;

	// get disk-disk forces
	for (ci=0; ci<NCELLS; ci++){
		for (cj=ci+1; cj<NCELLS; cj++){
			// contact distance
			contactDistance = lenscales.at(ci) + lenscales.at(cj);

			// center-to-center distance
			centerDistance = 0.0;
			for (d=0; d<NDIM; d++){
				// vectorial quantity
				distanceVec.at(d) = cell(ci).cellDistance(cell(cj),d);

				// add to distance
				centerDistance += pow(distanceVec.at(d),2);
			}

			// check for contact
			if (contactDistance*contactDistance > centerDistance){
				// add to contact checking
				addContact(ci,cj);

				// get true distance
				centerDistance = sqrt(centerDistance);

				// overlap scale
				overlap = centerDistance/contactDistance;

				// force scale
				ftmp = (1 - overlap);

				// add to potential energy (energy should increase because particles are growing)
				utmp = 0.5*contactDistance*pow(1 - overlap,2);
				for (vi=0; vi<cell(ci).getNV(); vi++)
					cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());
				for (vi=0; vi<cell(cj).getNV(); vi++)
					cell(cj).setUInt(vi,cell(cj).uInt(vi) + utmp/cell(cj).getNV());

				// add to forces
				for (d=0; d<NDIM; d++){
					// unit vector
					uv = distanceVec.at(d)/centerDistance;

					// add to forces (MIND FORCE DIRECTION; rij points from i -> j, so need extra minus sign)
					cell(ci).setCForce(d,cell(ci).cforce(d) - ftmp*uv);
					cell(cj).setCForce(d,cell(cj).cforce(d) + ftmp*uv);

					// add to virial stresses
					if (d == 0){
						sigmaXX += ftmp*uv*distanceVec.at(0);
						sigmaXY += ftmp*uv*distanceVec.at(1);
					}
					else{
						sigmaYX += ftmp*uv*distanceVec.at(0);
						sigmaYY += ftmp*uv*distanceVec.at(1);
					}
				}
			}
		}
	}
}

void cellPacking2D::spPosVerlet(){
	// local variables
	int ci, vi, d;
	double postmp, acctmp, dpos;

	// update com position
	for (ci=0; ci<NCELLS; ci++){
		// loop over positions
		for (d=0; d<NDIM; d++){
			// calculate com acceleration
			acctmp = 0.0;
			for (vi=0; vi<cell(ci).getNV(); vi++)
				acctmp += cell(ci).vacc(vi,d);

			// update new position based on acceleration
			postmp = cell(ci).cpos(d) + dt*cell(ci).cvel(d) + 0.5*dt*dt*acctmp;

			// translate vertices based on cpos change
			for (vi=0; vi<cell(ci).getNV(); vi++)
				cell(ci).setVPos(vi,d,cell(ci).vpos(vi,d) + (postmp - cell(ci).cpos(d)));

			// update new positions
			cell(ci).setCPos(d,postmp);

			// set forces to 0
			cell(ci).setCForce(d,0.0);
		}

		// set interaction energy to 0
		for (vi=0; vi<cell(ci).getNV(); vi++)
			cell(ci).setUInt(vi,0.0);
	}
}

void cellPacking2D::spVelVerlet(vector<double>& radii){
	// local variables
	int ci, vi, d;
	double veltmp, aold, anew, diskMass;

	// update com velocity
	for (ci=0; ci<NCELLS; ci++){
		// get disk mass
		diskMass = PI*pow(radii.at(ci),2);

		// loop over velocities
		for (d=0; d<NDIM; d++){
			// get current velocity
			veltmp = cell(ci).cvel(d);

			// calculate old com acceleration
			aold = 0.0;
			for (vi=0; vi<cell(ci).getNV(); vi++)
				aold += cell(ci).vacc(vi,d);

			// get new accelation
			anew = cell(ci).cforce(d)/diskMass;

			// update velocity
			veltmp += 0.5*dt*(anew + aold);

			// set new velocity and acceleration
			cell(ci).setCVel(d,veltmp);
			for (vi=0; vi<cell(ci).getNV(); vi++)
				cell(ci).setVAcc(vi,d,anew/cell(ci).getNV());
		}
	}
}

void cellPacking2D::shrinkSP(vector<double>& lenscales){
	// local variables
	int ci, cj, vi, d;
	double Utmp, Utol;
	double contactDistance, centerDistance, overlap, ftmp, utmp;
	vector<double> distanceVec(NDIM,0.0);

	// tolerance for potential energy
	Utol = 1e-16;
	Utmp = Utol*10;

	// loop over sizes while potential energy too large
	while(Utmp > Utol){
		// reset potential to 0
		Utmp = 0.0;

		// reset forces to 0
		for (ci=0; ci<NCELLS; ci++){
			for (d=0; d<NDIM; d++)
				cell(ci).setCForce(d,0.0);
		}

		// loop over particles, get total energy, update lenscales based on forces
		// (STORE FORCES IN d = 0 ARRAY ENTRY OF CFORCE)
		for (ci=0; ci<NCELLS; ci++){
			for (cj=ci+1; cj<NCELLS; cj++){
				// contact distance
				contactDistance = lenscales.at(ci) + lenscales.at(cj);

				// center-to-center distance
				centerDistance = 0.0;
				for (d=0; d<NDIM; d++){
					// vectorial quantity
					distanceVec.at(d) = cell(ci).cellDistance(cell(cj),d);

					// add to distance
					centerDistance += pow(distanceVec.at(d),2);
				}

				// check for contact
				if (contactDistance*contactDistance > centerDistance){
					// get true distance
					centerDistance = sqrt(centerDistance);

					// overlap scale
					overlap = centerDistance/contactDistance;

					// force scale (force on radii, NOT positions)
					ftmp = -(overlap/contactDistance)*(1 - overlap);

					// add to forces
					cell(ci).setCForce(0,cell(ci).cforce(0) + ftmp);
					cell(cj).setCForce(0,cell(cj).cforce(0) - ftmp);

					// add to potential energy (energy should increase because particles are growing)
					utmp = 0.5*contactDistance*pow(1 - overlap,2);
					for (vi=0; vi<cell(ci).getNV(); vi++)
						cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());
					for (vi=0; vi<cell(cj).getNV(); vi++)
						cell(cj).setUInt(vi,cell(cj).uInt(vi) + utmp/cell(cj).getNV());

					// update monitor U
					Utmp += utmp;
				}
			}
		}

		// update sizes based on radii force (assume overdamped motion, so eq of motion is first order)
		for (ci=0; ci<NCELLS; ci++)
			lenscales.at(ci) += dt*cell(ci).cforce(0);
	}
}






// NON-EQUILIBRIUM MD FUNCTIONS


// quasistaic gel forming function
void cellPacking2D::isoExtensionQS(int plotIt, int& frameCount, double phiTarget, double dphi){
	// local variables
	int t = 0;
	int isRelaxed = 0;
	double Ftol = 1e-6;
	double Ktol = 1e-18;

	// calculate initial packing fraction
	phi = packingFraction();

	// loop over packing fraction, only decrease 
	while (phi > phiTarget){
		// check if relaxed
		if (isRelaxed){
			// reset relaxed check
			isRelaxed = 0;

			// decrement packing fraction
			phi = packingFraction();
			setPackingFraction(phi-dphi);

			// print statement
			cout << "===================================================" << endl << endl;
			cout << " 	CHANGING PACKING FRACTION FROM " << phi + dphi << " to " << phi << ": t = " << t << ", frame = " << frameCount << endl << endl;
			cout << "===================================================" << endl;
			
			cout << "	* Run data:" << endl;
			cout << "	* old phi 	= " << phi + dphi << endl;
			cout << "	* new phi 	= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << "	* nc 		= " << totalNumberOfContacts() << endl;
			cout << endl << endl;

			// increment frame count
			frameCount++;
		}
		// run relaxation until code relaxed
		else{
			fireMinimizeF(Ftol, Ktol, plotIt, frameCount);
			isRelaxed = 1;
		}

		// increment t
		t++;
	}

	if (phi < phiTarget)
		cout << "	Target phi = " << phiTarget << "reached! Ending isoExtensionQS protocol" << endl;
}






// TUMOR TISSUE FUNCTIONS

void cellPacking2D::tumorForce(int NTUMORCELLS, double forceScale, double adiposeDamping){
	// local variables
	int t, ci, frameCount, forceUnit;
	double U,K;

	// check that NT, NPRINT have been assigned values
	if (NT <= 0){
		cout << "	* ERROR: tried to run tumorNVE() function without proper NT, ending." << endl;
		exit(1);
	}
	else if (NPRINT <= 0){
		cout << "	* ERROR: tried to run tumorNVE() function without proper NPRINT, ending." << endl;
		exit(1);
	}

	// reinitialize frame counter
	frameCount = 0;

	// run NVE for allotted time
	for (t=0; t<NT; t++){
		// calculate energies
		U = totalPotentialEnergy();
		K = totalKineticEnergy();

		// print data first to get the initial condition
		if (t % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 	TUMOR NVE, t = " << t << ", frame = " << frameCount << endl << endl;
			cout << "===================================================" << endl;

			if (packingPrintObject.is_open()){
				cout << "	* Printing vetex positions to file" << endl;
				printSystemPositions(frameCount);
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing cell energy to file" << endl;
				printSystemEnergy(frameCount,U,K);
			}
			frameCount++;
			
			cout << "	* Run data:" << endl;
			cout << "	* U 		= " << U << endl;
			cout << "	* K 		= " << K << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << "	* nc 		= " << totalNumberOfContacts() << endl;
			cout << endl << endl;
		}

		// use velocity verlet to advance time

		// update positions
		for (ci=0; ci<NCELLS; ci++){
			cell(ci).verletPositionUpdate(dt);
			cell(ci).updateCPos();
		}

		// reset contacts before force calculation
		resetContacts();

		// calculate forces
		calculateForces();

		// add dragging force to tumor cell
		for (ci=0; ci<NTUMORCELLS; ci++){
			// calculate force unit
			forceUnit = cell(ci).getkint()/cell(ci).getdel();

			// add to x-direction of center-of-mass force
			cell(ci).setCForce(0, cell(ci).cforce(0) + forceScale*forceUnit );
		}

		// update velocities, with damping on adipose tissue
		for (ci=0; ci<NCELLS; ci++){
			if (ci < NTUMORCELLS)
				cell(ci).verletVelocityUpdate(dt,0.0);
			else
				cell(ci).verletVelocityUpdate(dt,adiposeDamping);
		}
	}

	// print something to conclude
	cout << "t = NT, which = " << t << ", so ending NVE run!" << endl;
}















// RELAXATION FUNCTIONS


// relieve overlap between particles
void cellPacking2D::overlapRelief(double phiT){
	// local variables
	int k = 0;
	int kmax = 1e5;
	int np = 0;
	int ci, cj, inContact;
	double alpha = 0.25;
	double rmsF, U, Ftol, Utol, energyScale, forceScale, nvScale;
	double FtolScale, UtolScale;
	double phi0;

	// vector to store radii
	vector<double> radii(NCELLS,0.0);
	double cellArea;
	int cellNV;

	// calculate radii, phi0 based on particle sizes
	phi0 = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		// get radii
		cellArea = cell(ci).area();
		cellNV = cell(ci).getNV();
		radii.at(ci) = sqrt((2.0*cellArea)/(cellNV*sin(2*PI/cellNV)));

		// update current packing fraction
		phi0 += PI*pow(radii.at(ci),2);
	}
	phi0 /= L*L;

	// scale radii to fit new packing fraction
	scaleLengths(sqrt(phiT/phi0));
	for (ci=0; ci<NCELLS; ci++)
		radii.at(ci) *= sqrt(phiT/phi0);

	// scale kint to make forces greater
	// for (ci=0; ci<NCELLS; ci++)
	// 	cell(ci).setkint(cell(ci).getkint()*cell(ci).getNV());

	// update time step
	setdt(0.1);

	// get initial forces
	diskForces(radii);

	// initialize force and potential
	rmsF = forceRMS();
	U = interactionPotentialEnergy();
	cout << "initial rmsF = " << rmsF << " and U = " << U << endl;

	// force and energy scales
	energyScale = cell(0).getkint();
	forceScale = energyScale/cell(0).getdel();
	nvScale = cell(0).getNV();

	// force and potential energy tolerance
	FtolScale = 1e-12;
	UtolScale = 1e-12;

	Ftol = NDIM*NCELLS*nvScale*forceScale*FtolScale;
	Utol = NCELLS*nvScale*energyScale*UtolScale;

	// loop over MD until force/U is minimized
	while (rmsF > Ftol && U > Utol && k < kmax){
		// do verlet update
		for (ci=0; ci<NCELLS; ci++){
			cell(ci).verletPositionUpdate(dt);
			cell(ci).updateCPos();
		}

		// reset contacts before force calculation
		resetContacts();

		// calculate forces
		diskForces(radii);

		// update velocities
		for (ci=0; ci<NCELLS; ci++)
			cell(ci).verletVelocityUpdate(dt,0.0);

		// fire step
		fireStep(np, alpha);

		// update forces
		rmsF = forceRMS();
		U = interactionPotentialEnergy();

		// update iterator
		k++;
	}

	if (k == kmax){
		cout << "	** ERROR: overlap relief could not converge in kmax = " << kmax << " iterations " << endl;
		exit(1);
	}

	// print ending information
	cout << "Starting phi as disks = " << phi0 << ", cell phi now = " << phiT << endl;
	cout << "ending rmsF/Ftol = " << rmsF/Ftol << ", U/Utol = " << U/Utol << endl;
	cout << "first radii = " << radii.at(0) << endl;
	cout << "rel radii = " << sqrt(cell(0).vrel(1,0)*cell(0).vrel(1,0) + cell(0).vrel(1,1)*cell(0).vrel(1,1)) << endl;
	cout << "meas radii = " << sqrt(pow(cell(0).vpos(1,0)-cell(0).cpos(0),2) + pow(cell(0).vpos(1,1)-cell(0).cpos(1),2)) << endl;
}



// calculate forces as if cells were disks
void cellPacking2D::diskForces(vector<double>& radii){
	// local variables
	int ci, cj, vi, d;
	double contactDistance, centerDistance, distTmp;
	double ftmp, utmp;
	vector<double> distVec(NDIM,0.0);
	double forceScale, energyScale, x, ur;

	// set force and energy scales
	energyScale = cell(0).getkint();

	// reset forces, interaction energy to 0
	for (ci=0; ci<NCELLS; ci++){
		for (vi=0; vi<cell(ci).getNV(); vi++){
			cell(ci).setUInt(vi,0.0);
			for (d=0; d<NDIM; d++)
				cell(ci).setVForce(vi,d,0.0);
		}
	}

	// loop over pairwise forces
	for (ci=0; ci<NCELLS; ci++){
		for (cj=ci+1; cj<NCELLS; cj++){
			// get contact distance
			contactDistance = radii.at(ci) + radii.at(cj) + cell(ci).getdel();

			// get actual distance
			centerDistance = 0.0;
			for (d=0; d<NDIM; d++){
				// get distance component
				distTmp = cell(ci).cellDistance(cell(cj),d);

				// add to vector
				distVec.at(d) = distTmp;

				// add to magnitude
				centerDistance += distTmp*distTmp;
			}
			centerDistance = sqrt(centerDistance);

			// add forces if overlapping
			if (centerDistance < contactDistance){
				// normalized distance x
				x = centerDistance/contactDistance;

				// force scale
				forceScale = energyScale/contactDistance;

				// add to interaction potential
				utmp = 0.5*energyScale*pow(1 - x,2);
				for (vi=0; vi<cell(ci).getNV(); vi++)
					cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());
				for (vi=0; vi<cell(cj).getNV(); vi++)
					cell(cj).setUInt(vi,cell(cj).uInt(vi) + utmp/cell(cj).getNV());

				// add to center-of-mass forces
				for (d=0; d<NDIM; d++){
					// unit vector pointing in OPPOSITE direction of force
					ur = distVec.at(d)/centerDistance;

					// component of force
					ftmp = -forceScale*(1 - x)*ur;

					// add to forces
					cell(ci).setCForce(d,cell(ci).cforce(d) + ftmp);
					cell(cj).setCForce(d,cell(cj).cforce(d) - ftmp);

					// add to virial stresses
					if (d == 0){
						sigmaXX += ftmp*distVec.at(0);
						sigmaXY += ftmp*distVec.at(1);
					}
					else{
						sigmaYX += ftmp*distVec.at(0);
						sigmaYY += ftmp*distVec.at(1);
					}
				}
			}
		}
	}
}

// ramp shape changes
void cellPacking2D::shapeRamp(double fixedPhi, double calATarget, double dCalA, double kbTarget, double dkb){
	// local variables
	int k, ci, delSgn, frameCount;
	double da, currCalA, decScale, Ptol, Ktol;
	double calACheck, tol;
	double kb = cell(0).getkb();

	// set check and tolerance
	tol = 1e-8;
	calACheck = 10*tol;

	// set tolerance for relaxation
	Ptol = 1e-8;
	Ktol = 1e-24;

	// get current attraction (assuming all attraction the same)
	currCalA = cell(0).calA0();

	// set decrease scale to be 1
	decScale = 1.0;

	// check whether to increase or decrease
	if (currCalA < calATarget)
		delSgn = 1;
	else
		delSgn = -1;

	// loop until attraction is the correct value
	k=0;
	frameCount=0;
	while (calACheck > tol && k < NT){
		// output to console
		cout << "===================================================" << endl << endl << endl;
		cout << " 	Ramping cal A to calATarget = " << calATarget << endl << endl;
		cout << "===================================================" << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* current CalA 	= " << currCalA << endl;
		cout << "	* current kb 	= " << kb << endl;
		cout << endl << endl;

		// decide whether attraction is too large or too small
		if (delSgn*currCalA < delSgn*calATarget)
			da = delSgn*dCalA;
		else{
			decScale *= 0.9;
			da = -delSgn*decScale*dCalA;
		}

		// update current attraction parameter
		currCalA += da;

		// update attraction in all of the cells
		for (ci=0; ci<NCELLS; ci++)
			cell(ci).setAsphericityConstA(currCalA);

		// relax potential energy
		// potentialRelaxFire(Ktol, Utol, 1, frameCount);
		fireMinimizeP(Ptol, Ktol);

		// update check
		calACheck = abs(currCalA - calATarget);

		// update packing fraction (keep fixed)
		phi = packingFraction();
		setPackingFraction(fixedPhi);

		// increase bending energy (ASSUME DISTANCE TO kbTarget IS NOT PART OF WHILE LOOP)
		if (kb < kbTarget){
			kb += dkb;

			for (ci=0; ci<NCELLS; ci++)
				cell(ci).setkb(kb);
		}

		// update iterator
		k++;
	}

	if (k == NT){
		cout << "	ERROR: particle shapes could not relax to desired calA in allotted time, ending." << endl;
		exit(1);
	}
}



/************************

	Printers

*************************/


// print positions to file (with frame)
void cellPacking2D::printSystemPositions(int frame){
	// local variables
	int w1 = 12;
	int w2 = 6;

	// check to see if file is open
	if (!packingPrintObject.is_open()) {
		cout << "	ERROR: packingPrintObject is not open in printSystemPositions(), ending." << endl;
		exit(1);
	}

	// print header if frame == 0
	if (frame == 0){
		packingPrintObject << setw(w1) << left << "START" << " " << endl;
		packingPrintObject << setw(w1) << left << "NUMCL" << setw(w2) << right << NCELLS << endl;
		packingPrintObject << setw(w1) << left << "NUMFR" << setw(w2) << right << nframes() << endl;
		packingPrintObject << setw(w1) << left << "BOXSZ" << setw(w2) << right << L << endl;
	}

	// print cell information for 0, to get new frame header
	cell(0).printVertexPositions(packingPrintObject,0,frame);

	// print info for rest of the cells
	for (int ci=1; ci<NCELLS; ci++)
		cell(ci).printVertexPositions(packingPrintObject,ci);

	// print terminal command if on final frame
	if (frame == NPRINT*(nframes()-1))
		packingPrintObject << setw(12) << left << "STERM" << " " << endl;
}

// print energies to file (with frame)
void cellPacking2D::printSystemEnergy(int frame, double Pval, double Kval){
	// local variables
	int ci;

	// check to see if file is open
	if (!energyPrintObject.is_open()) {
		cout << "	ERROR: energyPrintObject is not open in printSystemEnergy(), ending." << endl;
		exit(1);
	}

	// loop over particles, print cell energy
	energyPrintObject << setw(6) << right << frame;
	energyPrintObject << setw(30) << setprecision(16) << right << Kval;
	energyPrintObject << setw(30) << setprecision(16) << right << interactionPotentialEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << totalPotentialEnergy();

	// ADD ON STUFF FOR DEBUGGING JAMMING FINDER
	double calA0, meanCalA;

	calA0 = (pow(cell(0).getl0(),2)*pow(cell(0).getNV(),2))/(4*PI*cell(0).geta0());
	meanCalA = cell(0).asphericity();

	energyPrintObject << setw(30) << setprecision(16) << right << Pval;
	energyPrintObject << setw(30) << setprecision(16) << right << calA0;
	energyPrintObject << setw(30) << setprecision(16) << right << meanCalA;
	energyPrintObject << endl;

	// NOTE: HOW TO ADD ENERGIES OF INDIVIDUAL CELLS?
}


void cellPacking2D::printSystemPositions(){
	// local variables
	int w1 = 12;
	int w2 = 6;
	int w3 = 30;

	// check to see if file is open
	if (!packingPrintObject.is_open()) {
		cout << "	ERROR: packingPrintObject is not open in printSystemPositions(), ending." << endl;
		exit(1);
	}

	// print information starting information
	packingPrintObject << setw(w1) << left << "NEWFR" << " " << endl;
	packingPrintObject << setw(w1) << left << "NUMCL" << setw(w2) << right << NCELLS << endl;

	// print hopper information
	packingPrintObject << setw(w1) << left << "BOXSZ";
	packingPrintObject << setw(w2) << right << L;
	packingPrintObject << endl;

	// print stress information
	packingPrintObject << setw(w1) << left << "VRIAL";
	packingPrintObject << setw(w3) << right << sigmaXX;
	packingPrintObject << setw(w3) << right << sigmaXY;
	packingPrintObject << setw(w3) << right << sigmaYX;
	packingPrintObject << setw(w3) << right << sigmaYY;
	packingPrintObject << endl;

	// print info for rest of the cells
	for (int ci=0; ci<NCELLS; ci++)
		cell(ci).printVertexPositions(packingPrintObject,ci);

	// print end frame
	packingPrintObject << setw(w1) << left << "ENDFR" << " " << endl;
}

void cellPacking2D::printSystemEnergy(){

	// check to see if file is open
	if (!energyPrintObject.is_open()) {
		cout << "	ERROR: energyPrintObject is not open in printSystemEnergy(), ending." << endl;
		exit(1);
	}

	// loop over particles, print cell energy
	energyPrintObject << setw(30) << setprecision(16) << right << interactionPotentialEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << totalPotentialEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << totalKineticEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaXX;
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaXY;
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaYX;
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaYY;
	energyPrintObject << setw(30) << setprecision(16) << right << packingFraction();
	energyPrintObject << endl;

	// NOTE: HOW TO ADD ENERGIES OF INDIVIDUAL CELLS?
}

void cellPacking2D::printSystemContacts(){
	// check to see if file is open
	if (!statPrintObject.is_open()) {
		cout << "	ERROR: statPrintObject is not open in printSystemContacts(), ending." << endl;
		exit(1);
	}

	// loop over contacts
	for (int ci=0; ci<NCELLS; ci++){
		for (int cj=ci+1; cj<NCELLS; cj++)
			statPrintObject << setw(6) << contacts(ci,cj);
	}
	statPrintObject << endl;
}


// print stats to a file
void cellPacking2D::printSystemStats(){
	// local variables
	int p = 16;
	int w = 6;
	int ci, cj;

	// print information
	statPrintObject << NCELLS << endl;
	statPrintObject << L << endl;
	statPrintObject << setprecision(p) << phi << endl;

	// print contact matrix
	for (ci=0; ci<NCELLS; ci++){
		for (cj=ci+1; cj<NCELLS; cj++)
			statPrintObject << setw(w) << contacts(ci,cj);
	}
	statPrintObject << endl;
}



