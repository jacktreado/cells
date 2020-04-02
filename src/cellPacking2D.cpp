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
	T 			= -1.0;
	phi 		= -1.0;
	sigmaXX 	= 0.0;
	sigmaXY 	= 0.0;
	sigmaYX 	= 0.0;
	sigmaYY 	= 0.0;

	// set box lengths to 1.0
	L.resize(NDIM);
	for (int d=0; d<NDIM; d++)
		L.at(d) = 1.0;

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
	int i, d, NC;

	// set initial seed
	seed = s;

	// set random number generator
	srand(seed);

	// first use starting variabls
	defaultvars();

	// set member variables based on inputs
	NCELLS 		= ncells;
	NT 			= nt;
	NPRINT 		= nprint;

	// set box lengths to be square
	for (d=0; d<NDIM; d++)
		L.at(d) = l;

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
	else if (L.at(0) <= 0.0 || L.at(1) <= 0.0){
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
	int i, ci, d, NC;

	// set initial seed
	seed = s;

	// first use starting variabls
	defaultvars();

	// set member variables based on inputs
	NCELLS 		= ncells;
	
	// set box lengths to be square
	for (d=0; d<NDIM; d++)
		L.at(d) = 5.0*NCELLS;

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

		// set box length for each cell ( WITH PBCs = 1 )
		for (d=0; d<NDIM; d++){
			cell(ci).setL(d,L.at(d));
			cell(ci).setpbc(d,1);
		}

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
	double ltmp, l0tmp, l0, a0, posTmp;

	// read in simulation information
	inputFileObject >> NCELLS;
	inputFileObject >> ltmp;
	inputFileObject >> l0tmp;

	cout << "NCELLS = " << NCELLS << endl;
	cout << "L = " << ltmp << endl;

	// set box lengths to be square, reset lengths to be in units of l0
	for (d=0; d<NDIM; d++){
		L.at(d) = ltmp;
		L.at(d) /= l0tmp;
	}

	// test for error in inputs
	if (NCELLS <= 0){
		cout << "	ERROR: in overloaded operator, NCELLS <= 0 so cannot initialize arrays, ending code here." << endl;
		exit(1);
	}
	else if (L.at(0) <= 0.0 || L.at(1) <= 0.0){
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

		// initialize cell objects ( set PBCs = 1 )
		for (d=0; d<NDIM; d++){
			cell(ci).setL(d,L.at(d));
			cell(ci).setpbc(d,1);
		}

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
		l0 = (1.0/nv)*sqrt(4*PI*cell(ci).polygonArea()*asphericity);
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
	T 		= onTheRight.T;
	phi 	= onTheRight.phi;

	// box length
	for (int d=0; d<NDIM; d++)
		L.at(d) = onTheRight.getL(d);

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
		cell(ci).setNV(onTheRight.cell(ci).getNV());
		cell(ci).initializeVertices();
		cell(ci).initializeCell();
		cell(ci) = onTheRight.cell(ci);

		// copy elements from contact matrix
		for (cj=ci+1; cj<NCELLS; cj++){
			if (onTheRight.contacts(ci,cj))
				addContact(ci,cj);
		}
	}
}


// save current state into saveObject
void cellPacking2D::saveState(cellPacking2D& saveObject){
	// local variables
	int ci, cj;

	// test if save object is properly initialized
	if (saveObject.phi < 0)
		// use overloaded assignment to store current state into object
		saveObject = *this;
	else{
		// if saveObject has been initialized, copy data, cell objects and contact network
		saveObject.dt = dt;
		saveObject.dt0 = dt0;
		saveObject.phi = phi;

		saveObject.sigmaXX = sigmaXX;
		saveObject.sigmaXY = sigmaXY;
		saveObject.sigmaYX = sigmaYX;
		saveObject.sigmaYY = sigmaYY;

		for (ci=0; ci<NCELLS; ci++){
			// copy cell objects (using overloaded operator in deformableParticle2D class)
			saveObject.cell(ci) = cell(ci);

			// copy elements from contact matrix
			for (cj=ci+1; cj<NCELLS; cj++){
				if (contacts(ci,cj))
					saveObject.addContact(ci,cj);
			}
		}
	}
}

// load saved state from saveObject
void cellPacking2D::loadState(cellPacking2D& saveObject){
	// local variables
	int ci, cj;

	// test that saveObject has been initialized
	if (saveObject.phi < 0){
		cout << "	** ERROR: trying to load from saveObject that has no saved data, ending." << endl;
		exit(1);
	}

	// load saved packing fraction, dt
	dt = saveObject.dt;
	dt0 = saveObject.dt0;
	phi = saveObject.phi;

	sigmaXX = saveObject.sigmaXX;
	sigmaXY = saveObject.sigmaXY;
	sigmaYX = saveObject.sigmaYX;
	sigmaYY = saveObject.sigmaYY;

	// load cell and contact data from saveObject
	for (ci=0; ci<NCELLS; ci++){
		// copy cell objects (using overloaded operator in deformableParticle2D class)
		cell(ci) = saveObject.cell(ci);

		// copy elements from contact matrix
		for (cj=ci+1; cj<NCELLS; cj++){
			if (saveObject.contacts(ci,cj))
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
	int i, d;
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
		// set box length for each cell ( set PBCs = 1 )
		for (d=0; d<NDIM; d++){
			cell(i).setL(d,L.at(d));
			cell(i).setpbc(d,1);
		}

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
	int i,nv,nvLarge,smallIndex,d;
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
		// set box length for each cell ( set PBCs = 1 )
		for (d=0; d<NDIM; d++){
			cell(i).setL(d,L.at(d));
			cell(i).setpbc(d,1);
		}

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
	int i,nv,nvLarge,smallIndex,d;
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
		// set box length for each cell ( set PBC = 1 )
		for (d=0; d<NDIM; d++){
			cell(i).setL(d,L.at(d));
			cell(i).setpbc(d,1);
		}

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
		cell(i).setdel(del);
		cell(i).seta(a);
	}
}


// set initial positions to populate a square lattice
void cellPacking2D::squareLattice(){
	// local variables
	int ci, xIndex, yIndex;
	int gridPoints = round(1.5*NCELLS);

	double xpos, ypos;
	double buffer = 0.05*L.at(0);
	double spacing = (L.at(0) - 2.0*buffer)/(gridPoints-1);

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
		cout << " and calA0 = " << cell(ci).calA0() << endl;
	}

	// calculate packing fraction
	phi = packingFraction();

	// print statement
	cout << "Particles initialized on square lattice with initial packing fraction phi = " << phi << endl;
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
	vscale = sqrt(T/(cell(0).polygonArea()*ek));
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
	vscale = sqrt(tmp0/(cell(ci).polygonArea()*ek));
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
	double atot = 0.0;
	double val = 0.0;

	// loop over full area of cells (not just polygon area)
	for (ci=0; ci<NCELLS; ci++)
		atot += cell(ci).area();

	// divide by box area
	val = atot/(L.at(0)*L.at(1));

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

// set vertex DPM time step mag
void cellPacking2D::vertexDPMTimeScale(double timeStepMag){
	// local variables
	int ci;
	double a0mean, l0mean, kamean;
	double t0;

	// initialize means
	a0mean = 0.0;
	l0mean = 0.0;
	kamean = 0.0;

	// loop over cells, get mean a0, l0, r0
	for(ci=0; ci<NCELLS; ci++){
		a0mean += cell(ci).geta0()/NCELLS;
		l0mean += cell(ci).getl0()/NCELLS;
		kamean += cell(ci).getka()/NCELLS;
	}

	// get fundamental time step
	t0 = sqrt((2.0*l0mean)/(kamean*a0mean));

	// set dt and dt0
	dt = timeStepMag*t0;
	dt0 = dt;
}


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

	// calculate val to scale lengths with
	scaleFactor	= pow(val/phi,1.0/NDIM);

	// scale all lengths by scale factor
	scaleLengths(scaleFactor);
}


// scale all length scales in the system by scaleFactor
void cellPacking2D::scaleLengths(double scaleFactor){
	// local variables
	int i;

	// loop over cells, use scale to change lengths
	for (i=0; i<NCELLS; i++)
		cell(i).scale(scaleFactor);
}

void cellPacking2D::setAsphericity(double val){
	// local variables
	int ci, nvtmp;
	double calAMin = 0.0;

	// set all cells to specified asphericity values (in units of calAMin)
	for (ci=0; ci<NCELLS; ci++){
		// number of vertices for cell ci
		nvtmp = cell(ci).getNV();

		// determine calAmin for given cell
		calAMin = nvtmp*tan(PI/nvtmp)/PI;

		// set asphericity in units of calAMin
		cell(ci).setAsphericityConstA(val*calAMin);
	}
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

			// compute virial stresses
			sigmaXX += fij.at(0)*rij.at(0);
			sigmaXY += fij.at(1)*rij.at(0);
			sigmaYX += fij.at(0)*rij.at(1);
			sigmaYY += fij.at(1)*rij.at(1);
		}

		// forces on vertices due to shape
		cell(ci).balancedShapeForces();
	}
}

void cellPacking2D::gelationForces(){
	// local variables
	int ci,cj,d,dd,inContact,numACtmp;
	double aij;

	// vector to store pairwise forces
	vector<double> fij(NDIM,0.0);
	vector<double> rij(NDIM,0.0);

	// vector to store number of attractive contacts per cell
	vector<int> nac(NCELLS,0);

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

	// get number of attractive contacts
	// loop over cell pairs
	for (ci=0; ci<NCELLS; ci++){
		for (cj=0; cj<NCELLS; cj++){

			// find number of pairwise attractive contacts between vertices
			// on ci and cj
			numACtmp = cell(ci).pwAttractiveContacts(cell(cj));
			nac.at(ci) += numACtmp;

		}
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

			// if attractive contact engaged
			if (nac.at(ci) > 0 && nac.at(cj) > 0){
				// get effective attraction scale (max is 0.5*a)
				aij = 0.25*(cell(ci).geta()/nac.at(ci) + cell(cj).geta()/nac.at(cj));

				// calculate forces
				inContact = cell(ci).vertexForce(cell(cj),fij,rij,aij);
			}
			// else, use normal force routine
			else
				inContact = cell(ci).vertexForce(cell(cj),fij,rij);
			
			if (inContact == 1)
				addContact(ci,cj);

			// compute virial stresses
			sigmaXX += fij.at(0)*rij.at(0);
			sigmaXY += fij.at(1)*rij.at(0);
			sigmaYX += fij.at(0)*rij.at(1);
			sigmaYY += fij.at(1)*rij.at(1);
		}

		// forces on vertices due to shape
		cell(ci).balancedShapeForces();
	}

	// normalize virial stresses by the area
	// sigmaXX /= L.at(0)*L.at(1);
	// sigmaXY /= L.at(0)*L.at(1);
	// sigmaYX /= L.at(0)*L.at(1);
	// sigmaYY /= L.at(0)*L.at(1);
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
		Pvirial = 0.5*(sigmaXX + sigmaYY)/(L.at(0)*L.at(1));
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


// FIRE 2.0 pressure minimzation with backstepping if P < 0
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
	Pvirial = 0.5*(sigmaXX + sigmaYY)/(L.at(0)*L.at(1));

	// initialize energy and force tracking, pressure
	Knew = totalKineticEnergy();

	// scale P and K for convergence checking
	Pcheck = Pvirial/NCELLS;
	Kcheck = Knew/NCELLS;

	// reset velocities to 0
	for (ci=0; ci<NCELLS; ci++){
		for (vi=0; vi<cell(ci).getNV(); vi++){
			for (d=0; d<NDIM; d++)
				cell(ci).setVVel(vi,d,0.0);
		}
	}

	// iterate through MD time until system converged
	kmax = 5e6;
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
		Pvirial = 0.5*(sigmaXX + sigmaYY)/(L.at(0)*L.at(1));

		// scale P and K for convergence checking
		Pcheck = Pvirial;
		Kcheck = Knew/NCELLS;

		// update if Pvirial under tol
		if (abs(Pcheck) < Ptol)
			npPMIN++;
		else
			npPMIN = 0;

		// check for convergence
		converged = (abs(Pcheck) < Ptol && npPMIN > NMIN && Kcheck < 100*Ktol);
		converged = (converged || (abs(Pcheck) > Ptol && Kcheck < Ktol));

		if (converged){
			cout << "	** FIRE has converged!" << endl;
			cout << "	** Kcheck = " << Kcheck << endl;
			cout << "	** Pcheck = " << Pcheck << endl;
			cout << "	** k = " << k << ", t = " << t << endl;
			cout << "	** Breaking out of FIRE protocol." << endl;

			// print minimized config, energy and contact network
			if (packingPrintObject.is_open()){
				cout << "	* Printing vetex positions to file" << endl;
				printSystemPositions();
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing cell energy to file" << endl;
				printSystemEnergy(k);
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
	if (k == kmax){
		cout << "	** ERROR: FIRE not converged in kmax = " << kmax << " force evaluations, ending code" << endl;
		exit(1);
	}
}


// FIRE 2.0 force minimization with backstepping
void cellPacking2D::fireMinimizeF(double Ftol, double Ktol, double& Fcheck, double& Kcheck){
	// HARD CODE IN FIRE PARAMETERS
	const double alpha0 	= 0.1;
	const double finc 		= 1.05;
	const double fdec 		= 0.5;
	const double falpha 	= 0.99;
	const double dtmax 		= 5*dt0;
	const double dtmin 		= 0.01*dt0;
	const int NMIN 			= 20;
	const int NNEGMAX 		= 2000;
	const int NDELAY 		= 1000;
	int npPos				= 0;
	int npNeg 				= 0;
	int npPMIN				= 0;
	double alpha 			= alpha0;
	double alphat 			= alpha;
	double t 				= 0.0;
	double energyScale 		= cell(0).getka()*pow(cell(0).geta0()/cell(0).getl0(),2.0);
	double forceScale		= cell(0).getka()*cell(0).geta0()/cell(0).getl0();
	double P 				= 0;

	// local variables
	int ci,vi,d,k,kmax;
	double vstarnrm,fstarnrm,vtmp,ftmp;
	double K, F, Pcheck;

	// variable to test for potential energy minimization
	bool converged = false;

	// reset time step
	dt = dt0;

	// initialize forces
	resetContacts();
	calculateForces();

	// norm of total force vector, kinetic energy
	F = forceRMS();
	K = totalKineticEnergy();
	Pcheck = 0.5*(sigmaXX + sigmaYY)/(L.at(0)*L.at(1));
	Pcheck /= forceScale;

	// scale P and K for convergence checking
	Fcheck = F/forceScale;
	Kcheck = K/(NCELLS*energyScale);

	// reset velocities to 0
	for (ci=0; ci<NCELLS; ci++){
		for (vi=0; vi<cell(ci).getNV(); vi++){
			for (d=0; d<NDIM; d++)
				cell(ci).setVVel(vi,d,0.0);
		}
	}

	// iterate through MD time until system converged
	kmax = 1e6;
	for (k=0; k<kmax; k++){

		// output some information to console
		if (k % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 	FIRE MINIMIZATION, k = " << k << endl << endl;
			cout << "===================================================" << endl;			
			cout << "	* Run data:" << endl;
			cout << "	* Kcheck 	= " << Kcheck << endl;
			cout << "	* Fcheck 	= " << Fcheck << endl;
			cout << "	* Pcheck 	= " << Pcheck << endl;
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
		F = forceRMS();
		K = totalKineticEnergy();
		Pcheck = 0.5*(sigmaXX + sigmaYY)/(L.at(0)*L.at(1));
		Pcheck /= forceScale;

		// scale P and K for convergence checking
		Fcheck = F/forceScale;
		Kcheck = K/(NCELLS*energyScale);

		// update if Fcheck under tol
		if (abs(Fcheck) < Ftol)
			npPMIN++;
		else
			npPMIN = 0;

		// check that P is not crazy
		if (abs(P) > 800){
			cout << "	ERROR: P = " << P << ", ending." << endl;
			cout << "	** Kcheck = " << Kcheck << endl;
			cout << "	** Fcheck = " << Fcheck << endl;
			cout << "	** Pcheck = " << Pcheck << endl;
			cout << "	** k = " << k << ", t = " << t << endl;

			// print minimized config, energy and contact network
			if (packingPrintObject.is_open()){
				cout << "	* Printing vetex positions to file" << endl;
				printSystemPositions();
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing cell energy to file" << endl;
				printSystemEnergy(k);
			}

			exit(1);
		}

		// check for convergence
		converged = (abs(Fcheck) < Ftol && npPMIN > NMIN && Kcheck < Ktol);

		if (converged){
			cout << "	** FIRE has converged!" << endl;
			cout << "	** F = " << F << endl;
			cout << "	** K = " << K << endl;
			cout << "	** Kcheck = " << Kcheck << endl;
			cout << "	** Fcheck = " << Fcheck << endl;
			cout << "	** Pcheck = " << Pcheck << endl;
			cout << "	** k = " << k << ", t = " << t << endl;
			cout << "	** Breaking out of FIRE protocol." << endl;

			// print minimized config, energy and contact network
			if (packingPrintObject.is_open()){
				cout << "	* Printing vetex positions to file" << endl;
				printSystemPositions();
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing cell energy to file" << endl;
				printSystemEnergy(k);
			}

			break;
		}
	}

	// reset dt to be original value before ending function
	dt = dt0;

	// if no convergence, just stop
	if (k == kmax){
		cout << "	** ERROR: FIRE not converged in kmax = " << kmax << " force evaluations, ending code" << endl;
		exit(1);
	}
}

void cellPacking2D::fireMinimizeGel(double Ptol, double Ktol){
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
	Pvirial = 0.5*(sigmaXX + sigmaYY)/(L.at(0)*L.at(1));

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
		gelationForces();

		// update velocities
		for (ci=0; ci<NCELLS; ci++)
			cell(ci).verletVelocityUpdate(dt);

		// update t
		t += dt;

		// track energy and forces
		Knew = totalKineticEnergy();
		Pvirial = 0.5*(sigmaXX + sigmaYY)/(L.at(0)*L.at(1));

		// scale P and K for convergence checking
		Pcheck = Pvirial/(energyScale*NCELLS);
		Kcheck = Knew/(energyScale*NCELLS);

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

			// print minimized config, energy and contact network
			if (packingPrintObject.is_open()){
				cout << "	* Printing vetex positions to file" << endl;
				printSystemPositions();
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing cell energy to file" << endl;
				printSystemEnergy(k);
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




// GELATION FUNCTIONS


// function to test isotropic extension protocol
void cellPacking2D::twoParticleContact(int NV){
	// local variables
	int ci, vi, d, nvtmp;
	double rtmp, l0tmp, a0tmp, p0tmp;

	// box length from 3 particle diameters (each has unit radius)
	L.at(0) = 6.0;
	L.at(1) = 6.0;

	// number of vertices on both particles
	nvtmp = NV;

	// generate cells and sizes
	for (ci=0; ci<NCELLS; ci++){
		// boundary information ( SET PBCS TO 1 for plant cells )
		for (d=0; d<NDIM; d++){
			cell(ci).setL(d,L.at(d));
			cell(ci).setpbc(d,1);
		}

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
	cell(0).setCPos(0,0.5*L.at(0) - 1.0 - 0.5*cell(0).getl0());
	cell(0).setCPos(1,0.5*L.at(0));

	cell(1).setCPos(0,0.5*L.at(1) + 1.0 + 0.5*cell(0).getl0());
	cell(1).setCPos(1,0.5*L.at(1));

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

	cout << "		-- Initial packing fraction = " << phi << endl;
}


// initialize plant cell particles as disks at input packing fraction
void cellPacking2D::initializeGel(int NV, double phiDisk, double sizeDispersion, double delval, double ka){
	// local variables
	int ci, vi, d, nvtmp;
	int lx, ly;
	double calA;
	double xpos, ypos, dx, dy;
	double xmin, xmax, ymin, ymax;
	double r1, r2, g1, radsum;
	double rtmp, a0tmp, l0tmp, p0tmp, calA0tmp;
	double a0mean, l0mean, r0mean;
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
	for (d=0; d<NDIM; d++)
		L.at(d) = sqrt(PI*radsum/phiDisk);

	// set phi to input
	phi = phiDisk;

	// reseed rng
	srand48(56835698*seed);

	// calculate mean a0, l0, r0 so correct time scale can be calculated
	a0mean = 0.0;
	l0mean = 0.0;
	r0mean = 0.0;

	// initialize cell information
	cout << "		-- Ininitializing plant cell objects" << endl;
	for (ci=0; ci<NCELLS; ci++){
		// boundary information ( SET PBCS TO 1 for plant cells )
		for (d=0; d<NDIM; d++){
			cell(ci).setL(d,L.at(d));
			cell(ci).setpbc(d,1);
		}

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

		// calculate means
		a0mean += a0tmp;
		l0mean += l0tmp;
		r0mean += rtmp;

		// set preferred area and length 
		cell(ci).seta0(a0tmp);
		cell(ci).setl0(l0tmp);
		cell(ci).setdel(delval);
	}

	// take means
	a0mean /= NCELLS;
	l0mean /= NCELLS;
	r0mean /= NCELLS;

	// initialize particle positions
	cout << "		-- Ininitializing cell positions" << endl;

	// set min and max values of positions
	xmin = 0;
	xmax = L.at(0);
	ymin = 0;
	ymax = L.at(1);

	// set lattice spacing


	for (ci=0; ci<NCELLS; ci++){

		// map onto random x and y position
		xpos = (xmax-xmin)*drand48() + xmin;
		ypos = (ymax-ymin)*drand48() + ymin;

		// set as initial position of com
		cell(ci).setCPos(0,xpos);
		cell(ci).setCPos(1,ypos);

		// initialize vertices as a regular polygon
		cell(ci).regularPolygon();

		// perturb vertex positions a little bit
		cell(ci).vertexPerturbation(0.1);
	}

	// initial time scales
	cout << "		-- Ininitializing time scale" << endl;
	dt = 0.01*sqrt((2.0*l0mean)/(ka*a0mean));
	dt0 = dt;

	// use FIRE in PBC box to relax overlaps
	cout << "		-- Using FIRE to relax overlaps..." << endl;
	fireMinimizeSP(lenscales);
}

// set forces values for plant cells
void cellPacking2D::gelForceVals(double calA0, double kl, double ka, double gam, double kb, double kint, double del, double a){
	// local variables
	int ci;

	// set asphericity for all particles in sim
	setAsphericity(calA0);

	// loop over cells, add values
	for (ci=0; ci<NCELLS; ci++){
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
void cellPacking2D::qsIsoCompression(double phiTarget, double deltaPhi, double Ftol, double Ktol){
	// local variables
	double phi0, phiNew, dphi, Fcheck, Kcheck;
	int NSTEPS, k;

	// get initial packing fraction
	phi = packingFraction();
	phi0 = phi;

	// determine number of steps to target
	NSTEPS = floor((phiTarget - phi0)/deltaPhi);
	if (NSTEPS == 0)
		NSTEPS = 1;

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
		cout << " 	quasistatic isotropic compression with NSTEPS = " << NSTEPS << " and dphi = " << dphi << endl << endl;
		cout << "===================================================" << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* NSTEPS 		= " << NSTEPS << endl;		
		cout << "	* dphi 			= " << dphi << endl << endl;
		cout << "	AFTER LAST MINIMIZATION:" << endl;
		cout << "	* phi 			= " << phi << endl;
		cout << "	* Fcheck 		= " << Fcheck << endl;
		cout << "	* Kcheck 		= " << Kcheck << endl;
		cout << endl << endl;

		// increase packing fraction to new phi value
		phiNew = phi0 + k*dphi;
		setPackingFraction(phiNew);

		// calculate phi after minimization
		phi = packingFraction();

		// relax shapes (energies calculated in relax function)
		fireMinimizeF(Ftol, Ktol, Fcheck, Kcheck);

		// calculate phi after minimization
		phi = packingFraction();
	}
}


// compress isostatically to jamming
void cellPacking2D::findJamming(double dphi0, double Ktol, double Ftol, double Ptol){
	// local variables
	double phiNew, dphi, dphiH, Ptest, Ktest, Ftest;
	double forceScale = cell(0).getka()*cell(0).geta0()/cell(0).getl0();
	int NSTEPS, k, kmax, kr, nc, nr;
	cellPacking2D savedState;

	double phiBefore = phi;

	// calculate phi before initial minimization
	phi = packingFraction();

	// relax shapes (energies calculated in relax function)
	cout << "	** IN findJamming, performing initial relaxation" << endl;
	fireMinimizeF(Ftol, Ktol, Ftest, Ktest);

	// calculate Ptest for testing
	Ptest = 0.5*(sigmaXX + sigmaYY)/(L.at(0)*L.at(1));
	Ptest /= forceScale;

	// update number of contacts
	nc = totalNumberOfContacts();

	// get initial packing fraction
	phi = packingFraction();

	// save last state before packing fraction is changed
	saveState(savedState);

	// iterator
	k = 0;
	kmax = 1e5;

	// jamming variables
	bool jammed, overcompressed, undercompressed;
	double phiH, phiL;

	// initialize to unjammed
	jammed = false;

	// phiJ bounds
	phiH = -1;
	phiL = -1;

	// loop until phi is the correct value
	while (!jammed && k < kmax){
		// update iterator
		k++;

		// relax shapes (energies/forces calculated during FIRE minimization)
		fireMinimizeF(Ftol, Ktol, Ftest, Ktest);

		// calculate Ptest for comparison
		Ptest = 0.5*(sigmaXX + sigmaYY)/(L.at(0)*L.at(1));
		Ptest /= forceScale;

		// remove rattlers
		kr = 0;
		nr = removeRattlers(kr);

		// update number of contacts
		nc = totalNumberOfContacts();

		// boolean checks (NOTE: Ktest < Ktol NEEDS TO BE INCLUDED IN FIRE MINIMIZATION, OTHERWISE INCLUDE IT HERE)
		undercompressed = ((Ptest < 2.0*Ptol && phiH < 0) || (Ptest < Ptol && phiH > 0));
		overcompressed = (Ptest > 2.0*Ptol && nc > 0);
		jammed = (Ptest < 2.0*Ptol && Ptest > Ptol && nc > 0 && phiH > 0);

		// output to console
		cout << "===================================================" << endl << endl << endl;
		cout << " 	quasistatic isotropic compression to jamming " << endl << endl;
		cout << "===================================================" << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* dphi 			= " << dphi << endl;
		cout << "	* phi 			= " << phi << endl;
		cout << "	* phiH 			= " << phiH << endl;
		cout << "	* phiL 			= " << phiL << endl;
		cout << "	* Ftest 		= " << Ftest << endl;
		cout << "	* Ktest 		= " << Ktest << endl;
		cout << "	* Ptest 		= " << Ptest << endl;
		cout << "	* # of contacts = " << nc << endl;
		cout << "	* # of rattlers = " << nr << endl;
		cout << "	* undercompressed = " << undercompressed << endl;
		cout << "	* overcompressed = " << overcompressed << endl;
		cout << "	* jammed = " << jammed << endl;
		cout << endl << endl;

		// update packing fraction based on jamming check
		dphi = 0.0;
		if (phiL < 0){
			// if still undercompressed, then grow until overjammed found
			if (undercompressed)
				dphi = dphi0;
			// if first overcompressed, return to pre-overcompression state, to midpoint between phi and phiH
			else if (overcompressed){
				phiH = phi;
				loadState(savedState);
				phiL = phi;
				dphi = (0.5*(phiH + phiL)) - phi;
				cout << "	-- -- overcompressed for first time, resetting to phi0 = " << phi << ", phiH = " << phiH << ", compressing by dphi = " << dphi << endl;
			}
		}
		else{
			// if found undercompressed state, go to state between undercompressed and last overcompressed states (from saved state)
			if (undercompressed){
				phiL = phi;
				loadState(savedState);
				dphi = (0.5*(phiH + phiL)) - phi;
				cout << "	-- -- now undercompressed, resetting to phi0 = " << phi << ", phiH = " << phiH << ", compressing by dphi = " << dphi << endl;

			}
			else if (overcompressed){
				phiH = phi;
				loadState(savedState);
				dphi = (0.5*(phiH + phiL)) - phi;
				cout << "	-- -- overcompressed (phiL > 0), rresetting to phi0 = " << phi << ", phiH = " << phiH << ", compressing by dphi = " << dphi << endl;
			}
			else if (jammed){
				cout << "	** At k = 0, jamming found!" << endl;
				cout << "	** phiJ = " << phi << endl;
				cout << "	** F = " << Ftest << endl;
				cout << "	** P = " << Ptest << endl;
				cout << "	** K = " << Ktest << endl;
				cout << "	** nc = " << nc << endl;
				cout << " WRITING JAMMED CONFIG TO .jam FILE" << endl;
				cout << " ENDING COMPRESSION SIMULATION" << endl;
				printJammedConfig();
				break;
			}
		}

		// save last state before packing fraction is changed
		phi = packingFraction();
		saveState(savedState);

		// change packing fraction to new phi value (decided on above)
		phiNew = phi + dphi;
		setPackingFraction(phiNew);

		// update new phi (only update here, do NOT calculate relaxed phi value)
		phi = packingFraction();
	}

	if (k == kmax){
		cout << "	** ERROR: IN 2d cell jamming finding, k reached kmax without finding jamming. Ending." << endl;
		exit(1);
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


// decrease packing fraction quasistatically (ratchet perimeter based on applied forces)
void cellPacking2D::qsIsoGelRatchet(double phiGel, double deltaPhi, double plThresh, double dl0, double calA0max, double timeStepMag){
	// local variables
	double phi0, phiNew, dphi, lenScale;
	int NSTEPS, k, ci;

	// tolerances
	const double Ktol = 1e-10;
	const double Ptol = 1e-6;

	// calculate phi before initial minimization
	phi = packingFraction();

	// relax shapes (energies calculated in relax function)
	cout << "	** IN qsIsoCompression, performing initial relaxation" << endl;
	fireMinimizeGel(Ptol, Ktol);

	// get initial packing fraction
	phi = packingFraction();
	phi0 = phi;

	// determine number of steps to target
	NSTEPS = ceil((phi0 - phiGel)/deltaPhi) + 1;

	// update new dphi to make steps even
	dphi = (phi0 - phiGel)/NSTEPS;

	// update time step
	lenScale = sqrt(4.0*L.at(0)*L.at(1)*phi0/(NCELLS*PI));
	dt = timeStepMag*lenScale;
	dt0 = dt;

	// iterator
	k = 0;

	// loop until phi is the correct value
	while (k < NSTEPS){
		// update iterator
		k++;

		// output to console
		cout << "===================================================" << endl << endl << endl;
		cout << " 	quasistatic isotropic compression to target phi = " << phiGel << endl << endl;
		cout << "===================================================" << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* NSTEPS 		= " << NSTEPS << endl;		
		cout << "	* dphi 			= " << dphi << endl << endl;
		cout << "	AFTER LAST MINIMIZATION:" << endl;
		cout << "	* phi 			= " << phi << endl;
		cout << "	CELL calA0 VALUES:" << endl;
		for (ci=0; ci<NCELLS; ci++){
			if (ci % 5 == 0) 
				cout << endl;
			cout << setw(12) << cell(ci).calA0();			
		}
		cout << endl << endl;

		// increase packing fraction to new phi value
		phiNew = phi0 - k*dphi;
		setPackingFraction(phiNew);

		// calculate phi after minimization
		phi = packingFraction();

		// relax shapes (energies calculated in relax function)
		fireMinimizeGel(Ptol, Ktol);

		// calculate phi after minimization
		phi = packingFraction();

		// update l0 based on segment forces
		ratchetPerimeter(plThresh, dl0, calA0max);
	}
}


// ratchet perimeter based on force
void cellPacking2D::ratchetPerimeter(double plThresh, double dl0, double calA0max){
	// local variables
	int ci, vi;
	double kltmp, l0tmp, pl;
	double fltmp = 0.0;

	// calculate current forces on perimeter 
	for (ci=0; ci<NCELLS; ci++){
		// spring constant
		kltmp = cell(ci).getkl();
		l0tmp = cell(ci).getl0();

		// loop over vertices
		fltmp = 0.0;
		for (vi=0; vi<cell(ci).getNV(); vi++)
			fltmp += kltmp*(cell(ci).segmentLength(vi) - l0tmp);

		// calculate line pressure
		pl = fltmp/cell(ci).perimeter();

		// increment l0 if pl is above threshold and calA0 is below max
		if (pl > plThresh && cell(ci).calA0() < calA0max)
			cell(ci).setl0(l0tmp + dl0);
	}
}





// decrease packing fraction at fixed rate
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
	phi = packingFraction();
	phi0 = phi;
	phitmp = phi0;

	// iterator
	k = 0;

	// loop time until packing fraction below threshold
	while (phitmp > phiGel){
		// // DEBUG: PRINT OUT EVERY TIME WHEN NAN POPS UP
		// if (k > 10000){
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

		// 	if (k > 10500){
		// 		cout << "should have found a nan by now, ending" << endl;
		// 		exit(1);
		// 	}

		// 	cout << endl;
		// 	cout << "k = " << k << endl;
		// 	cout << "phi = " << phi << endl;
		// 	for (ci=0; ci<NCELLS; ci++){
		// 		cout << cell(ci).cpos(0) << "; " << cell(ci).cpos(1) << "; ";
		// 		cout << cell(ci).cvel(0) << "; " << cell(ci).cvel(1) << "; ";
		// 		cout << cell(ci).cforce(0) << "; " << cell(ci).cforce(1);
		// 		cout << endl;
		// 	}
		// }

		// NOTE: IF PARTICLES ARE OVERLY CONCAVE, THEN ERROR OCCURS IN PF UPDATE, ONLY 
		// SCALE BY APPARENT PHI

		// set new target phi
		phi = phitmp;
		phitmp = phi0/(gelRate*t + 1.0);

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
			
			cout << "	* old phi 	= " << phi << endl;
			cout << "	* new phi 	= " << phitmp << endl;
			cout << "	* del phi 	= " << phitmp - phi << endl;
			cout << "	* K 		= " << totalKineticEnergy() << endl;
			cout << "	* U 		= " << totalPotentialEnergy() << endl;
			cout << "	* nc 		= " << totalNumberOfContacts() << endl;
			cout << endl << endl;
		}

		// scale system to decrease packing fraction
		setPackingFraction(phitmp);

		// do verlet update
		for (ci=0; ci<NCELLS; ci++){
			cell(ci).verletPositionUpdate(dt);
			cell(ci).updateCPos();
		}

		// reset contacts before force calculation
		resetContacts();

		// calculate forces
		gelationForces();

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

// decrease packing fraction at fixed rate
void cellPacking2D::gelVarPerimRate(double phiGel, double gelRate, double varPerimRate, double timeStepMag){
	// local variables
	int ci, vi, d, k;
	double t = 0.0;
	double phi0, postmp, phitmp, veltmp;

	// max deformation
	const double calA0max = 1.4;
	double l0max = 0.0;
	double l0tmp, currperim;

	// update time step
	dt = timeStepMag*sqrt(PI);
	dt0 = dt;

	// get initial packing fraction
	phi = packingFraction();
	phi0 = phi;
	phitmp = phi0;

	// iterator
	k = 0;

	// loop time until packing fraction below threshold
	while (phitmp > phiGel){
		// set new target phi
		phi = phitmp;
		phitmp = phi0/(gelRate*t + 1.0);

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
			
			cout << "	* old phi 	= " << phi << endl;
			cout << "	* new phi 	= " << phitmp << endl;
			cout << "	* del phi 	= " << phitmp - phi << endl;
			cout << "	* K 		= " << totalKineticEnergy() << endl;
			cout << "	* U 		= " << totalPotentialEnergy() << endl;
			cout << "	* nc 		= " << totalNumberOfContacts() << endl;
			cout << endl << endl;
		}

		// scale system to decrease packing fraction
		setPackingFraction(phitmp);

		// update positions based on forces (EULER)
		for (ci=0; ci<NCELLS; ci++){
			for (vi=0; vi<cell(ci).getNV(); vi++){
				for (d=0; d<NDIM; d++){
					// updated velocity from forces
					veltmp = cell(ci).vvel(vi,d);

					// update positions (EULER STEP)
					postmp = cell(ci).vpos(vi,d) + dt*veltmp;

					// update positions 
					cell(ci).setVPos(vi,d,postmp);

					// update velocities
					cell(ci).setVVel(vi,d,veltmp);

					// reset forces and energies
					cell(ci).setVForce(vi,d,0.0);
					cell(ci).setUInt(vi,0.0);
				}
			}

			// update cpos
			cell(ci).updateCPos();
		}

		// update preferred perimeter
		for (ci=0; ci<NCELLS; ci++){
			// determine l0 max for this cell
			l0max = sqrt(4.0*PI*cell(ci).geta0()*calA0max)/cell(ci).getNV();

			// current and target l0
			l0tmp = cell(ci).getl0();
			currperim = cell(ci).perimeter()/cell(ci).getNV();

			// set new l0
			l0tmp += dt*varPerimRate*(currperim - l0tmp);

			// check to see whether or not perimeter is extended by too much
			if (l0tmp < l0max && currperim > l0tmp)
				cell(ci).setl0(l0tmp);
		}

		// reset contacts before force calculation
		resetContacts();

		// calculate forces
		gelationForces();

		// update velocities based on forces
		for (ci=0; ci<NCELLS; ci++){
			for (vi=0; vi<cell(ci).getNV(); vi++){
				for (d=0; d<NDIM; d++)
					cell(ci).setVVel(vi,d,cell(ci).vforce(vi,d));
			}
		}
		

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
	const double alpha0 	= 0.1;
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
	double Ptol 			= 1e-12;
	double Ktol 			= 1e-16;
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
	Pvirial = 0.5*(sigmaXX + sigmaYY)/(L.at(0)*L.at(1));

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
		Pvirial = 0.5*(sigmaXX + sigmaYY)/(L.at(0)*L.at(1));

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
		converged = (abs(Pcheck) < Ptol && npPMIN > NMIN);
		converged = (converged || (abs(Pcheck) > 2*Ptol && Kcheck < Ktol));

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

void cellPacking2D::spNVE(vector<double>& lenscales, int nt){
	for (int t=0; t<nt; t++){
		// verlet position update
		spPosVerlet();

		// reset contacts before force calculation
		resetContacts();

		// calculate forces between disks (with door closed)
		spForces(lenscales);

		// verlet velocity update
		spVelVerlet(lenscales);
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
	double Fcheck, Kcheck;

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
			fireMinimizeF(Ftol, Ktol, Fcheck, Kcheck);
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




// VDOS FUNCTIONS

// incremement by dphi, relax, calculate vdos, print to file
void cellPacking2D::cellVDOS(ofstream& vdosOutObj, double dphi, double Ftol, double Ktol){
	// local variables
	int vi, kx, ky, lx, ly;
	double Ftest, Ktest;
	double phiNew;

	// vdos variables

	// calculate current packing fraction
	phi = packingFraction();
	cout << "	** IN cellVDOS, current packing fraction = " << phi << endl;

	// increment packing fraction
	phiNew = phi + dphi;
	setPackingFraction(phiNew);
	cout << "	** IN cellVDOS, incrementing phi by dphi = " << dphi << endl;

	// recalculate new packing fraction
	phi = packingFraction();
	cout << "	** IN cellVDOS, new packing fraction = " << phi << endl;

	// relax shapes (energies calculated in relax function)
	cout << "	** IN cellVDOS, performing relaxation" << endl;
	fireMinimizeF(Ftol, Ktol, Ftest, Ktest);

	// now loop over cells, compute shape stiffness and stress matrices

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
		packingPrintObject << setw(w1) << left << "BOXSZ" << setw(w2) << right << L.at(0) << setw(w2) << right << L.at(1) << endl;
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
	packingPrintObject << setw(w1) << left << "PACKF" << setw(w3) << right << packingFraction() << endl;


	// print hopper information
	packingPrintObject << setw(w1) << left << "BOXSZ";
	packingPrintObject << setw(w3) << right << L.at(0);
	packingPrintObject << setw(w3) << right << L.at(1);
	packingPrintObject << endl;

	// print info for rest of the cells
	for (int ci=0; ci<NCELLS; ci++)
		cell(ci).printVertexPositions(packingPrintObject,ci);

	// print end frame
	packingPrintObject << setw(w1) << left << "ENDFR" << " " << endl;
}

void cellPacking2D::printJammedConfig(){
	// local variables
	int w1 = 12;
	int w2 = 6;
	int w3 = 30;

	// check to see if file is open
	if (!jamPrintObject.is_open()) {
		cout << "	ERROR: jamPrintObject is not open in printJammedConfig(), ending." << endl;
		exit(1);
	}

	// print information starting information
	jamPrintObject << setw(w1) << left << "NEWFR" << " " << endl;
	jamPrintObject << setw(w1) << left << "NUMCL" << setw(w2) << right << NCELLS << endl;
	jamPrintObject << setw(w1) << left << "PACKF" << setw(w3) << right << packingFraction() << endl;

	// print box size information
	jamPrintObject << setw(w1) << left << "BOXSZ";
	jamPrintObject << setw(w3) << right << L.at(0);
	jamPrintObject << setw(w3) << right << L.at(1);
	jamPrintObject << endl;	

	// print contact information
	jamPrintObject << setw(w1) << left << "CTCTS";
	for (int ci=0; ci<NCELLS; ci++){
		for (int cj=ci+1; cj<NCELLS; cj++)
			jamPrintObject << setw(w2) << contacts(ci,cj);
	}
	jamPrintObject << endl;

	// print info for rest of the cells
	for (int ci=0; ci<NCELLS; ci++)
		cell(ci).printVertexPositions(jamPrintObject,ci);

	// print end frame
	jamPrintObject << setw(w1) << left << "ENDFR" << " " << endl;
}

void cellPacking2D::printSystemEnergy(){
	// local variables
	int ci;

	// check to see if file is open
	if (!energyPrintObject.is_open()) {
		cout << "	ERROR: energyPrintObject is not open in printSystemEnergy(), ending." << endl;
		exit(1);
	}

	// print system energies, stress and packing fraction
	energyPrintObject << setw(30) << setprecision(16) << right << interactionPotentialEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << totalPotentialEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << totalKineticEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaXX;
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaXY;
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaYY;
	energyPrintObject << setw(30) << setprecision(16) << right << packingFraction();

	// loop over particles, print cell calA and calA0
	for (ci=0; ci<NCELLS; ci++){
		energyPrintObject << setw(30) << setprecision(16) << right << cell(ci).calA0();
		energyPrintObject << setw(30) << setprecision(16) << right << cell(ci).asphericity();
	}

	// print new line
	energyPrintObject << endl;
}

void cellPacking2D::printSystemEnergy(int intVal){

	// local variables
	int ci;
	double forceScale = cell(0).getka()*cell(0).geta0()/cell(0).getl0();


	// check to see if file is open
	if (!energyPrintObject.is_open()) {
		cout << "	ERROR: energyPrintObject is not open in printSystemEnergy(), ending." << endl;
		exit(1);
	}

	// loop over particles, print cell energy
	energyPrintObject << setw(30) << setprecision(16) << right << interactionPotentialEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << totalPotentialEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << totalKineticEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << forceRMS()/forceScale;
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaXX;
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaXY;
	energyPrintObject << setw(30) << setprecision(16) << right << sigmaYY;
	energyPrintObject << setw(30) << setprecision(16) << right << packingFraction();
	energyPrintObject << setw(30) << setprecision(16) << right << intVal;

	// loop over particles, print cell calA and calA0
	for (ci=0; ci<NCELLS; ci++){
		energyPrintObject << setw(30) << setprecision(16) << right << cell(ci).calA0();
		energyPrintObject << setw(30) << setprecision(16) << right << cell(ci).asphericity();
	}

	// print new line
	energyPrintObject << endl;
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
	statPrintObject << L.at(0) << setw(w) << L.at(1) << endl;
	statPrintObject << setprecision(p) << phi << endl;

	// print contact matrix
	for (ci=0; ci<NCELLS; ci++){
		for (cj=ci+1; cj<NCELLS; cj++)
			statPrintObject << setw(w) << contacts(ci,cj);
	}
	statPrintObject << endl;
}



