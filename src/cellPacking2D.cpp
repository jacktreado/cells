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

// overloaded constryctor with file stream object as argument
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

	// copy scalar variables
	NCELLS 	= onTheRight.NCELLS;
	NT 		= onTheRight.NT;
	NPRINT 	= onTheRight.NPRINT;

	dt 		= onTheRight.dt;
	dt0 	= onTheRight.dt0;
	L 		= onTheRight.L;
	T 		= onTheRight.T;
	phi 	= onTheRight.phi;

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

		// get random direction
		rv = drand48();

		// add to velocities and mean
		for (d=0; d<NDIM; d++){
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
	int ci;
	double val = 0.0;

	// loop over cells, packing fraction is : triangular area + 0.5*delta*perimeter area + area of circular corners
	for (ci=0; ci<NCELLS; ci++){
		val += cell(ci).area() + (0.5*cell(ci).getdel())*cell(ci).perimeter() + PI*0.25*cell(ci).getdel()*cell(ci).getdel();
		// val += cell(ci).area();		// needs to be this if vertex forces only
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
	double scaleFactor, dphi;

	// check if packing fraction has not yet been set
	if (phi <= 0.0)
		phi = packingFraction();

	// calculate val to scale lengths with
	dphi 			= val - phi;
	scaleFactor		= pow((phi+dphi)/phi,1.0/NDIM);

	// scale all lengths by scale factor
	scaleLengths(scaleFactor);

	// update new phi
	phi = val;
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


/************************

	Printers

*************************/


// print positions to file
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

// print energies to file
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
	energyPrintObject << setw(30) << setprecision(16) << right << interactionPotentialEnergy()/cell(0).getkint();
	energyPrintObject << setw(30) << setprecision(16) << right << totalPotentialEnergy()/cell(0).getkint();

	// ADD ON STUFF FOR DEBUGGING JAMMING FINDER
	double calA0, meanCalA;

	calA0 = (pow(cell(0).getl0(),2)*pow(cell(0).getNV(),2))/(4*PI*cell(0).geta0());
	meanCalA = cell(0).asphericity();

	energyPrintObject << setw(30) << setprecision(16) << right << Pval;
	energyPrintObject << setw(30) << setprecision(16) << right << calA0;
	energyPrintObject << setw(30) << setprecision(16) << right << meanCalA;
	energyPrintObject << endl;

	// NOTE: HOW TO ADD ENERGIES OF INDIVIDUAL CELLS?
	// for (i=0; i<NCELLS; i++)
	// 	cell(i).printCellEnergy(energyPrintObject,frame);
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




/**************************

	Forces and position 
		updates

***************************/


// calculate all forces, both shape and pairwise
void cellPacking2D::calculateForces(){
	// local variables
	int ci,cj,inContact;

	// vector to store pairwise forces
	vector<double> fij(NDIM,0.0);
	vector<double> rij(NDIM,0.0);

	// reset virial stresses to 0
	sigmaXX = 0.0;
	sigmaXY = 0.0;
	sigmaYX = 0.0;
	sigmaYY = 0.0;

	// loop over cells and cell pairs, calculate shape and interaction forces
	for (ci=0; ci<NCELLS; ci++){
		// loop over pairs, add info to contact matrix
		for (cj=ci+1; cj<NCELLS; cj++){
			inContact = cell(ci).vertexForce(cell(cj),fij,rij);
			if (inContact == 1)
				addContact(ci,cj);

			// compute virial stresses
			sigmaXX += fij.at(0)*rij.at(0);
			sigmaXY += fij.at(0)*rij.at(1);
			sigmaYX += fij.at(1)*rij.at(0);
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
	const double alpha0 	= 0.1;
	const double finc 		= 1.1;
	const double fdec 		= 0.5;
	const double falpha 	= 0.999;
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
		fireMinimize(Ptol, Ktol, plotIt, k);
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
					fireMinimize(Ptol, Ktol, plotIt, k);
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


// fire energy minimzation with backstepping if P < 0
void cellPacking2D::fireMinimize(double Ptol, double Ktol, int plotIt, int& frameCount){
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
	double Knew, Pvirial;

	// variable to test for potential energy minimization
	bool converged = false;

	// reset time step
	dt = dt0;

	// initialize virial pressure from pressure from last time
	Pvirial = 0.5*(sigmaXX + sigmaYY);

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
					printSystemEnergy(frameCount,Pvirial/energyScale,Knew/energyScale);
				}
				frameCount++;
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << Knew << endl;
			cout << "	* Pvirial 	= " << Pvirial/(Ptol*energyScale) << endl;
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
		Pvirial = 0.5*(sigmaXX + sigmaYY);
		energyScale = cell(0).getkint();

		// update if Pvirial under tol
		if (abs(Pvirial) < Ptol*energyScale)
			npPMIN++;
		else
			npPMIN = 0;

		// check for convergence
		converged = (abs(Pvirial) < Ptol*energyScale && npPMIN > NMIN);
		converged = (converged || (abs(Pvirial) > 2*Ptol*energyScale && Knew < NCELLS*cell(0).getNV()*energyScale*Ktol));

		if (converged){
			cout << "	** FIRE has converged!" << endl;
			cout << "	** k = " << k << ", t = " << t << endl;
			cout << "	** virial P = " << Pvirial << endl;
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






// NON-EQUILIBRIUM MD FUNCTIONS


// quasistaic gel forming function
void cellPacking2D::isoExtensionQS(int plotIt, int& frameCount, double phiTarget, double dphi){
	// local variables
	int t = 0;
	int isRelaxed = 0;
	double Ptol = 1e-8;
	double Ktol = 1e-24;

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
			fireMinimize(Ptol, Ktol, plotIt, frameCount);
			isRelaxed = 1;
		}

		// increment t
		t++;
	}

	if (phi < phiTarget)
		cout << "	Target phi = " << phiTarget << "reached! Ending isoExtensionQS protocol" << endl;
}










// TUMOR TISSUE FUNCTIONS


// NVE dynamics
void cellPacking2D::tumorNVE(){
	// local variables
	int t, ci, frameCount;
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

		// update velocities
		for (ci=0; ci<NCELLS; ci++)
			cell(ci).verletVelocityUpdate(dt,0.0);
	}

	// print something to conclude
	cout << "t = NT, which = " << t << ", so ending NVE run!" << endl;
}


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
void cellPacking2D::overlapRelief(){
	// local variables
	int ci,cj,k,kmax,fcheck,fconst,inContact;
	int frameCount = 0;
	double Fnew,Fold,dF,Ftol,dampingParameter;

	// interaction variables
	double dist2, cdist2;

	// initial variables for checking for constant U
	fcheck = 0;
	fconst = 0;
	Fold = 0.0;
	Ftol = 1e-8;

	// damping
	dampingParameter = 0.3;

	// max number of iterations
	kmax = 1e5;

	// loop over system and push cells aways from one another 
	// as damped frictionless disks until force has plateau'd constant
	while (fconst == 0 && k < kmax){
		// update particle positions
		for (ci=0; ci<NCELLS; ci++){
			cell(ci).verletPositionUpdate(20*dt);
			cell(ci).updateCPos();
		}

		// update forces based on center-to-center distances only
		for (ci=0; ci<NCELLS; ci++){
			// loop over pairs, add info to contact matrix
			for (cj=ci+1; cj<NCELLS; cj++){
				inContact = cell(ci).radialForce(cell(cj));
				if (inContact == 1)
					addContact(ci,cj);
			}
		}

		// update damped velocities
		for (ci=0; ci<NCELLS; ci++)
			cell(ci).verletVelocityUpdate(20*dt,dampingParameter);

		// check that F is below threshold and constant
		Fnew = maxForceMagnitude();
		dF = Fnew - Fold;
		Fold = Fnew;
		if (abs(dF) < Ftol){
			fcheck++;
			if (fcheck > 50){
				fconst = 1;
				fcheck = 0;
			}
			
		}
		else{
			fcheck = 0;
			fconst = 0;
		}

		// increment k
		k++;
	}

	if (k == kmax){
		cout << "	ERROR: overlapRelief() function failed, could not reduce potential energy, ending." << endl;
		exit(1);
	}
}

// increase attraction quasi-statically, relax after each add
void cellPacking2D::attractionRamp(double attractionTarget, double dAttraction, int plotIt, int initalFrame){
	// local variables
	int k, ci, delSgn, frameCount;
	double da, currAttraction, decScale, Ptol, Ktol;
	double attractionCheck, tol;

	// set check and tolerance
	tol = 1e-8;
	attractionCheck = 10*tol;

	// set tolerance for relaxation
	Ktol = 1e-24;
	Ptol = 1e-6;

	// get current attraction (assuming all attraction the same)
	currAttraction = cell(0).geta();

	// set decrease scale to be 1
	decScale = 1.0;

	// check whether to increase or decrease
	if (currAttraction < attractionTarget)
		delSgn = 1;
	else
		delSgn = -1;

	// loop until attraction is the correct value
	k=0;
	frameCount=initalFrame;
	while (attractionCheck > tol && k < NT){
		// output to console
		cout << "===================================================" << endl << endl << endl;
		cout << " 	Ramping attraction to aTarget = " << attractionTarget << endl << endl;
		cout << "===================================================" << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* current a 	= " << currAttraction << endl;
		cout << endl << endl;

		// decide whether attraction is too large or too small
		if (delSgn*currAttraction < delSgn*attractionTarget)
			da = delSgn*dAttraction;
		else{
			decScale *= 0.9;
			da = -delSgn*decScale*dAttraction;
		}

		// update current attraction parameter
		currAttraction += da;

		// update attraction in all of the cells
		for (ci=0; ci<NCELLS; ci++)
			cell(ci).seta(currAttraction);

		// relax potential energy
		// potentialRelaxFire(Ktol, Utol, 1, frameCount);
		fireMinimize(Ptol, Ktol, plotIt, frameCount);

		// update check
		attractionCheck = abs(currAttraction - attractionTarget);

		// update iterator
		k++;
	}

	if (k == NT){
		cout << "	ERROR: particle shapes could not relax to desired attraction in allotted time, ending." << endl;
		exit(1);
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
		fireMinimize(Ptol, Ktol, 1, frameCount);

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






