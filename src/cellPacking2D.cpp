/*

	Methods file for cellPacking2D class

*/


// include file
#include "deformableParticles2D.h"
#include "cellPacking2D.h"

// namespace
using namespace Eigen;
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
	NDIM 			= 2;

	// scalar variables set to 0
	NCELLS 			= 0;
	NT 				= 0;
	NPRINT 			= 0;
	dt 				= 0.0;
	dt0 			= 0.0;
	T 				= -1.0;
	phi 			= -1.0;
	sigmaXX 		= 0.0;
	sigmaXY 		= 0.0;
	sigmaYX 		= 0.0;
	sigmaYY 		= 0.0;
	shearStrain 	= 0.0;

	// set box lengths to 1.0
	L.resize(NDIM);
	for (int d=0; d<NDIM; d++)
		L.at(d) 	= 1.0;

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
}

// overloaded constructor with file stream object as argument
cellPacking2D::cellPacking2D(string& inputFile, double T0, double s){
	// set initial seed
	seed = s;
	srand(seed);

	// set variables to default
	defaultvars();

	// local variables
	int NC, ci, vi, nv, d;
	double lxtmp, lytmp, x, y;
	double l0tmp, a0tmp, deltmp;

	// open input file
	ifstream inputFileObject;
	inputFileObject.open(inputFile.c_str());
	if (!inputFileObject.is_open()) {
		cout << "	ERROR: inputFileObject could not open " << inputFile << "..." << endl;
		exit(1);
	}

	// string readin
	string inputStr;

	// LINE 1: should be NEWFR
	getline(inputFileObject, inputStr);
	if (inputStr.compare(0,5,"NEWFR") != 0){
		cout << "	** first line of input file NOT NEWFR, inputStr = " << inputStr << ". Ending." << endl;
		exit(1);
	}

	// read in simulation information
	getline(inputFileObject, inputStr);
	sscanf(inputStr.c_str(),"NUMCL %d",&NCELLS);
	
	getline(inputFileObject, inputStr);
	sscanf(inputStr.c_str(),"PACKF %lf",&phi);

	getline(inputFileObject, inputStr);
	sscanf(inputStr.c_str(),"BOXSZ %lf %lf",&lxtmp,&lytmp);

	getline(inputFileObject, inputStr);
	sscanf(inputStr.c_str(),"NCONTS %d %d",&Ncc,&Nvv);


	cout << "Finished reading header in read-in operator, values so far are: " << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cout << "	** phi = " << phi << endl;
	cout << "	** L = " << lxtmp << endl;
	cout << "	** Ncc = " << Ncc << ", Nvv = " << Nvv << endl;

	// set box lengths to be square, reset lengths to be in units of l0
	L.at(0) = lxtmp;
	L.at(1) = lytmp;

	// test for error in inputs
	if (NCELLS <= 0){
		cout << "	ERROR: in overloaded read-in operator, NCELLS <= 0 so cannot initialize arrays, ending code here." << endl;
		exit(1);
	}
	else if (L.at(0) <= 0.0 || L.at(1) <= 0.0){
		cout << "	ERROR: in overloaded read-in operator, L <= 0.0 so cannot run simulations, ending code here." << endl;
		exit(1);
	}

	// test that memory has not yet been initialized
	if (contactMatrix){
		cout << "	ERROR: in overloaded read-in operator, contactMatrix ptr already initialized, ending code here." << endl;
		exit(1);
	}

	// initialize cell array
	cellArray = new deformableParticles2D[NCELLS];

	// loop over cells
	for (ci=0; ci<NCELLS; ci++){
		// first read number of vertices
		getline(inputFileObject, inputStr);
		sscanf(inputStr.c_str(),"NVERT %d",&nv);

		// initialize cell objects ( set PBCs = 1 )
		for (d=0; d<NDIM; d++){
			cell(ci).setL(d,L.at(d));
			cell(ci).setpbc(d,1);
		}

		// initialize cell object
		cell(ci).setNV(nv);
		cell(ci).initializeVertices();
		cell(ci).initializeCell();

		// read in cell com information
		getline(inputFileObject, inputStr);
		sscanf(inputStr.c_str(),"CELLP %*d %lf %lf %lf %lf %lf %*lf %*lf %*lf", &x, &y, &l0tmp, &a0tmp, &deltmp);

		// set cell com position
		cell(ci).setCPos(0,x);
		cell(ci).setCPos(1,y);

		// set cell shape information
		cell(ci).seta0(a0tmp);
		cell(ci).setl0(l0tmp);
		cell(ci).setdel(deltmp);

		// print to console
		cout << "	** on cell " << ci << ", NV = " << nv << ";\t com at x = " << x << ", y = " << y;
		cout << ";\t l0 = " << l0tmp << ", a0 = " << a0tmp << ", so calA0 = " << (nv*nv*l0tmp*l0tmp)/(4.0*PI*a0tmp) << endl;

		// read in descriptive string
		getline(inputFileObject, inputStr);
		cout << "header check : " << inputStr << endl;

		// set vertex positions
		for (vi=0; vi<nv; vi++){
			// read info for each vertex
			getline(inputFileObject, inputStr);
			sscanf(inputStr.c_str(),"VERTP %*d %lf %lf %*lf %*lf %*lf %*lf", &x, &y);

			// set vertex positions
			cell(ci).setVPos(vi,0,x);
			cell(ci).setVPos(vi,1,y);

			// print to console
			cout << "	-- -- cell " << ci << ", vertex " << vi << ": x = " << x << ", y = " << y << endl;
		}
	}

	// initialize contact matrix
	NC = NCELLS*(NCELLS-1)/2;
	contactMatrix = new int[NC];

	// intialize contact matrix to 0
	for (ci=0; ci<NC; ci++)
		contactMatrix[ci] = 0;

	// initialize velocities
	initializeVelocities(T0);

	// initialize forces
	calculateForces();

	// print cell information
	for (ci=0; ci<NCELLS; ci++){
		cout << "CELLP\t" << ": NV = " << cell(ci).getNV() << setw(6) << "; cx = " << cell(ci).cpos(0) << setw(6) << "; cy = " << cell(ci).cpos(1) << endl;
		for (vi=0; vi<cell(ci).getNV(); vi++){
			cout << "VERTP\t";
			cout << cell(ci).vpos(vi,0) << "; ";
			cout << cell(ci).vpos(vi,1) << "; ";
			cout << cell(ci).vvel(vi,0) << "; ";
			cout << cell(ci).vvel(vi,1) << "; ";
			cout << cell(ci).vforce(vi,0) << "; ";
			cout << cell(ci).vforce(vi,1);
			cout << endl;
		}
		cout << endl << endl << endl;
	}

	// close inputFile object
	inputFileObject.close();
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

	sigmaXX 	= onTheRight.sigmaXX;
	sigmaXY 	= onTheRight.sigmaXY;
	sigmaYX 	= onTheRight.sigmaYX;
	sigmaYY 	= onTheRight.sigmaYY;
	shearStrain = onTheRight.shearStrain;

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

		saveObject.shearStrain = shearStrain;

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

// load saved state from loadObject
void cellPacking2D::loadState(cellPacking2D& loadObject){
	// local variables
	int ci, cj;

	// test that loadObject has been initialized
	if (loadObject.phi < 0){
		cout << "	** ERROR: trying to load from loadObject that has no saved data, ending." << endl;
		exit(1);
	}

	// load saved packing fraction, dt
	dt = loadObject.dt;
	dt0 = loadObject.dt0;
	phi = loadObject.phi;

	sigmaXX = loadObject.sigmaXX;
	sigmaXY = loadObject.sigmaXY;
	sigmaYX = loadObject.sigmaYX;
	sigmaYY = loadObject.sigmaYY;

	shearStrain = loadObject.shearStrain;

	// load cell and contact data from loadObject
	for (ci=0; ci<NCELLS; ci++){
		// copy cell objects (using overloaded operator in deformableParticle2D class)
		cell(ci) = loadObject.cell(ci);

		// copy elements from contact matrix
		for (cj=ci+1; cj<NCELLS; cj++){
			if (loadObject.contacts(ci,cj))
				addContact(ci,cj);
		}
	}
}



/************************

	Initialization

*************************/

// initialize bidisperse packing 
void cellPacking2D::initializeBidisperse(int NV, double phiDisk, double sizeRatio, double sizeFraction, double delval){
	// local variables
	int ci, vi, d, nvtmp;
	int lx, ly;
	double calA0tmp;
	double xpos, ypos, dx, dy;
	double xmin, xmax, ymin, ymax;
	double areaSum;
	double rtmp, a0tmp, l0tmp;

	// minimum number of vertices
	const int nvmin = 12;

	// output to console
	cout << "		-- In bidisperse initialization, initializing cells and relaxing initial overlaps as repulsive SP particles" << endl;

	// initialize length scales
	areaSum = 0.0;
	vector<double> lenscales(NCELLS,0.0);
	vector<double> diskradii(NCELLS,0.0);
	for (ci=0; ci<NCELLS; ci++){
		// store length scale depending on size fraction
		if (ci < round(sizeFraction*NCELLS))
			lenscales.at(ci) = 1.0;
		else
			lenscales.at(ci) = sizeRatio;

		// effective particle area
		a0tmp = lenscales.at(ci)*lenscales.at(ci);

		// add to lenscales sum for boundary size
		areaSum += a0tmp;

		// save disk radius based on area square root
		diskradii.at(ci) = 1.05*sqrt(a0tmp/(NV*sin(PI/NV)*cos(PI/NV)));
	}

	// determine box length from particle sizes and input packing fraction
	for (d=0; d<NDIM; d++)
		L.at(d) = sqrt(areaSum/phiDisk);

	// set phi to input
	phi = phiDisk;

	// seed rng
	srand48(56835698*seed);

	// initialize cell information
	cout << "		-- Ininitializing cell objects" << endl;
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

		// initialize cells as regular polygons
		calA0tmp = nvtmp*tan(PI/nvtmp)/PI;

		// preferred area is lenscales squared
		a0tmp = lenscales.at(ci)*lenscales.at(ci);

		// initial length of polygon side
		l0tmp = 2.0*lenscales.at(ci)*sqrt(PI*calA0tmp)/nvtmp;

		// set preferred area and length 
		cell(ci).seta0(a0tmp);
		cell(ci).setl0(l0tmp);
		cell(ci).setdel(delval);
	}

	// initialize particle positions
	cout << "		-- Ininitializing cell positions" << endl;

	// set min and max values of positions
	xmin = 0;
	xmax = L.at(0);
	ymin = 0;
	ymax = L.at(1);

	// initialize positions of each cell
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

	// set time step
	vertexDPMTimeScale(0.1);

	// use FIRE in PBC box to relax overlaps
	cout << "		-- Using FIRE to relax overlaps..." << endl;
	fireMinimizeSP(diskradii);

	// print cell centers of mass as a check
	cout << "		-- FIRE relax overlap complete, cell centers of mass = " << endl;
	for (ci=0; ci<NCELLS; ci++){
		cout << cell(ci).cpos(0) << "  " << cell(ci).cpos(1) << endl;
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
			if (NCELLS > 1)
				cell(ci).setCVel(d,cell(ci).cvel(d) - vmean.at(d));

			// calc ek
			ek += 0.5*pow(cell(ci).cvel(d),2);
		}
	}

	// get vscale
	vscale = sqrt(T/ek);
	for (ci=0; ci<NCELLS; ci++){
    	for (d=0; d<NDIM; d++)
        	cell(ci).setCVel(d,cell(ci).cvel(d)*vscale);
    }
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
	frms = sqrt(frms)/(NDIM*NVTOTAL);

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

// set forces values for all cells to be same
void cellPacking2D::forceVals(double calA0, double ka, double kl, double gam, double kb, double kint, double del, double a){
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

// set vertex DPM time step mag ( which is 1, as area force = 1)
void cellPacking2D::vertexDPMTimeScale(double timeStepMag){
	dt = timeStepMag;
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
	int ci,cj,vi,d,dd,inContact;

	// reset virial stresses to 0
	sigmaXX = 0.0;
	sigmaXY = 0.0;
	sigmaYX = 0.0;
	sigmaYY = 0.0;

	// reset contacts before force calculation
	resetContacts();
	Ncc = 0;
	Nvv = 0;

	// reset forces
	for (ci=0; ci<NCELLS; ci++){
		// reset center of mass forces
		for (d=0; d<NDIM; d++)
			cell(ci).setCForce(d,0.0);

		// reset vertex forces and interaction energy
		for (vi=0; vi<cell(ci).getNV(); vi++){
			// forces
			for (d=0; d<NDIM; d++)
				cell(ci).setVForce(vi,d,0.0);

			// energies
			cell(ci).setUInt(vi,0.0);
		}
	}


	// loop over cells and cell pairs, calculate shape and interaction forces
	for (ci=0; ci<NCELLS; ci++){
		// loop over pairs, add info to contact matrix
		for (cj=ci+1; cj<NCELLS; cj++){
			// calculate forces, add to number of vertex-vertex contacts
			inContact = cell(ci).vertexForce(cell(cj),sigmaXX,sigmaXY,sigmaYX,sigmaYY);
			if (inContact > 0){
				// add to cell-cell contacts
				addContact(ci,cj);
				Ncc++;

				// increment vertex-vertex contacts
				Nvv += inContact;
			}
		}

		// forces on vertices due to shape
		cell(ci).shapeForces();
	}
}

void cellPacking2D::gelationForces(){
	// local variables
	int ci,cj,vi,d,dd,inContact,numACtmp;
	double aij;

	// vector to store number of attractive contacts per cell
	vector<int> nac(NCELLS,0);

	// reset virial stresses to 0
	sigmaXX = 0.0;
	sigmaXY = 0.0;
	sigmaYX = 0.0;
	sigmaYY = 0.0;

	// reset contacts before force calculation
	resetContacts();
	Ncc = 0;
	Nvv = 0;

	// reset forces
	for (ci=0; ci<NCELLS; ci++){
		// reset center of mass forces
		for (d=0; d<NDIM; d++)
			cell(ci).setCForce(d,0.0);

		// reset vertex forces and interaction energy
		for (vi=0; vi<cell(ci).getNV(); vi++){
			// forces
			for (d=0; d<NDIM; d++)
				cell(ci).setVForce(vi,d,0.0);

			// energies
			cell(ci).setUInt(vi,0.0);
		}
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
			// if attractive contact engaged
			if (nac.at(ci) > 0 && nac.at(cj) > 0){
				// get effective attraction scale (max is 0.5*a)
				aij = 0.25*(cell(ci).geta()/nac.at(ci) + cell(cj).geta()/nac.at(cj));

				// calculate forces
				inContact = cell(ci).vertexForce(cell(cj),sigmaXX,sigmaXY,sigmaYX,sigmaYY,aij);
			}
			// else, use normal force routine
			else
				inContact = cell(ci).vertexForce(cell(cj),sigmaXX,sigmaXY,sigmaYX,sigmaYY);
			
			if (inContact > 0){
				// add to cell-cell contacts
				addContact(ci,cj);
				Ncc++;

				// increment vertex-vertex contacts
				Nvv += inContact;
			}
		}

		// forces on vertices due to shape
		cell(ci).shapeForces();
	}
}


/**************************

	FIRE energy minimzation

***************************/

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
void cellPacking2D::fireMinimizeF(double Ftol, double& Fcheck, double& Kcheck){
	// HARD CODE IN FIRE PARAMETERS
	const double alpha0 	= 0.3;
	const double finc 		= 1.1;
	const double fdec 		= 0.5;
	const double falpha 	= 0.99;
	const double dtmax 		= 10*dt0;
	const double dtmin 		= 1e-8*dt0;
	const double Trescale 	= 1e-8*NCELLS;
	const int NMIN 			= 20;
	const int NNEGMAX 		= 2000;
	const int NDELAY 		= 1000;
	int npPos				= 0;
	int npNeg 				= 0;
	int npPMIN				= 0;
	double alpha 			= 0.0;
	double t 				= 0.0;
	double P 				= 0.0;

	// local variables
	int ci,vi,d,k,kmax;
	double vstarnrm,fstarnrm,vtmp,ftmp;
	double K, F, Pcheck;
	double xold, xnew, vold, vnew, fold;

	// variable to test for potential energy minimization
	bool converged = false;

	// reset time step
	dt = dt0;

	// initialize forces
	calculateForces();

	// rescale velocities
	rescaleVelocities(Trescale);

	// norm of total force vector, kinetic energy
	F = forceRMS();
	K = totalKineticEnergy();
	Pcheck = 0.5*(sigmaXX + sigmaYY)/(NCELLS*L.at(0)*L.at(1));

	// scale P and K for convergence checking
	Fcheck = F;
	Kcheck = K/NCELLS;

	// iterate until system converged
	kmax = 1e6;
	for (k=0; k<kmax; k++){
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
			cout << "	* P 		= " << P << endl;
			cout << "	* Pdir 		= " << P/(vstarnrm*fstarnrm) << endl;
			cout << endl << endl;
		}


		// Step 2. Adjust simulation based on net motion of system
		if (P > 0){
			// increment pos counter
			npPos++;

			// reset neg counter
			npNeg = 0;

			// alter sim if enough positive steps taken
			if (npPos > NMIN){
				// change time step
				if (dt*finc < dtmax)
					dt *= finc;

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
			if (k > NDELAY){
				// decrease time step 
				if (dt*fdec > dtmin)
					dt *= fdec;

				// change alpha
				alpha = alpha0;
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
						vtmp = (1 - alpha)*cell(ci).vvel(vi,d) + alpha*(cell(ci).vforce(vi,d)/fstarnrm)*vstarnrm;
						cell(ci).setVVel(vi,d,vtmp);
					}
				}
			}
		}

		// VV update in FIRE 2.0: position update
		for (ci=0; ci<NCELLS; ci++){
			cell(ci).verletPositionUpdate(dt);
			cell(ci).updateCPos();
		}

		// calculate forces
		calculateForces();

		// VV update in FIRE 2.0: Velocity update 2
		for (ci=0; ci<NCELLS; ci++)
			cell(ci).verletVelocityUpdate(dt);

		// update t
		t += dt;

		// track energy and forces
		F = forceRMS();
		K = totalKineticEnergy();
		Pcheck = 0.5*(sigmaXX + sigmaYY)/(NCELLS*L.at(0)*L.at(1));

		// scale P and K for convergence checking
		Fcheck = F;
		Kcheck = K/NCELLS;

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
		converged = (abs(Fcheck) < Ftol && npPMIN > NMIN);

		if (converged){
			cout << "	** FIRE has converged!" << endl;
			cout << "	** Fcheck = " << Fcheck << endl;
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
				printSystemEnergy(1);
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
		cellRK4();
	}

	// print something to conclude
	cout << "t = NT, which = " << t << ", so ending NVE run!" << endl;
}



// compress isotropically to jamming
void cellPacking2D::findJamming(double dphi0, double Ftol, double Ptol){
	// local variables
	double Ptest, Ktest, Ftest;
	int NSTEPS, k, kmax, kr, nc, nr, ci, cj;
	cellPacking2D savedState;
	int NDOF;

	// get total number of degrees of freedom
	NDOF = 0;
	for (ci=0; ci<NCELLS; ci++)
		NDOF += NDIM*cell(ci).getNV();

	// get initial packing fraction
	phi = packingFraction();

	// iterator
	k = 0;
	kmax = 1e5;

	// jamming variables
	bool jammed, overcompressed, undercompressed;
	double rH, r0, rL, dr0, scaleFactor;

	// compute first dr0 based on current phi (i.e. non root search)
	dr0 = sqrt((phi+dphi0)/phi);

	// save initial state
	r0 = sqrt(cell(0).geta0());
	saveState(savedState);

	// initialize as unjammed
	jammed = false;

	// phiJ bounds
	rH = -1;
	rL = -1;

	// initialize velocities
	double Tinit = 1e-6;
	initializeVelocities(Tinit);

	// loop until phi is the correct value
	while (!jammed && k < kmax){
		// update iterator
		k++;

		// relax shapes (energies/forces calculated during FIRE minimization)
		fireMinimizeF(Ftol, Ftest, Ktest);

		// update new phi after minimization
		phi = packingFraction();

		// calculate Ptest for comparison
		Ptest = 0.5*(sigmaXX + sigmaYY)/(NDOF*L.at(0)*L.at(1));

		// remove rattlers
		kr = 0;
		nr = removeRattlers(kr);

		// update number of contacts
		nc = totalNumberOfContacts();

		// boolean checks
		undercompressed = ((Ptest < 2.0*Ptol && rH < 0) || (Ptest < Ptol && rH > 0));
		overcompressed = (Ptest > 2.0*Ptol && nc > 0);
		jammed = (Ptest < 2.0*Ptol && Ptest > Ptol && nc > 0 && rH > 0);

		// output to console
		cout << "===================================================" << endl << endl << endl;
		cout << " 	quasistatic isotropic compression to jamming " << endl << endl;
		cout << "===================================================" << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* dphi 			= " << dphi0 << endl;
		cout << "	* phi 			= " << phi << endl;
		cout << "	* r0 			= " << r0 << endl;
		cout << "	* rH 			= " << rH << endl;
		cout << "	* rL 			= " << rL << endl;
		cout << "	* Ftest 		= " << Ftest << endl;
		cout << "	* Ktest 		= " << Ktest << endl;
		cout << "	* Ptest 		= " << Ptest << endl;
		cout << "	* # of contacts = " << nc << endl;
		cout << "	* # of rattlers = " << nr << endl;
		cout << "	* undercompressed = " << undercompressed << endl;
		cout << "	* overcompressed = " << overcompressed << endl;
		cout << "	* jammed = " << jammed << endl << endl;
		cout << "	* contact matrix:" << endl;

		// print contact matrix to console
		for (ci=0; ci<NCELLS; ci++){
			for (cj=0; cj<NCELLS; cj++){
				if (cj == ci)
					cout << "0" << "  ";
				else
					cout << contacts(ci,cj) << "  ";
			}
			cout << endl;
		}

		// print final two lines
		cout << endl << endl;

		// update packing fraction based on jamming check
		if (rL < 0){
			// if still undercompressed, then grow until overjammed found
			if (undercompressed){
				// set scale to normal compression
				scaleFactor = dr0;

				// save state
				r0 = sqrt(cell(0).geta0());
				saveState(savedState);
			}
			// if first overcompressed, return to pre-overcompression state, to midpoint between phi and phiH
			else if (overcompressed){
				// current = upper bound length scale r
	            rH = sqrt(cell(0).geta0());
	            
	            // old = old length scale
	            rL = r0;

	            // save overcompressed state
	            loadState(savedState);

	            // compute new scale factor
	            scaleFactor = 0.5*(rH + rL)/r0;

	            // print to console
				cout << "	-- -- overcompressed for first time, scaleFactor = " << scaleFactor << endl;
			}
		}
		else{
			// if found undercompressed state, go to state between undercompressed and last overcompressed states (from saved state)
			if (undercompressed){
				// current = new lower bound length scale r
				rL = sqrt(cell(0).geta0());

				// load state
				loadState(savedState);

				// compute new scale factor
	            scaleFactor = 0.5*(rH + rL)/r0;

				// print to console
				cout << "	-- -- undercompressed, scaleFactor = " << scaleFactor << endl;

			}
			else if (overcompressed){
				// current = new upper bound length scale r
	            rH = sqrt(cell(0).geta0());

				// load state
				loadState(savedState);

				// compute new scale factor
	            scaleFactor = 0.5*(rH + rL)/r0;

				// print to console
				cout << "	-- -- overcompressed, scaleFactor = " << scaleFactor << endl;
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

		// grow or shrink particles by scale factor
		scaleLengths(scaleFactor);
	}

	if (k == kmax){
		cout << "	** ERROR: IN 2d cell jamming finding, k reached kmax without finding jamming. Ending." << endl;
		exit(1);
	}
}


// compress isotropically to fixed packing fraction by ~deltaPhi steps (use length scaler instead for robustness)
void cellPacking2D::qsIsoCompression(double phiTarget, double deltaPhi, double Ftol){
	// local variables
	double dr, phiNew, dphi, Fcheck, Kcheck;
	int kmax, k;

	// get initial packing fraction
	phi = packingFraction();

	// compute length scaler based on deltaPhi
	dr = sqrt((phi + deltaPhi)/phi);

	// iterator
	k = 0;
	kmax = 1e6;

	// loop until phi is the correct value
	while (phi < phiTarget && k < kmax){
		// update iterator
		k++;

		// scale lengths
		scaleLengths(dr);

		// relax shapes (energies calculated in relax function)
		fireMinimizeF(Ftol, Fcheck, Kcheck);

		// update packing fraction
		phi = packingFraction();

		// output to console
		cout << "===================================================" << endl << endl << endl;
		cout << " 	quasistatic isotropic compression " << endl;
		cout << "===================================================" << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* dphi 			= " << deltaPhi << endl;
		cout << "	* phi 			= " << phi << endl;
		cout << "	* Fcheck 		= " << Fcheck << endl;
		cout << "	* Kcheck 		= " << Kcheck << endl;
		cout << endl;
		cout << "	* * distance to target : " << "phiTarget = " << phiTarget << ", distance = " << phiTarget - phi << endl;
		cout << endl << endl;
	}
	if (k == kmax){
		cout << "	** IN qsIsoCompression, compression iteration did not converge in kmax = " << kmax << " iterations. Ending. " << endl;
		exit(1);
	}
}


// compute instantaneous shear modulus with energy minimization
// NOTE: ASSUME STARTING CONFIGURATION IS JAMMED
double cellPacking2D::shearModulus(){
	// local variables
	int ci, vi, d, im;
	double x0, y0, x1;
	double sxy0, sxy, G;
	cellPacking2D savedState;

	// strain step variables
	const int nstrains = 20;
	double gamma = 0.0;
	double dgamma = shearStrain;
	int straini;

	// tolerances for relaxation
	const double Ftol = 1e-9;
	const double Ktol = 1e-11;

	// shear modulus variables
	G = 0.0;
	sxy0 = sigmaXY;

	// check variables for relaxation
	double Fcheck, Kcheck;

	// save initial state
	saveState(savedState);

	// take small strain steps, G is average slope of shear stress
	for (straini=0; straini < nstrains; straini++){
		// increment strain
		gamma += dgamma;

		// set strain parameter in individual cells, will activate LEbc
		for (ci=0; ci<NCELLS; ci++){
			// set shear strain parameter for each cell for bc's
			cell(ci).setstrain(gamma);

			// implement first shear step
			for (vi=0; vi<NCELLS; vi++){
				// initial particle positions
				x0 = cell(ci).vpos(vi,0);
				y0 = cell(ci).vpos(vi,1);

				// affine displacement in shear profile
				x1 = x0 + gamma*y0;

				// update xposition
				cell(ci).setVPos(vi,0,x1);
			}
		}

		// relax using fire minimization
		fireMinimizeF(Ftol, Fcheck, Kcheck);

		// add to shear modulus average
		sxy 	= sigmaXY;
		G 		+= -(sxy - sxy0)/dgamma;
		sxy0 	= sxy;
	}

	// reset positions back to initial unsheared state
	for (ci=0; ci<NCELLS; ci++)
		cell(ci).setstrain(0);

	// load initial state
	loadState(savedState);

	// take average of G
	G /= (double)nstrains;

	// return shear modulus
	return G;
}


// compute the vibrational density of states of a given configuration of cells
void cellPacking2D::vdos(){
	// first check to see if stat object is open; if not, do not compute VDOS
	if (!statPrintObject.is_open()) {
		std::cout << "	ERROR: statPrintObject is not open in VDOS function, ending. " << endl;
		exit(1);
	}

	// LOCAL VARIABLES

	// integers
	int ci, cj, vi, vj, nv, k, l;
	int vim2, vim1, vip1, vip2, vjm1, vjp1;
	int kx, kxp1, kxp2, ky, kyp1, kyp2, lx, ly;
	int mxi, myi, mxj, myj;
	int NDOF = 0;
	vector<int> Mu(NCELLS,0);
	cellPacking2D relaxedState;

	// doubles
	double calA0, l0, fl, kl, kb, eb, fb;
	double lim1, li, delim1, deli, delA;
	double lim2x, lim1x, lix, lip1x;
	double lim2y, lim1y, liy, lip1y;
	double dlim1_dxi, dlim1_dyi, dli_dxi, dli_dyi, dli_dxip1, dli_dyip1;
	double da_dxi, da_dyi, da_dxj, da_dyj;
	double kapim1, kapi, kapip1;
	double dkapi_dxi, dkapi_dyi, dkapip1_dxi, dkapip1_dyi, dkapim1_dxi, dkapim1_dyi;
	double dkapi_dxip1, dkapi_dyip1, dkapip1_dxip1, dkapip1_dyip1, dkapip1_dxip2, dkapip1_dyip2;
	double eij, kij, sij, dr, dx, dy, h;
	double dr_dxi, dr_dyi;

	// compute total number of degrees of freedom
	for (ci=0; ci<NCELLS; ci++){
		// save cumulative number of vertices
		if (ci == 0)		
			Mu.at(ci) = 0;
		else
			Mu.at(ci) = Mu.at(ci-1) + cell(ci-1).getNV();

		// total number of degrees of freedom
		NDOF += cell(ci).getNV()*NDIM;
	}


	// initialize matrices
	// NOTE: surface tension VDOS not yet supported
	Eigen::MatrixXd Ha(NDOF,NDOF);	// stiffness matrix for area term
	Eigen::MatrixXd Sa(NDOF,NDOF);	// stress matrix for area term
	Eigen::MatrixXd Hl(NDOF,NDOF);	// stiffness matrix for perimeter term
	Eigen::MatrixXd Sl(NDOF,NDOF);	// stress matrix for perimeter term
	Eigen::MatrixXd Hb(NDOF,NDOF);	// stiffness matrix for bending energy
	Eigen::MatrixXd Sb(NDOF,NDOF);	// stress matrix for bending term
	Eigen::MatrixXd Dvv(NDOF,NDOF);	// dynamical matrix for vertex-vertex interactions only
	Eigen::MatrixXd D(NDOF,NDOF);	// full dynamical matrix

	// initialize all matrices to be 0 initially
	for (k=0; k<NDOF; k++){
		for (l=0; l<NDOF; l++){
			Ha(k,l) = 0.0;
			Sa(k,l) = 0.0;
			Hl(k,l) = 0.0;
			Sl(k,l) = 0.0;
			Hb(k,l) = 0.0;
			Sb(k,l) = 0.0;
			Dvv(k,l) = 0.0;
			D(k,l) = 0.0;
		}
	}

	// compute initial forces to have upodate contact network
	calculateForces();

	// Loop over cells, compute shape forces for each individual cell and contributions from
	// vertex-vertex interactions
	for (ci=0; ci<NCELLS; ci++){

		// print statement
		cout << "	-- Computing dynamical matrix elements for cell ci = " << ci << endl;
		

		// ------------------------------------------
		// 
		// 				SHAPE
		// 					CONTRIBUTIONS
		//
		// ------------------------------------------

		// number of vertices of cell ci
		nv 			= cell(ci).getNV();

		// interaction constants
		calA0 		= cell(ci).calA0();
		l0 			= cell(ci).getl0();
		fl 			= cell(ci).getkl();
		kl 			= fl/l0;
		kb 			= cell(ci).getkb();
		eb 			= (kb*nv*calA0)/(4.0*PI*PI);
		fb 			= eb/(l0*l0);
		delA 		= (cell(ci).polygonArea()/cell(ci).geta0()) - 1.0;

		// loop over vertices, compute each DM element
		for (vi=0; vi<nv; vi++){

			// wrap vertices
			vip2 		= (vi + 2) % nv;
			vip1 		= (vi + 1) % nv;
			vim1 		= (vi - 1 + nv) % nv;
			vim2 		= (vi - 2 + nv) % nv;

			// dof elements
			kx 			= NDIM*(Mu.at(ci) + vi);
			ky 			= NDIM*(Mu.at(ci) + vi) + 1;

			kxp1 		= NDIM*(Mu.at(ci) + vip1);
			kyp1 		= NDIM*(Mu.at(ci) + vip1) + 1;

			kxp2 		= NDIM*(Mu.at(ci) + vip2);
			kyp2 		= NDIM*(Mu.at(ci) + vip2) + 1;

			// segment length vector components
			lim2x 		= cell(ci).segment(vim2,0);
			lim2y 		= cell(ci).segment(vim2,1);

			lim1x 		= cell(ci).segment(vim1,0);
			lim1y 		= cell(ci).segment(vim1,1);

			lix 		= cell(ci).segment(vi,0);
			liy 		= cell(ci).segment(vi,1);

			lip1x 		= cell(ci).segment(vip1,0);
			lip1y 		= cell(ci).segment(vip1,1);

			// segment lengths
			lim1 		= sqrt(lim1x*lim1x + lim1y*lim1y);
			li 			= sqrt(lix*lix + liy*liy);


			// -- PERIMETER SPRINGS

			// derivatives of lim1
			dlim1_dxi 	= lim1x/lim1;
			dlim1_dyi 	= lim1y/lim1;

			// derivatives of li
			dli_dxip1 	= lix/li;
			dli_dyip1 	= liy/li;

			dli_dxi 	= -dli_dxip1;
			dli_dyi 	= -dli_dyip1;

			// spring strains
			delim1 		= (lim1/l0) - 1.0;
			deli 		= (li/l0) - 1.0;

			// 	STIFFNESS MATRIX

			// main diagonal
		    Hl(kx,kx)       = kl*(dlim1_dxi*dlim1_dxi + dli_dxi*dli_dxi);
		    Hl(ky,ky)       = kl*(dlim1_dyi*dlim1_dyi + dli_dyi*dli_dyi);
		    
		    Hl(kx,ky)       = kl*(dlim1_dxi*dlim1_dyi + dli_dxi*dli_dyi);
		    Hl(ky,kx)       = Hl(kx,ky);
		    
		    // 1off diagonal
		    Hl(kx,kxp1)     = kl*dli_dxi*dli_dxip1;
		    Hl(ky,kyp1)     = kl*dli_dyi*dli_dyip1;
		    
		    Hl(kx,kyp1)     = kl*dli_dxi*dli_dyip1;
		    Hl(ky,kxp1)     = kl*dli_dyi*dli_dxip1;
		    
		    // enforce symmetry in lower triangle
		    Hl(kxp1,kx)     = Hl(kx,kxp1);
		    Hl(kyp1,ky)     = Hl(ky,kyp1);
		    
		    Hl(kyp1,kx)     = Hl(kx,kyp1);
		    Hl(kxp1,ky)     = Hl(ky,kxp1);


		    // 	STRESS MATRIX

		    // main diagonal block
		    Sl(kx,kx) 		= fl*( (delim1/lim1)*(1.0 - (dlim1_dxi*dlim1_dxi)) + (deli/li)*(1.0 - (dli_dxi*dli_dxi)) );
		    Sl(ky,ky) 		= fl*( (delim1/lim1)*(1.0 - (dlim1_dyi*dlim1_dyi)) + (deli/li)*(1.0 - (dli_dyi*dli_dyi)) );

		    Sl(kx,ky) 		= -fl*( (delim1/lim1)*dlim1_dxi*dlim1_dyi + (deli/li)*dli_dxi*dli_dyi );
		    Sl(ky,kx) 		= Sl(kx,ky);

		    // 1off diagonal
		    Sl(kx,kxp1) 	= fl*(deli/li)*((dli_dxip1*dli_dxip1) - 1.0);
		    Sl(ky,kyp1)		= fl*(deli/li)*((dli_dyip1*dli_dyip1) - 1.0);

		    Sl(kx,kyp1) 	= fl*(deli/li)*dli_dxip1*dli_dyip1;
		    Sl(ky,kxp1)		= fl*(deli/li)*dli_dyip1*dli_dxip1;

		    // enforce symmetry in lower triangle
		    Sl(kxp1,kx)     = Sl(kx,kxp1);
    		Sl(kyp1,ky)     = Sl(ky,kyp1);
    
    		Sl(kyp1,kx)     = Sl(kx,kyp1);
    		Sl(kxp1,ky)     = Sl(ky,kxp1);




    		// -- CURVATURE SPRINGS


    		// dimensionless curvatures
    		kapim1 			= sqrt(pow(lim1x - lim2x,2.0) + pow(lim1y - lim2y,2.0))/l0;
    		kapi 			= sqrt(pow(lix - lim1x,2.0) + pow(liy - lim1y,2.0))/l0;
    		kapip1 			= sqrt(pow(lip1x - lix,2.0) + pow(lip1y - liy,2.0))/l0;

    		// curvature derivatives

    		// derivatives of kapim1
    		dkapim1_dxi 	= (lim1x - lim2x)/(kapim1*l0*l0);
    		dkapim1_dyi 	= (lim1y - lim2y)/(kapim1*l0*l0);

    		// derivatives of kapi
    		dkapi_dxip1 	= (lix - lim1x)/(kapi*l0*l0);
    		dkapi_dyip1 	= (liy - lim1y)/(kapi*l0*l0);
    		dkapi_dxi 		= -2.0*dkapi_dxip1;
    		dkapi_dyi 		= -2.0*dkapi_dyip1;	

    		// derivatives of kapip1
    		dkapip1_dxi 	= (lip1x - lix)/(kapip1*l0*l0);
    		dkapip1_dyi 	= (lip1y - liy)/(kapip1*l0*l0);
    		dkapip1_dxip1 	= -2.0*dkapip1_dxi;
    		dkapip1_dyip1 	= -2.0*dkapip1_dyi;
    		dkapip1_dxip2	= dkapip1_dxi;
    		dkapip1_dyip2 	= dkapip1_dyi;


    		// 	STIFFNESS MATRIX

    		// block-diagonal terms
		    Hb(kx,kx)       = eb*(dkapim1_dxi*dkapim1_dxi + dkapi_dxi*dkapi_dxi + dkapip1_dxi*dkapip1_dxi);
		    Hb(ky,ky)       = eb*(dkapim1_dyi*dkapim1_dyi + dkapi_dyi*dkapi_dyi + dkapip1_dyi*dkapip1_dyi);
		    
		    Hb(kx,ky)       = eb*(dkapim1_dxi*dkapim1_dyi + dkapi_dxi*dkapi_dyi + dkapip1_dxi*dkapip1_dyi);
		    Hb(ky,kx)       = Hb(kx,ky);
		    
		    // 1off block-diagonal terms
		    Hb(kx,kxp1)     = eb*(dkapi_dxi*dkapi_dxip1 + dkapip1_dxi*dkapip1_dxip1);
		    Hb(ky,kyp1)     = eb*(dkapi_dyi*dkapi_dyip1 + dkapip1_dyi*dkapip1_dyip1);
		    
		    Hb(kx,kyp1)     = eb*(dkapi_dxi*dkapi_dyip1 + dkapip1_dxi*dkapip1_dyip1);
		    Hb(ky,kxp1)     = eb*(dkapi_dyi*dkapi_dxip1 + dkapip1_dyi*dkapip1_dxip1);
		    
		    // 2off block-diagonal terms
		    Hb(kx,kxp2)     = eb*dkapip1_dxi*dkapip1_dxip2;
		    Hb(ky,kyp2)     = eb*dkapip1_dyi*dkapip1_dyip2;
		    
		    Hb(kx,kyp2)     = eb*dkapip1_dxi*dkapip1_dyip2;
		    Hb(ky,kxp2)     = eb*dkapip1_dyi*dkapip1_dxip2;
		    
		    // enforce symmetry in lower triangle
		    Hb(kxp1,kx)     = Hb(kx,kxp1);
		    Hb(kyp1,ky)     = Hb(ky,kyp1);
		    
		    Hb(kxp1,ky)     = Hb(ky,kxp1);
		    Hb(kyp1,kx)     = Hb(kx,kyp1);
		    
		    Hb(kxp2,kx)     = Hb(kx,kxp2);
		    Hb(kyp2,ky)     = Hb(ky,kyp2);
		    
		    Hb(kyp2,kx)     = Hb(kx,kyp2);
		    Hb(kxp2,ky)     = Hb(ky,kxp2);
		    
		    
		    // 	STRESS MATRIX
		    
		    // block diagonal
		    Sb(kx,kx)       = fb*(6.0 - (l0*dkapim1_dxi)*(l0*dkapim1_dxi) - (l0*dkapi_dxi)*(l0*dkapi_dxi) - (l0*dkapip1_dxi)*(l0*dkapip1_dxi));
		    Sb(ky,ky)       = fb*(6.0 - (l0*dkapim1_dyi)*(l0*dkapim1_dyi) - (l0*dkapi_dyi)*(l0*dkapi_dyi) - (l0*dkapip1_dyi)*(l0*dkapip1_dyi));
		    
		    Sb(kx,ky)       = -eb*(dkapim1_dxi*dkapim1_dyi + dkapi_dxi*dkapi_dyi + dkapip1_dxi*dkapip1_dyi);
		    Sb(ky,kx)       = Sb(kx,ky);
		    
		    // 1off block diagonal
		    Sb(kx,kxp1)     = -2*fb*(2.0 - (l0*dkapi_dxip1)*(l0*dkapi_dxip1) - (l0*dkapip1_dxi)*(l0*dkapip1_dxi));
		    Sb(ky,kyp1)     = -2*fb*(2.0 - (l0*dkapi_dyip1)*(l0*dkapi_dyip1) - (l0*dkapip1_dyi)*(l0*dkapip1_dyi));
		    
		    Sb(kx,kyp1)     = -eb*(dkapi_dxi*dkapi_dyip1 + dkapip1_dxi*dkapip1_dyip1);
		    Sb(ky,kxp1)     = -eb*(dkapi_dyi*dkapi_dxip1 + dkapip1_dyi*dkapip1_dxip1);
		    
		    // 2off block diagonal
		    Sb(kx,kxp2)     = fb*(1.0 - (l0*dkapip1_dxi)*(l0*dkapip1_dxi));
		    Sb(ky,kyp2)     = fb*(1.0 - (l0*dkapip1_dyi)*(l0*dkapip1_dyi));
		    
		    Sb(kx,kyp2)     = -eb*dkapip1_dxi*dkapip1_dyip2;
		    Sb(ky,kxp2)     = -eb*dkapip1_dyi*dkapip1_dxip2;
		    
		    // enforce symmetry in lower triangle
		    Sb(kxp1,kx)     = Sb(kx,kxp1);
		    Sb(kyp1,ky)     = Sb(ky,kyp1);
		    
		    Sb(kxp1,ky)     = Sb(ky,kxp1);
		    Sb(kyp1,kx)     = Sb(kx,kyp1);
		    
		    Sb(kxp2,kx)     = Sb(kx,kxp2);
		    Sb(kyp2,ky)     = Sb(ky,kyp2);
		    
		    Sb(kxp2,ky)     = Sb(ky,kxp2);
		    Sb(kyp2,kx)     = Sb(kx,kyp2);


    		

		    // -- AREA SPRING (stress matrix)
		    Sa(kx,kyp1) = 0.5*delA;
    		Sa(ky,kxp1) = -0.5*delA;

    		Sa(kyp1,kx) = Sa(kx,kyp1);
    		Sa(kxp1,ky) = Sa(ky,kxp1);

    		// area derivatives (for stiffness matrix)
    		da_dxi      = 0.5*(cell(ci).vrel(vim1,1) - cell(ci).vrel(vip1,1));
    		da_dyi      = 0.5*(cell(ci).vrel(vip1,0) - cell(ci).vrel(vim1,0));

    		// loop over other vertices, for area elasticity stiffness matrix
    		for (vj=vi; vj<nv; vj++){

    			// wrap jp1 and jm1
    			vjp1 		= (vj + 1) % nv;
    			vjm1 		= (vj - 1 + nv) % nv;

    			// dof elements
    			lx 			= NDIM*(Mu.at(ci) + vj);
    			ly 			= NDIM*(Mu.at(ci) + vj) + 1;

    			// area derivatives
    			da_dxj      = 0.5*(cell(ci).vrel(vjm1,1) - cell(ci).vrel(vjp1,1));
    			da_dyj      = 0.5*(cell(ci).vrel(vjp1,0) - cell(ci).vrel(vjm1,0));

    			// 	STIFFNESS MATRIX
    			Ha(kx,lx) = da_dxi*da_dxj;
		        Ha(kx,ly) = da_dxi*da_dyj;
		        
		        Ha(ky,lx) = da_dyi*da_dxj;
		        Ha(ky,ly) = da_dyi*da_dyj;
		        
		        Ha(lx,kx) = Ha(kx,lx);
		        Ha(ly,kx) = Ha(kx,ly);
		        
		        Ha(lx,ky) = Ha(ky,lx);
		        Ha(ly,ky) = Ha(ky,ly);
    		}
		}






		// ------------------------------------------
		// 
		// 			INTERACTION
		// 					CONTRIBUTIONS
		//
		// ------------------------------------------

		// off-diagonal components
		for (cj=ci+1; cj<NCELLS; cj++){

			// interaction energy scale
			eij = 0.5*(cell(ci).getkint() + cell(cj).getkint());

			// loop over pairs of vertices on both cells, check for overlap, compute matrix elements
			for (vi=0; vi<nv; vi++){

				// matrix element indices (cell ci, vertex vi)
				mxi = NDIM*(Mu.at(ci) + vi);
				myi = NDIM*(Mu.at(ci) + vi) + 1;

				// contact distance
				sij = 0.5*(cell(ci).getdel()*l0 + cell(cj).getdel()*cell(cj).getl0());

				for (vj=0; vj<cell(cj).getNV(); vj++){

					// get distance between vertices
					dx = cell(ci).distance(cell(cj),vj,vi,0);
					dy = cell(ci).distance(cell(cj),vj,vi,1);
					dr = sqrt(dx*dx + dy*dy);

					// check for overlap
					if (dr < sij){

						// spring constant
						kij = eij/(sij*dr);

						// dimensionless overlap
						h = dr/sij;

						// matrix element indices (cell cj, vertex vj)
						mxj = NDIM*(Mu.at(cj) + vj);
						myj = NDIM*(Mu.at(cj) + vj) + 1;

						// derivatives of distance w.r.t. coordinates
						dr_dxi = -dx/dr;
						dr_dyi = -dy/dr;

						// set off diagonals, enforce symmetry in lower triangle
						Dvv(mxi,mxj) = -kij*(dr_dxi*dr_dxi + h - 1.0);
		                Dvv(myi,myj) = -kij*(dr_dyi*dr_dyi + h - 1.0);
		                Dvv(mxi,myj) = -kij*dr_dxi*dr_dyi;
		                Dvv(myi,mxj) = -kij*dr_dxi*dr_dyi;
		                
		                Dvv(mxj,mxi) = Dvv(mxi,mxj);
		                Dvv(myj,myi) = Dvv(myi,myj);
		                Dvv(mxj,myi) = Dvv(myi,mxj);
		                Dvv(myj,mxi) = Dvv(mxi,myj);
		                
		                // add to diagonal, using off diagonals and reciprocity
		                Dvv(mxi,mxi) -= Dvv(mxi,mxj);
		                Dvv(myi,myi) -= Dvv(myi,myj);
		                Dvv(mxi,myi) -= Dvv(mxi,myj);
		                Dvv(myi,mxi) -= Dvv(myi,mxj);
		                
		                Dvv(mxj,mxj) -= Dvv(mxi,mxj);
		                Dvv(myj,myj) -= Dvv(myi,myj);
		                Dvv(mxj,myj) -= Dvv(mxi,myj);
		                Dvv(myj,mxj) -= Dvv(myi,mxj);
					}
				}
			}
		}
	}

	// compute D from sum of other dynamical matrices
	// initialize all matrices to be 0 initially
	for (k=0; k<NDOF; k++){
		for (l=0; l<NDOF; l++)
			D(k,l) = Ha(k,l) + Sa(k,l) + Hl(k,l) + Sl(k,l) + Hb(k,l) + Sb(k,l) + Dvv(k,l);
	}

	// compute eigenvalues
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> normalModes(D);

	// print eigenvalues to print object
	statPrintObject << NDOF << endl;
	statPrintObject << normalModes.eigenvalues() << endl;
	statPrintObject << normalModes.eigenvectors() << endl;
	Eigen::MatrixXd evecs = normalModes.eigenvectors();

	// print to console
	double U0 = totalPotentialEnergy();
	cout << "	** Finished printing evals and evecs, now printing change from U0 = " << U0 << " along evecs" << endl;

	// loop over eigenvectors, perturb by small increments along eigenvector, compute change in potential energy
	int NSTEPS = 30;
	double p0 = 1e-6;
	double p1 = 1e0;
	vector<double> deList(NSTEPS,0.0);
	double dp = (log10(p1) - log10(p0))/(NSTEPS - 1);
	double logp, linp;
	logp = log10(p0);
	deList.at(0) = p0;

	// print number of steps
	statPrintObject << NSTEPS << endl;
	statPrintObject << setw(30) << setprecision(12) << U0 << endl;

	// loop over vector, populate with points
	statPrintObject << setw(30) << setprecision(12) << deList.at(0);
	for (k=1; k<NSTEPS; k++){
		// add to vector
		logp = logp + dp;
		linp = pow(10.0,logp);
		deList.at(k) = linp;

		// print vector value to file
		statPrintObject << setw(30) << setprecision(12) << deList.at(k);
	}
	statPrintObject << endl;

	// store packing in relaxed state
	saveState(relaxedState);

	// loop over modes
	int m, s, d;
	double de, ptmp;
	for (m=0; m<NDOF; m++){
		// loop over steps
		for (s=0; s<NSTEPS; s++){
			// load step size
			de = deList.at(s);

			// move system along eigenvector (note: always start from jammed packing)
			k = 0;
			for (ci=0; ci<NCELLS; ci++){
				for (vi=0; vi<cell(ci).getNV(); vi++){
					for (d=0; d<NDIM; d++){
						// initial position
						ptmp = cell(ci).vpos(vi,d);

						// increment from initial jammed state to perturbed state
						ptmp += de*evecs(k,m);
						k++;

						// store new position
						cell(ci).setVPos(vi,d,ptmp);
					}
				}
			}

			// output energy after all degrees of freedom have been perturbed
			statPrintObject << setw(30) << setprecision(12) << totalPotentialEnergy();

			// load relaxed state for next iteration (and next compression)
			loadState(relaxedState);
		}
		statPrintObject << endl;
	}
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
void cellPacking2D::initializeGel(int NV, double phiDisk, double sizeDispersion, double delval){
	// local variables
	int ci, vi, d, nvtmp;
	int lx, ly;
	double calA0tmp;
	double xpos, ypos, dx, dy;
	double xmin, xmax, ymin, ymax;
	double r1, r2, g1, areaSum;
	double rtmp, a0tmp, l0tmp, p0tmp;

	// seed random number generator
	srand48(23562457*seed);

	// minimum number of vertices
	const int nvmin = 12;

	// output to console
	cout << "		-- In gelation initialization, initializing cells and relaxing initial overlaps as repulsive SP particles" << endl;

	// initialize length scales as gaussian random variables (will becomes area square roots)
	areaSum = 0.0;
	vector<double> lenscales(NCELLS,0.0);
	vector<double> diskradii(NCELLS,0.0);
	for (ci=0; ci<NCELLS; ci++){
		// generate random numbers
		r1 = drand48();
		r2 = drand48();

		// calculate gaussian random variable using Box-Muller transform
		g1 = sqrt(-2.0*log(r1))*cos(2*PI*r2);

		// get root area
		lenscales.at(ci) = g1*sizeDispersion + 1.0;

		// effective particle area
		a0tmp = lenscales.at(ci)*lenscales.at(ci);

		// add to lenscales sum for boundary size
		areaSum += a0tmp;

		// save disk radius based on area square root
		diskradii.at(ci) = sqrt(a0tmp/PI);
	}

	// determine box length from particle sizes and input packing fraction
	for (d=0; d<NDIM; d++)
		L.at(d) = sqrt(areaSum/phiDisk);

	// set phi to input
	phi = phiDisk;

	// reseed rng
	srand48(56835698*seed);

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

		// initialize cells as regular polygons
		calA0tmp = nvtmp*tan(PI/NV)/PI;

		// preferred area is lenscales squared
		a0tmp = lenscales.at(ci)*lenscales.at(ci);

		// initial length of polygon side
		l0tmp = 2.0*lenscales.at(ci)*sqrt(PI*calA0tmp)/NV;

		// set preferred area and length 
		cell(ci).seta0(a0tmp);
		cell(ci).setl0(l0tmp);
		cell(ci).setdel(delval);
	}

	// initialize particle positions
	cout << "		-- Ininitializing cell positions" << endl;

	// set min and max values of positions
	xmin = 0;
	xmax = L.at(0);
	ymin = 0;
	ymax = L.at(1);

	// initialize positions of each cell
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

	// set time step
	vertexDPMTimeScale(0.1);

	// use FIRE in PBC box to relax overlaps
	cout << "		-- Using FIRE to relax overlaps..." << endl;
	fireMinimizeSP(diskradii);
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


// rk4 update
void cellPacking2D::cellRK4(){
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



// CANCER SIMULATION FUNCTIONS

// function to simulate a single active tumor cell (particle 0) 
// void cellPacking2D::singleActiveTumorCell(double tumorCalA0, double tv0, double tDr)





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
	double Ftol 			= 1e-10;
	double Ktol 			= 1e-20;
	bool converged 			= false;

	// local variables
	int ci,vi,d,itr,itrMax;
	double P,vstarnrm,fstarnrm,vtmp,ftmp,ptmp;
	double F, K, Fcheck, Kcheck, Pvirial;

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

	// update force RMS
	F = forceRMS();

	// update kinetic energy based on com velocity
	K = 0.0;
	for (ci=0; ci<NCELLS; ci++)
		K += 0.5*(PI*pow(lenscales.at(ci),2))*sqrt(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2));

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
			cout << "	* K 		= " << K << endl;
			cout << "	* F 		= " << F << endl;
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

		// calculate forces between disks
		spForces(lenscales);

		// verlet velocity update
		spVelVerlet(lenscales);

		// update t
		t += dt;

		// track energy and forces
		F = forceRMS();
		K = totalKineticEnergy();

		// scale P and K for convergence checking
		Fcheck = F;
		Kcheck = K/NCELLS;

		// update if Pvirial under tol
		if (Fcheck < Ftol)
			npPMIN++;
		else
			npPMIN = 0;

		// check for convergence
		converged = (Fcheck < Ftol && npPMIN > NMIN && Kcheck < Ktol);

		if (converged){
			cout << "	** FIRE has converged!" << endl;
			cout << "	** Kcheck of Sp particles = " << Kcheck << endl;
			cout << "	** Fcheck of Sp particles = " << Fcheck << endl;
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
		for (d=0; d<NDIM; d++){
			// update new position based on acceleration
			postmp = cell(ci).cpos(d) + dt*cell(ci).cvel(d) + 0.5*dt*dt*cell(ci).cforce(d);

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
	double veltmp, aold, anew;

	// update com velocity
	for (ci=0; ci<NCELLS; ci++){
		// loop over velocities
		for (d=0; d<NDIM; d++){
			// get current velocity
			veltmp = cell(ci).cvel(d);

			// calculate old com acceleration
			aold = 0.0;
			for (vi=0; vi<cell(ci).getNV(); vi++)
				aold += cell(ci).vacc(vi,d);

			// get new accelation
			anew = cell(ci).cforce(d);

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
	jamPrintObject << setw(w1) << left << "NCONTS";
	jamPrintObject << setw(w1) << right << Ncc;
	jamPrintObject << setw(w1) << right << Nvv;
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

	// check to see if file is open
	if (!energyPrintObject.is_open()) {
		cout << "	ERROR: energyPrintObject is not open in printSystemEnergy(), ending." << endl;
		exit(1);
	}

	// loop over particles, print cell energy
	energyPrintObject << setw(30) << setprecision(16) << right << interactionPotentialEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << totalPotentialEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << totalKineticEnergy();
	energyPrintObject << setw(30) << setprecision(16) << right << forceRMS();
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



