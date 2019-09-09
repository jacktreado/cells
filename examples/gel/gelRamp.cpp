/*

	Test .cpp file to test DP extensional
	simulation

*/

// include file
#include "cellPacking2D.h"
#include "deformableParticles2D.h"

// use std name space
using namespace std;

// define PI
const double PI = 4.0*atan(1);

// simulation constants
const int NT 					= 1e6;			// number of time steps
const int NPRINT 				= 500;			// number of time steps between prints
const double T0 				= 0.0;			// temperature scale
const double timeStepMag 		= 0.01;		// time step in MD units
const double dampingParameter 	= 0.9;			// damping parameter
const double phiTarget 			= 0.6;			// target packing fraction
const double dphi				= 0.001;		// decrease in packing fraction increment

// force parameters
const double kl 			= 1.0;			// perimeter force constant
const double ka 			= 1.0;			// area force constant
const double gam 			= 0.0;			// surface tension force constant
const double kb 			= 0.0;			// bending energy constant
const double kint 			= 5.0;			// interaction energy constant
const double del 			= 0.1;			// width of circulo lines
const double C 				= 0.01;			// attraction parameter (strength)
const double l 				= 2.0;			// attraction parameter (distance)

int main(){
	// output file strings
	string positionsStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/gel/gelRampPositions.dat";
	string energyStr = "/Users/JackTreado/Jamming/Flowers/sim/cells/examples/gel/gelRampEnergy.dat";

	// input simulation params
	double asphericity 	= 1.09;

	// print initial greeting
	cout << "	** BEGINNING main code for gel simulation, initial asphericity = " << asphericity << endl;

	// open input object
	cout << "	** Opening file object from string " << inputStr << endl;
	ifstream fileInputObject(inputStr.c_str());
	if (!fileInputObject.is_open()){
		cout << "	** ERROR: could not open file " << inputStr << ", ending." << endl;
		return 1;
	}

	// instantiate object
	cout << "	** Instantiating object for gel simulation" << endl;
	cellPacking2D gelObject(fileInputObject, asphericity, 1);

	/***************************************

		----- RAMP PROTOCOL -----

			1. Ramp from initial shapes to desired asphericity
			with step size dCalA
	
			2. Take shapes from initial packing fraction phiInitial
			to final phiFinal (if phiFinal = -1, pack until jammed)
			with step size dPhiGrow

			3. Add in attraction slowly using ramp, from initially no
			attraction to desired attraction with C and l parameters,
			with step sizes dC and dl

			4. Quasistatically decrease packing fraction to phiTarget using 
			step size dPhiShrink

		----- RAMP PROTOCOL -----

	*****************************************/




	// end main successfully
	return 0;
}











