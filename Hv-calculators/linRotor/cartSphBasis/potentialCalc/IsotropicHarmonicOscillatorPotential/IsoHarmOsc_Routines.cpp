#include "IsoHarmOsc_Routines.h"

using namespace std;

//This calculates the potential energy of an Alavi Hydrogen at a particular centre of mass position and orientation
double IsoHarmOscEng(interfaceStor *interface, H2_orient *lin_mol) {
	
	double dist2, potential, factor, V_ceil;
	VECT distVec;
	
	distVec.DIM(3);
	
	distVec = *(lin_mol->CM);
	
	dist2 = distVec*distVec;
	
	factor = interface->potential->molCharge;
	V_ceil = interface->potential->coulombPotentialCeil;
	
	potential = factor * dist2;
	
	//cout << "HamrOscVal: " << factor << "; " << dist2 << "; " << potential<< "; " << V_ceil << endl;
	
	//potential =0.0;
	
	if (fabs(potential) > fabs(V_ceil)) {
		potential = V_ceil;
	} 
	
	return potential;
}


pointPotentialStorH2* preCalcPotential_IsoHarmOsc(int numDim, double *gridMax, int *gridPoints, string geometryFilename, interfaceStor *interface, int argc, char **argv) {
	
	
	double frequency, mass, factor, IsoHarmOscCeil;
	string line, junk, inputFilename;
	ifstream inputFile;
	
	inputFilename = argv[3];
	
	inputFile.open(inputFilename.c_str(), ios::in);
	if (inputFile.is_open()) {
		//Get rid of comment lines;
		getline(inputFile, line);
		getline(inputFile, line);
		
		//Gather data for input; each line has the name of the value separate from the value with a space
		//Always ignore the name of the value.
		inputFile >> junk;
		inputFile >> mass;		
		
		inputFile >> junk;
		inputFile >> frequency;
		
		inputFile >> junk;
		inputFile >> IsoHarmOscCeil;
		
	}
	else {
		cerr << "Isotropic Harmonic Oscillator Potential input file '" << inputFilename << "' could not be opened." << endl;
		exit(1);
	}
	
	cout << "Isotropic Harmonic Oscillator Potential Prepared." << endl;
	
	cout << "----------------------------------------------------------------------------" << endl;
	
	//Store everything for interface
	pointPotentialStorH2 *partialPotential = new pointPotentialStorH2();
	
	factor = 0.5 * mass * frequency * frequency;

	partialPotential->molCharge = factor;
	partialPotential->coulombPotentialCeil = IsoHarmOscCeil;
	
	return partialPotential;
}
