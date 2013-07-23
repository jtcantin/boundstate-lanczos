#include "Coulomb_Routines.h"

using namespace std;

//This calculates the potential energy of an Alavi Hydrogen at a particular centre of mass position and orientation
double CoulombPotential(interfaceStor *interface, H2_orient *lin_mol) {
	
	double potential, q_centre, q_mol, V_ceil;
	VECT zeroVec;
	
	potential =0.0;
	/*
	zeroVec.DIM(3);
	zeroVec.COOR(0) = 0.0;
	zeroVec.COOR(1) = 0.0;
	zeroVec.COOR(2) = 0.0;
	
	q_centre = interface->potential->centreCharge;
	q_mol = interface->potential->molCharge;
	V_ceil = interface->potential->coulombPotentialCeil;
	
	potential = CoulombEng(zeroVec, q_centre, *(lin_mol->CM), q_mol);
	
	if (fabs(potential) > fabs(V_ceil)) {
		potential = V_ceil;
	} */
	
	return potential;
}


//The following function precalculates the Alavi hydrogen site potentials by setting up the PES grid, getting the system geometry, 
//   and calculating the CM and H-atom potentials
pointPotentialStorH2* preCalcPotential_Coulomb(int numDim, double *gridMax, int *gridPoints, string geometryFilename, interfaceStor *interface, int argc, char **argv) {
	
	//There are 3 spatial dimensions, so numDim = 3 and the gridMax and gridPoints arrays are of length 3
	
	//Get the grid for the potential
	//universeProp *potentialUniverse;
	
	//potentialUniverse = generateGrid(numDim, gridMax, gridPoints);
	
	double centreCharge, molCharge, coulombPotentialCeil;
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
		inputFile >> centreCharge;		
		
		inputFile >> junk;
		inputFile >> molCharge;
		
		inputFile >> junk;
		inputFile >> coulombPotentialCeil;
		
	}
	else {
		cerr << "Coulomb Potential input file '" << inputFilename << "' could not be opened." << endl;
		exit(1);
	}
	
	cout << "Coulomb Potential Prepared." << endl;
	
	cout << "----------------------------------------------------------------------------" << endl;
	
	//Store everything for interface
	pointPotentialStorH2 *partialPotential = new pointPotentialStorH2();

	partialPotential->centreCharge = centreCharge;
	partialPotential->molCharge = molCharge;
	partialPotential->coulombPotentialCeil = coulombPotentialCeil;
	
	return partialPotential;
}
