#include <iostream>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <string.h>
#include "gridFcns.h"
#include "energyRoutines.h"
#include "vectClass.h" //This is to select which linear algebra class should be used.
#include "TIP4P_AH2_EngRoutines.h"
#include "rotmat.h"
#include "potentialUnits.h"
#include "boundStateContainers.h"

using namespace std;

universeProp generateGrid(int numDim, double *gridMax, int *gridPoints) {
	universeProp universe;
	int chunk, thread, nthreads, i, j;
	
	universe.numDim = numDim;
	universe.grid_max = gridMax;
	universe.grid_num = gridPoints;
	
	universe.d_i = new double [universe.numDim];
	
	universe.sysSize = 1;
	for (i=0; i<universe.numDim; i++) {
		universe.d_i[i] = (universe.grid_max[i]) / ((double) universe.grid_num[i]);
		universe.sysSize *= universe.grid_num[i];
	}
	
	universe.grid = new VECT [universe.sysSize];
	
	
	//Determine for-loop chunk for OpenMP
	if (universe.sysSize > CHUNK_RATIO) {
		chunk = (int) (universe.sysSize/CHUNK_RATIO);
	}
	else if (universe.sysSize > (CHUNK_RATIO/10)){
		chunk = (int) (universe.sysSize/(CHUNK_RATIO/10));
	}
	else {
		chunk = 1;
	}
	
	//Generate universe grid ###############################################################
#pragma omp parallel default(shared) private (i,j,thread)
	{
	thread = omp_get_thread_num();
	if (thread == 0) {
		nthreads = omp_get_num_threads();
		cout << "Grid Generation Parallel Section Number of Threads: " << nthreads << endl;
	}
	
#pragma omp for schedule(dynamic, chunk) 
	for (i=0; i<universe.sysSize; i++) {
		universe.grid[i].DIM(universe.numDim);
		for (j=0; j<universe.numDim; j++) {
			//universe.grid[i].COOR(j) = indexDim(j, i, universe.numDim, universe.grid_num)*universe.d_i[j]; //Seems to work, but not fully tested.
			universe.grid[i].COOR(j) = (indexDim(j, i, universe.numDim, universe.grid_num)*universe.d_i[j]) - (universe.grid_max[j]/2.0); //Seems to work, but not fully tested.; subtract 1/2*box_i_max to centre origin
		}		
	}
	}
	
	cout << "Grid generated." << endl;
	
	return universe;
}
	
	

/* To run:
 ./PointPotentials_TIP4P_AH2 [dim] [max(1) num(1) max(2) num(2) ... max(dim) num(dim)] [phi] [theta] [inputFilename] [OutputFilename]
 */
int main (int argc, char** argv)
{	
	
	//Universe Grid variables
	int *box_i_num, numDim;
	double *box_i_max;	
	
	//Alavi potential variables
	double *CMpotential, *Hpotential, *H2potential;
	
	//Gaussian Cube File Variables
	VECT *dimGenerator, origin;
	int *atomNum, unitFlag;
	double *atomCharges;
	string comments[2];
	string filename, filenameCM, filenameH, filenameH2;
	
	//System Geometry variables
	char *atomType;
	VECT *atomPos, pos;
	int nAtoms;
	
	int i, j;
	
	//H2 orientation variables
	double phi, theta;
	
	//Storage containers
	H2_orient H2_mol_orient;
	universeProp point_universe, H2_universe;
	
	//Check number of command line arguments
	if (argc <= (0+1)) {
		cout << "Invalid number of arguments, need at least one." << endl;
		exit (1);
	}
	
	// Get universe size ###############################################################
	
	//Get number of dimensions
	numDim = atoi(argv[1]);
	
	if (argc != (1+(numDim*2)+1+1+1+2)) { //runCommand, (dim_max dim_size), numDim, filename_in, filename_out, angles
		cout << "Invalid number of arguments, need " << (numDim*2)+1+1+1; 
		cout << " arguments for " << numDim << " dimensions." << endl;
		exit (1);
	}
	
	//Get dimension sizes
	box_i_max = new double [numDim];
	box_i_num = new int [numDim];
	
	//Gather in command line universe dimension information
	for (i=0; i<numDim; i++) {
		box_i_max[i] = atof(argv[2*i+2]);
		box_i_num[i] = atoi(argv[2*i+3]);
	}
	
	//Get orientation from command line
	phi = atof(argv[2*numDim+2]);
	theta = atof(argv[2*numDim+2+1]);
	
	H2_mol_orient.phi = phi;
	H2_mol_orient.theta = theta;
	
	cout << "Program Initialized." << endl;
	
	//Generate point universe grid ###############################################################
	
	point_universe = generateGrid(numDim, box_i_max, box_i_num);
	
	//Store information for the molecule universe (which may have a different grid than for the point universe, though it is the same in this case)
	H2_universe.grid = point_universe.grid;
	H2_universe.grid_max = point_universe.grid_max;
	H2_universe.d_i = point_universe.d_i;
	H2_universe.grid_num = point_universe.grid_num;
	H2_universe.numDim = point_universe.numDim;
	H2_universe.sysSize = point_universe.sysSize;
	
	//Initialize potentials
	CMpotential = new double [point_universe.sysSize];
	Hpotential = new double [point_universe.sysSize];
	H2potential = new double [point_universe.sysSize];
	
	//Acquire system geometry ###############################################################
	filename = argv[(numDim*2)+4];
	
	getTIP4Patoms(&atomType, &atomPos, &nAtoms, filename);
	
	
	//Store geometry information
	sysAtoms atomGeo;
	
	atomGeo.atomPos = atomPos;
	atomGeo.atomType = atomType;
	atomGeo.nAtoms = nAtoms;
	
	//Calculate energy for centre of mass and H atom ###############################################################
	Alavi_point_Eng(CMpotential, Hpotential, H2potential, &point_universe, &atomGeo);
	
	//Calculate Energy for H2 molecule ###############################################################
	Alavi_H2_Eng(H2potential, CMpotential, Hpotential, &H2_mol_orient, &point_universe);
	
	//##############################################################################################################################
	//##############################################################################################################################
	/*
	 //Write the data to a cube file ###############################################################
	 */
	//Make generators for each dimension
	dimGenerator = new VECT [point_universe.numDim];
	for (i=0; i<point_universe.numDim; i++) {
		dimGenerator[i].DIM(point_universe.numDim);
		
		for (j=0; j<point_universe.numDim; j++) {
			dimGenerator[i].COOR(j) = 0.0;
		}
		
		dimGenerator[i].COOR(i) = point_universe.d_i[i]/NM_PER_BOHR; //Convert to Bohr for Cube file
	}
	
	//Set up atom numbers
	atomNum = new int [nAtoms];
	
	for (i=0; i<nAtoms; i++) {
		switch (atomType[i]) {
			case 'H':
				atomNum[i] = 1;
				break;
			case 'M':
				atomNum[i] = 18; //Set to Argon
				break;
			case 'O':
				atomNum[i] = 8;
				break;
			default:
				cerr << "Atom type not recognized while determining atom numbers for Cube file." << endl;
				exit(1);
				break;
		}
	}
	
	//Set up atom charges (set all to zero)
	atomCharges = new double [nAtoms];
	
	for (i=0; i<nAtoms; i++) {
		atomCharges[i] = 0.0;
	}
	
	// Generate Origin Location, which is the centre of the universe
	origin.DIM(point_universe.numDim);
	for (j=0; j<point_universe.numDim; j++) {
		origin.COOR(j) = -1.0*point_universe.grid_max[j]/2/NM_PER_BOHR; //Convert to Bohr for Cube file
	}
	
	for (i=0; i<nAtoms; i++) {
		for (j=0; j<point_universe.numDim; j++) {
			atomPos[i].COOR(j) /= NM_PER_BOHR; //Convert to Bohr for Cube file
		}
	}
	
	unitFlag = +1;
	
	comments[0] = "Comment 1";
	comments[1] = "Comment 2";
	
	filename = argv[(numDim*2)+5];	
	
	filenameCM = filename + "_CM.cube";
	
	writeCubeFile(point_universe.numDim, point_universe.grid_num, dimGenerator, CMpotential, point_universe.sysSize, 
				  nAtoms, atomNum, atomCharges, atomPos, origin, comments, 
				  unitFlag, filenameCM);
	
	cout << "Centre of Mass data written to file '" << filenameCM << "'." << endl;
	
	filenameH = filename + "_H.cube";
	
	writeCubeFile(point_universe.numDim, point_universe.grid_num, dimGenerator, Hpotential, point_universe.sysSize, 
				  nAtoms, atomNum, atomCharges, atomPos, origin, comments, 
				  unitFlag, filenameH);
	
	cout << "H-atom data written to file '" << filenameH << "'." << endl;
	
	
	filenameH2 = filename + "_H2.cube";
	
	writeCubeFile(H2_universe.numDim, H2_universe.grid_num, dimGenerator, H2potential, H2_universe.sysSize, 
				  nAtoms, atomNum, atomCharges, atomPos, origin, comments, 
				  unitFlag, filenameH2);
	
	cout << "H2 data written to file '" << filenameH2 << "'." << endl;
	
	
	
	//Free memory
	delete [] atomNum;
	delete [] atomCharges;
	delete [] atomPos;	
	delete [] atomType;
	delete [] dimGenerator;
	delete [] box_i_max;
	delete [] box_i_num;
	delete [] CMpotential;
	delete [] Hpotential;
	delete [] H2potential;
	
	return 0;
}
