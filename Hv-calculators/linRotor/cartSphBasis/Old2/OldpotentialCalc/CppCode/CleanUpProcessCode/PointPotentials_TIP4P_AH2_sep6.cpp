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

/* To run:
 ./PointPotentials_TIP4P_AH2 [dim] [max(1) num(1) max(2) num(2) ... max(dim) num(dim)] [phi] [theta] [inputFilename] [OutputFilename]
 */
int main (int argc, char** argv)
{	
	
	//Universe Grid variables
	int *box_i_num, numDim;
	int systemSize;
	double *box_i_max, *d_i;
	VECT *universe_grid;
	
	
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
	
	//OpenMP Variables
	int thread, nthreads;
	int chunk;
	
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
	d_i = new double [numDim]; //ODDLY, this seems to have given a value of 0 once, though that was before the memory problems with origin and atompos were fixed.
	
	//Gather in command line universe dimension information
	systemSize = 1;
	for (i=0; i<numDim; i++) {
		box_i_max[i] = atof(argv[2*i+2]);
		box_i_num[i] = atoi(argv[2*i+3]);
		d_i[i] = (box_i_max[i]) / ((double) box_i_num[i]);
		systemSize *= box_i_num[i];
	}
	
	//Initialize grids
	universe_grid = new VECT [systemSize];
	CMpotential = new double [systemSize];
	Hpotential = new double [systemSize];
	H2potential = new double [systemSize];
	
	//Get orientation from command line
	phi = atof(argv[2*numDim+2]);
	theta = atof(argv[2*numDim+2+1]);
	
	H2_mol_orient.phi = phi;
	H2_mol_orient.theta = theta;
	
	//Determine for-loop chunk for OpenMP
	if ((systemSize*numDim) > CHUNK_RATIO) {
		chunk = (int) (systemSize/CHUNK_RATIO);
	}
	else if ((systemSize*numDim) > (CHUNK_RATIO/10)){
		chunk = (int) (systemSize/(CHUNK_RATIO/10));
	}
	else {
		chunk = 1;
	}
	
	cout << "Program Initialized." << endl;
	
	//Generate universe grid ###############################################################
#pragma omp parallel default(shared) private (i,j,thread)
	{
	thread = omp_get_thread_num();
	if (thread == 0) {
		nthreads = omp_get_num_threads();
		cout << "Grid Generation Parallel Section Number of Threads: " << nthreads << endl;
	}
	
#pragma omp for schedule(dynamic, chunk)  //Possibly try: http://software.intel.com/en-us/articles/openmp-loop-collapse-directive
	for (i=0; i<systemSize; i++) {
		universe_grid[i].DIM(numDim);
		for (j=0; j<numDim; j++) {
			//universe_grid[i].COOR(j) = indexDim(j, i, numDim, box_i_num)*d_i[j]; //Seems to work, but not fully tested.
			universe_grid[i].COOR(j) = (indexDim(j, i, numDim, box_i_num)*d_i[j]) - (box_i_max[j]/2.0); //Seems to work, but not fully tested.; subtract 1/2*box_i_max to centre origin
		}		
	}
	}
	
	cout << "Grid generated." << endl;
	
	
	//Store information for the point universe (ie. CM and H atom locations within H2)
	point_universe.grid = universe_grid;
	point_universe.grid_max = box_i_max;
	point_universe.d_i = d_i;
	point_universe.grid_num = box_i_num;
	point_universe.numDim = numDim;
	point_universe.sysSize = systemSize;
	
	//Store information for the molecule universe (which may have a different grid than for the point universe)
	H2_universe.grid = universe_grid;
	H2_universe.grid_max = box_i_max;
	H2_universe.d_i = d_i;
	H2_universe.grid_num = box_i_num;
	H2_universe.numDim = numDim;
	H2_universe.sysSize = systemSize;
	
	
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
	dimGenerator = new VECT [numDim];
	for (i=0; i<numDim; i++) {
		dimGenerator[i].DIM(numDim);
		
		for (j=0; j<numDim; j++) {
			dimGenerator[i].COOR(j) = 0.0;
		}
		
		dimGenerator[i].COOR(i) = d_i[i]/NM_PER_BOHR; //Convert to Bohr for Cube file
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
	origin.DIM(numDim);
	for (j=0; j<numDim; j++) {
		origin.COOR(j) = -1.0*box_i_max[j]/2/NM_PER_BOHR; //Convert to Bohr for Cube file
	}
	
	for (i=0; i<nAtoms; i++) {
		for (j=0; j<numDim; j++) {
			atomPos[i].COOR(j) /= NM_PER_BOHR; //Convert to Bohr for Cube file
		}
	}
	
	unitFlag = -1;
	
	comments[0] = "Comment 1";
	comments[1] = "Comment 2";
	
	filename = argv[(numDim*2)+5];	
	
	filenameCM = filename + "_CM.cube";
	
	writeCubeFile(numDim, box_i_num, dimGenerator, CMpotential, systemSize, 
				  nAtoms, atomNum, atomCharges, atomPos, origin, comments, 
				  unitFlag, filenameCM);
	
	cout << "Centre of Mass data written to file '" << filenameCM << "'." << endl;
	
	filenameH = filename + "_H.cube";
	
	writeCubeFile(numDim, box_i_num, dimGenerator, Hpotential, systemSize, 
				  nAtoms, atomNum, atomCharges, atomPos, origin, comments, 
				  unitFlag, filenameH);
	
	cout << "H-atom data written to file '" << filenameH << "'." << endl;
	
	
	filenameH2 = filename + "_H2.cube";
	
	writeCubeFile(numDim, box_i_num, dimGenerator, H2potential, systemSize, 
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
	delete [] d_i;
	delete [] universe_grid;
	delete [] CMpotential;
	delete [] Hpotential;
	delete [] H2potential;
	
	return 0;
}
