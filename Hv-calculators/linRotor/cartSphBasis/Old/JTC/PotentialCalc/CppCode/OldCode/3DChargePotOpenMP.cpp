#include <iostream>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <string.h>
#include "gridFcns.h"
#include "energyRoutines.h"
#include "vectClass.h" //This is to select which linear algebra class should be used.

using namespace std;

#define CHUNK_RATIO 100 //How large the chunks should be relative to the total data size; 10 or 100 seem to be reasonable values.

/* To run:
 mpirun -np N 3DChargePot [dim] [max(1) num(1) max(2) num(2) ... max(dim) num(dim) filename]
 */
int main (int argc, char** argv)
{	
	
	int *box_i_num, dim_num;
	int systemSize;
	double *box_i_max, *d_i;
	double *potential, *coor_temp;
	VECT *universe_grid;
	int *universe_grid_index;
	
	VECT *dimGenerator, origin, *atomPos;
	int nAtoms, *atomNum, unitFlag;
	double *atomCharges;
	char **comments, *filename;
	
	int i, j;
	
	VECT sourceCharge, testCharge;
	double sourceQ, testQ;
	
	//OpenMP Variables ###############################################################	
	int thread, nthreads;
	int chunk;
	
	if (argc <= (0+1)) {
		cout << "Invalid number of arguments, need at least one." << endl;
		exit (1);
	}
	
	// Get universe size ###############################################################
	
	//Get number of dimensions
	dim_num = atoi(argv[1]);
	
	if (argc != ((dim_num*2)+2+1)) { //dim_max dim_size, dim_num+runCommand, filename
		cout << "Invalid number of arguments, need " << (dim_num*2)+1; 
		cout << "arguments for " << dim_num << " dimensions." << endl;
		exit (1);
	}
	
	//Get dimension sizes
	box_i_max = new double [dim_num];
	box_i_num = new int [dim_num];
	d_i = new double [dim_num];
	systemSize = 1;
	for (i=0; i<dim_num; i++) {
		box_i_max[i] = atof(argv[2*i+2]);
		box_i_num[i] = atoi(argv[2*i+3]);
		d_i[i] = (box_i_max[i]) / ((double) box_i_num[i]);
		systemSize *= box_i_num[i];
	}
	universe_grid = new VECT [systemSize];
	universe_grid_index = new int [systemSize];
	potential = new double [systemSize];
	//coor_temp = new double [dim_num];
	
	//Fill up the universe_grid
	/*for (i=0; i<systemSize; i++) {
		universe_grid[i].setDim(dim_num);
		sizes = 1;
		for (j=0; j<dim_num; j++) { //SHOULD try to replace this with indexDim function, but it still needs to be tested.
			if (j>0) {
				sizes *= box_i_num[j-1];
			}
			quo = (int) (i/sizes);
			rem = quo % box_i_num[j];
			universe_grid[i].coor[j] = rem;
		}
		//Generate a list of the 1D representation indices
		universe_grid_index[i] = i;
		
	}*/
	if ((systemSize*dim_num) > CHUNK_RATIO) {
		chunk = (int) (systemSize/CHUNK_RATIO);
	}
	else if ((systemSize*dim_num) > (CHUNK_RATIO/10)){
		chunk = (int) (systemSize/(CHUNK_RATIO/10));
	}
	else {
		chunk = 1;
	}
	/*
	// Max stack size of 
	if (chunk > 200000) {
		chunk = 200000;
	}

	chunk = 1000;*/
	
	cout << "Program Initialized." << endl;
	
#pragma omp parallel default(shared) private (i,j,thread)
	{
	thread = omp_get_thread_num();
	if (thread == 0) {
		nthreads = omp_get_num_threads();
		cout << "1st Parallel Section Number of Threads: " << nthreads << endl;
	}
	
#pragma omp for schedule(dynamic, chunk)  //Possibly try: http://software.intel.com/en-us/articles/openmp-loop-collapse-directive
	for (i=0; i<systemSize; i++) {
		universe_grid[i].setDim(dim_num);
		for (j=0; j<dim_num; j++) {
			universe_grid[i].coor[j] = indexDim(j, i, dim_num, box_i_num)*d_i[j]; //Seems to work, but not fully tested.
		}
		//Generate a list of the 1D representation indices
		universe_grid_index[i] = i;
		
	}
	}
	
	// Display the contents of universe_grid
	/*
	cout << "i  x  y  z  t" << endl;
	for (i=0; i<systemSize; i++) {
		cout << i << " ";
		for (j=0; j<dim_num; j++) {
			cout << universe_grid[i].coor[j] << "  ";
		}
		cout << endl;
	}
	 */
	
	
	/*Test code for disp function
	 index_test = new int [dim_num];
	 cout << "Input Index: ";
	 for (i=0; i<dim_num; i++) {
	 index_test[i] = atof(argv[(2*(dim_num-1) + 4) + i]);
	 cout << index_test[i] << " ";
	 }
	 cout << endl;
	 cout << "Output Displacement and Index:" << endl;
	 cout << "i  x  y  z  t" << endl;
	 
	 cout << disp(index_test, dim_num, box_i_num) << " ";
	 for (j=0; j<dim_num; j++) {
	 cout << universe_grid[disp(index_test, dim_num, box_i_num)].coor[j] << "  ";
	 }
	 cout << endl;
	 */
	
	cout << "Grid generated." << endl;
	
	// Determine which portion of the universe_grid each processor gets #########################
	
	// Set Charges
	sourceQ = 1; //e
	testQ = 1; //e
	
	coor_temp = new double [dim_num];
	
	for (i=0; i<dim_num; i++) {
		coor_temp[i] = 0.0;
	}
	
	sourceCharge.setDim(dim_num);
	sourceCharge.setCor(coor_temp, dim_num);

#pragma omp parallel default(shared) private (i,j)
	{
	
	thread = omp_get_thread_num();
	if (thread == 0) {
		nthreads = omp_get_num_threads();
		cout << "2nd Parallel Section Number of Threads: " << nthreads << endl;
	}
	
	#pragma omp for schedule(dynamic, chunk)
	//Calculate energy
	for (i=0; i<systemSize; i++) {
		potential[i] = CoulombEng (sourceCharge, sourceQ, universe_grid[i], testQ);
		
	}
	}
	
	potential[0] = 100000.0;
	/*
	for (i=0; i<systemSize; i += 100) {
			cout << "System Energy: " << potential[i] << " kJ/mol" << endl;
	}*/
	
	
	cout << "Potential Calculated." << endl;
	
	//Write the data to a cube file
	
	//Make generators for each dimension
	dimGenerator = new VECT [dim_num];
	for (i=0; i<dim_num; i++) {
		dimGenerator[i].setDim(dim_num);
		
		for (j=0; j<dim_num; j++) {
			dimGenerator[i].coor[j] = 0.0;
		}
		
		dimGenerator[i].coor[i] = d_i[i];
	}
	
	
	comments = new char* [2]; //The location of these memory allocations is important;I NEED to find out why BEFORE any code is published.
	comments[0] = new char [100];
	comments[1] = new char [100];
	filename = new char [100];
	
	nAtoms = 1;
	
	atomNum = new int [nAtoms];
	atomNum[0] = 1;
	
	atomCharges = new double [nAtoms];
	atomCharges[0] = 0.0;
	
	atomPos = new VECT [nAtoms];
	for (j=0; j<dim_num; j++) {
		atomPos[0].coor[j] = 0.0;
	}
	
	//cout << atomPos[0].coor[2] << endl;
	
	for (j=0; j<dim_num; j++) {
		origin.coor[j] = 0.0;
	}
	
	//cout << atomPos[0].coor[2] << endl;
	
	strcpy(comments[0], "Comment 1");
	strcpy(comments[1], "Comment 2");
	
	unitFlag = -1;
	
	//strcpy(filename, "testGaussian.cube");
	filename = argv[(dim_num*2)+2];	
	
	writeCubeFile(dim_num, box_i_num, dimGenerator, potential, systemSize, 
				  nAtoms, atomNum, atomCharges, atomPos, origin, comments, 
				  unitFlag, filename);
	
	cout << "Data written to file '" << filename << "'." << endl;

	return 0;
}
