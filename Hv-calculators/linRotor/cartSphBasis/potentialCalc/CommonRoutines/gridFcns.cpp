#include "gridFcns.h"

using namespace std;

//Changed d_x = x_max/n to d_x = x_max/(n-1) to get a proper grid - JTC 2013/08/20
universeProp* generateGrid(int numDim, double *gridMax, int *gridPoints) {
	universeProp *universe = new universeProp();
	int chunk, thread, nthreads, i, j;
	
	universe->numDim = numDim;
	universe->grid_max = gridMax;
	universe->grid_num = gridPoints;
	
	universe->d_i = new double [universe->numDim];
	
	universe->sysSize = 1;
	for (i=0; i<universe->numDim; i++) {
		universe->d_i[i] = (universe->grid_max[i]) / (double(universe->grid_num[i] - 1));
		//universe->d_i[i] = (universe->grid_max[i]) / ((double) universe->grid_num[i]);
		universe->sysSize *= universe->grid_num[i];
	}
	
	universe->grid = new VECT [universe->sysSize];
	
	
	//Determine for-loop chunk for OpenMP
	if (universe->sysSize > CHUNK_RATIO) {
		chunk = (int) (universe->sysSize/CHUNK_RATIO);
	}
	else if (universe->sysSize > (CHUNK_RATIO/10)){
		chunk = (int) (universe->sysSize/(CHUNK_RATIO/10));
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
	for (i=0; i<universe->sysSize; i++) {
		universe->grid[i].DIM(universe->numDim);
		for (j=0; j<universe->numDim; j++) {
			//universe->grid[i].COOR(j) = indexDim(j, i, universe->numDim, universe->grid_num)*universe->d_i[j]; //Seems to work, but not fully tested.
			universe->grid[i].COOR(j) = (double(indexDim(j, i, universe->numDim, universe->grid_num))*universe->d_i[j]) - (universe->grid_max[j]/2.0); //Seems to work, but not fully tested.; subtract 1/2*box_i_max to centre origin
		}		
	}
	}
	
	cout << "Grid generated." << endl;
	
	return universe;
}

/*This function returns the displacement for an n-dimensional system
 mapped to a 1D vector, given an n-dimensional index, the number of 
 dimensions and the size of each dimension.
 i.e. (1,2,1,4) -> 52
 */
int disp(int *index, int dim, int *dim_size)
{
	/* index	- an array containing the different indices in the n-dimensional basis
	   dim		- the number of dimensions
	   dim_size - an array containing the size of each dimension
	 */
	/* The displacement is displ = SUMoverDimensions(index[i]*PROD(from 0 to (i-1))(dim_size[i]))*/
	
	int displ = 0;
	int size_tot = 1;
	int i, j;
	
	for (i=0; i<dim; i++) {
		for (j=0; j<i; j++) { size_tot *= dim_size[j]; } //Product of the sizes of all dimensions less than i
														 // This is also the displacement for each dimension
		
		//Catch incorrect indices
		if ((index[i] >= dim_size[i])||(index[i] < 0)) {
			cerr << "Error, index out of bounds; disp failed." << endl;
			cerr << "Index " << i << ": " << index[i] << endl;
			exit (1);
		}
		
		displ += index[i]*size_tot; //Sum of indvidual displacements
		size_tot = 1;
	}
	return displ;
}

//This function still needs to be tested thoroughly
int indexDim(int dim, int disp, int dim_num, int *dim_size){
	int j, sizes, quo, rem;
	sizes = 1;
	
	if (dim_num <= dim) {
		cerr << "There are not this many dimensions present; indexDim failed." << endl;
		exit (1);
	}
	
	for (j=0; j<=dim; j++) {
		if (j>0) {
			sizes *= dim_size[j-1];
		}
		quo = (int) (disp/sizes);
		rem = quo % dim_size[j];
	}
	return rem;
}

/* Test code for disp()
int main(int argc, char** argv)
{
	int dim_num, *index, *dim_size;
	int i, j;
	
	//Get number of dimensions
	dim_num = atoi(argv[1]);
	
	index = new int [dim_num];
	dim_size = new int [dim_num];
	
	cout << "Sizes:";
	
	for (i=0; i<dim_num; i++) {
		dim_size[i] = atoi(argv[i+2]);
		index[i] = atoi(argv[i+2+dim_num]);
		cout << " " << dim_size[i];
	}
	
	cout << endl;
	cout << "Index:";
	for (i=0; i<dim_num; i++) {
		cout << " " << index[i];
	}
	
	
	cout << endl;
	cout << "Displacement: " << disp(index, dim_num, dim_size) << endl;
	
	
	return 0;
}
 */

//Function to print out a Gaussian Cube File
// unitFlag > 0 means Bohr, unitFlag < 0 means Angstroms NOTE: this flag doesn't seem to work properly right now.

void writeCubeFile(int nDim, int *dimSize, VECT *dimGenerator, double *data, int dataSize, int nAtoms, int *atomNum, double *atomCharge, VECT *atomPos, const VECT& origin, string *comments, int unitFlag, string filename) {
	int i, j, k, count;
	int displ, NumCol;
	int index[3];
	char unitSign[2];
	ofstream dataFile;
	
	dataFile.open(filename.c_str(), ios::out);
	
	if (!(dataFile.is_open())) {
		cerr << "The file '" << filename << "' could not be opened for writing purposes." << endl;
		exit(1);
	}
	
	if (nDim != 3) {
		cerr << "The Cube File can only handle 3D data." << endl;
		dataFile.close();
		exit (1);
	}
	
	//Determine units
	if (unitFlag > 0) { //Bohr
		strcpy(unitSign, "");
	}
	else if (unitFlag < 0) { //Angstroms
		strcpy(unitSign, "");
	}
	else {
		cerr << "Undefined units. Aborting Gaussian cube file output." << endl;
		exit (1);
	}
	
	//cout << unitSign << endl;
	/*
	for (j=0; j<nDim; j++) {
		cout << origin.COOR(j) << endl;
	}*/
	
	// Write comments
	dataFile << "# " << comments[0] << endl;
	dataFile << "# " << comments[1] << endl;
	
	//Write number of atoms and origin location
	dataFile << nAtoms << " " << origin.COOR(0) << " " << origin.COOR(1) << " " << origin.COOR(2) << endl;
	/*
	for (i=0; i<nDim; i++) {
		dataFile << " " << fixed << setprecision(5) << origin.COOR(i);
	}
	dataFile << endl;
	*/
	//Write dimension sizes and generators
	for (i=0; i<nDim; i++) {
		dataFile << unitSign << dimSize[i] << " " << dimGenerator[i].COOR(0) << " " << dimGenerator[i].COOR(1) << " " << dimGenerator[i].COOR(2) << endl;
	}
	
	//Write atom numbers, charges, and positions
	for (i=0; i<nAtoms; i++) {
		dataFile << atomNum[i] << " " << atomCharge[i] << " " << atomPos[i].COOR(0) << " " << atomPos[i].COOR(1) << " " << atomPos[i].COOR(2) << endl;
	}
	
	//Write potential values, looping over x first, then y, then z and having only NumCol values per line.
	NumCol = 6;
	count = 0;
	
	for (i=0; i<dimSize[0]; i++) {
		index[0] = i;
		
		for (j=0; j<dimSize[1]; j++) {
			index[1] = j;
			
			for (k=0; k<dimSize[2]; k++) {
				index[2] = k;
				
				displ = disp(index, nDim, dimSize); //Get position of data in the 1D data array given the index (i, j, k).
				
				dataFile << data[displ] << " ";
				
				//Add a newline character if file is too wide.
				if ((count % NumCol) == (NumCol - 1)) {
					dataFile << endl;
				}
				count++;
			}
		}
	}

	dataFile << endl;
	
	dataFile.close();
}
	



