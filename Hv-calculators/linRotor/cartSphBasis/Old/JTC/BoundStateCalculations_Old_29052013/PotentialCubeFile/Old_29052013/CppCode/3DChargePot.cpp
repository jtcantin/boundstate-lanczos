
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include "gridFcns.h"
#include "Cvectors.h"
//#include "matvec.h"

using namespace std;

#define PI 3.14159265359 //From MMTKCNvect
#define EPS0 0.000572765638445 //From MMTK

double CoulombEng(const CNvect& q1, double q1_q, const CNvect& q2, double q2_q)
{
	double distance, k, Eng;
	CNvect diffVec;
	
	k = 1/(4*PI*EPS0);
	
	diffVec.setDim(q1.num_dim);
	
	cout << q1.num_dim << endl;
	cout << q2.num_dim << endl;
	
	diffVec = q2 - q1;
	
	distance = diffVec.eNorm();
	
	Eng = k*q1_q*q2_q/distance;
	
	return Eng;
}

/* To run:
 mpirun -np N 3DChargePot [dim] [max(1) num(1) max(2) num(2) ... max(dim) num(dim)]
 */
int main (int argc, char** argv)
{	
	
	int *box_i_num, dim_num;
	int *scatterSize, *scatterDisp;
	int localSystemSize, localSystemSizeRem, systemSize, quo, rem;
	double *box_i_max, *d_i, **i_grid, **local_i_grid;
	double *potential, *local_potential, *coor_temp;
	CNvect *universe_grid, *local_universe_grid;
	int *universe_grid_index, *local_universe_grid_index;
	
	int i, j, k, sizes, *index_test;
	
	CNvect sourceCharge, testCharge;
	double sourceQ, testQ;
	
	double Eng;
	
	//MPI Variables ###############################################################	
	int rank, size, root, last;
	MPI::Status status;
	
	if (argc <= (0+1)) {
		cout << "Invalid number of arguments, need at least one." << endl;
		exit (1);
	}
	
	//Initialize MPI ###############################################################	
	MPI::Init(argc, argv);
	rank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();
	root = 0; //Set process 0 to be the root process
	last = size - 1; //The last process
	
	// Get universe size ###############################################################
	
	//Get number of dimensions
	dim_num = atoi(argv[1]);
	
	if (argc != ((dim_num*2)+2)) {///////////////////////////////////////////////////////////////////
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
	universe_grid = new CNvect [systemSize];
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
	
	for (i=0; i<systemSize; i++) {
		universe_grid[i].setDim(dim_num);
		for (j=0; j<dim_num; j++) {
			universe_grid[i].coor[j] = indexDim(j, i, dim_num, box_i_num); //Seems to work, but not fully tested.
		}
		//Generate a list of the 1D representation indices
		universe_grid_index[i] = i;
		
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
	}*/
	 
	
	
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
	
	// Determine which portion of the universe_grid each processor gets #########################
	
	//Get local storage sizes; the remainder will be given to the last process
	localSystemSize = (int) (systemSize / size);
	localSystemSizeRem = systemSize % size;
	
	if (rank == 0) {
		cout << "Grid points per process: " << localSystemSize << endl;
	}
	
	// Alocate memory for each process
	
	scatterDisp = new int [size];
	scatterSize = new int [size];
	
	for (j=0; j<size; j++) { //j is for each process
		scatterSize[j] = localSystemSize;
		scatterDisp[j] = i*localSystemSize;
	}

	//Add remainder to last process
	scatterSize[last] += localSystemSizeRem;
	
	localSystemSize = scatterSize[rank];
	
	//Allocate local memory
	local_universe_grid = new CNvect [localSystemSize];
	local_potential = new double [localSystemSize];
	local_universe_grid_index = new int [localSystemSize];
	
	// Send each process its portion of the grid
	MPI::COMM_WORLD.Scatterv(universe_grid_index, scatterSize, scatterDisp, MPI::INT, 
							 local_universe_grid_index, localSystemSize, MPI::INT, root);

	//Fill up local universe grid
	for (i=0; i<localSystemSize; i++) {
		k = local_universe_grid_index[i];
		local_universe_grid[i].setDim(dim_num);
		for (j=0; j<dim_num; j++) {
			local_universe_grid[i].coor[j] = indexDim(j, k, dim_num, box_i_num); //Seems to work, but not fully tested.
		}
	}
	
	// Set Charges
	sourceQ = 1; //e
	testQ = 1; //e
	
	coor_temp = new double [dim_num];
	
	for (i=0; i<dim_num; i++) {
		coor_temp[i] = 0.0;
	}
	
	sourceCharge.setDim(dim_num);
	sourceCharge.setCor(coor_temp, dim_num);
	
	

	//Calculate energy
	for (i=0; i<localSystemSize; i++) {
		local_potential[i] = CoulombEng (sourceCharge, sourceQ, local_universe_grid[i], testQ);
		cout << "System Energy: " << Eng << " kJ/mol" << endl;
	}
	
	MPI::Finalize();	
	return 0;
}
