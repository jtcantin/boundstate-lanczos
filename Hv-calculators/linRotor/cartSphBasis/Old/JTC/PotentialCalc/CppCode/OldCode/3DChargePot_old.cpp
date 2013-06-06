
#include <iostream>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include "gridFcns.h"
#include "Cvectors.h"

using namespace std;

#define PI 3.14159265359 //From MMTK
#define EPS0 0.000572765638445 //From MMTK

double CoulombEng(const C3vect& q1, double q1_q, const C3vect& q2, double q2_q)
{
	double distance, k, Eng;
	C3vect diffVec;
	
	k = 1/(4*PI*EPS0);
	
	diffVec = q2 - q1;
	
	distance = diffVec.eNorm();
	
	Eng = k*q1_q*q2_q/distance;
	
	return Eng;
}


int main (int argc, char** argv)
{	
	int *box_i_num, *local_box_i_num, *local_box_i_num_rem, dim_num;
	int **scatterSize, **scatterDisp;
	double *box_i_max, *d_i, **i_grid, **local_i_grid;
	double **potential, ***local_potential;
	
	int i, j;
	
	C3vect sourceCharge, testCharge;
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
	
	i_grid = new double* [dim_num];
	potential = new double* [dim_num];
	
	if (argc != ((dim_num*2)+2)) {
		cout << "Invalid number of arguments, need " << (dim_num*2)+1; 
		cout << "arguments for " << dim_num << " dimensions." << endl;
		exit (1);
	}
	
	//Get dimension sizes
	for (i=0; i<dim_num; i++) {
		box_i_max[i] = atof(argv[2*i+2]);
		box_i_num[i] = atoi(argv[2*i+3]);
		d_i[i] = (2*box_i_max[i]) / ((double) box_i_num[i]);
	}

	
	// Determine which portion of the grid each processor gets #########################	
	if (rank == 0)
	{
		//Form grid from -box_i_max to box_i_max for each basis (could have cleaner and more flexible code with a double pointer).
		for (i=0; i<dim_num; i++) {
			
			//Allocate memory
			i_grid[i] = new double [box_i_num[i]];
			potential[i] = new double [box_i_num[i]];
			
			
			//Fill the grid
			for (j=0; j<box_i_num[i]; j++) {
				i_grid[i][j] = (j*d_i[i]) - box_i_max[i];
			}
			
			//Get local storage sizes; the remainder will be given to the last process
			local_box_i_num[i] = (int) box_i_num[i] / size;
			local_box_i_num_rem[i] = box_i_num[i] % size;
			
			cout << "Dimension " << i << " grid size per process: " << local_box_i_num[i] << endl;
			
		}		
	}
	
	//Broadcast local gridsizes and remainder
	for (i=0; i<dim_num; i++) {
		MPI::COMM_WORLD.Bcast((local_box_i_num + i), 1, MPI::INT, root);
		MPI::COMM_WORLD.Bcast((local_box_i_num_rem + i), 1, MPI::INT, root);
	}
	
	// Alocate memory for each dimension
	scatterSize = new int* [dim_num];
	scatterDisp = new int* [dim_num];
	
	
	for (i=0; i<dim_num; i++) { //i is for each dimension
		//Alocate memory for each process
		scatterSize[i] = new int [size];
		scatterDisp[i] = new int [size];
		
		for (j=0; j<size; j++) { //j is for each process
			scatterSize[i][j] = local_box_i_num[i];
			scatterDisp[i][j] = i*local_box_i_num[i];
		}
	}
	
	//Add remainder to last processor and ensure all processes have the right grid size
	for (i=0; i<dim_num; i++) {
		scatterSize[i][last] += local_box_i_num_rem[i];
		local_box_i_num[i] = scatterSize[i][rank];
		
		//Allocate local grid memory
		local_i_grid[i] = new double [local_box_i_num[i]];
		
	}
	
	// Send each process its portion of the grid
	for (i=0; i<dim_num; i++) {
		MPI::COMM_WORLD.Scatterv(i_grid[i], scatterSize[i], scatterDisp[i], MPI::DOUBLE, local_i_grid[i], local_box_i_num[i], MPI::DOUBLE, root);
	}
					 
							 
	// Set Charges
	sourceQ = 1; //e
	testQ = 1; //e
	
	sourceCharge.setCor(0.0, 0.0, 0.0);
	
	testCharge.setCor(1.0, 1.0, 1.0);
	
	Eng = CoulombEng (sourceCharge, sourceQ, testCharge, testQ);
	
	cout << "System Energy: " << Eng << " kJ/mol" << endl;
	
	
	
	MPI::Finalize();	
	return 0;
}
