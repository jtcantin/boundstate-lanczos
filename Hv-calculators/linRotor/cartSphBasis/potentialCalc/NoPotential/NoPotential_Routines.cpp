#include "Alavi_ContGrid_H2_Routines.h"

using namespace std;

//The following function calculates the energies for the Centre of Mass and Hydrogen sites of the Alavi hydrogen model
void ZeroPotential_Eng_contGrid(double *H2potential, double *H2potentialPI, universeProp *point_universe, interfaceStor *interface) {
	
	int i, thread, nthreads;
	
	//Quadrature Variables
	int na, nb;
	
	na = interface->quadrature->GLnum;
	
	nb = interface->quadrature->GCnum;
	
	int a, b, ind2;
	
	potentialGridSize = point_universe->sysSize * na * nb;
	
#pragma omp parallel default(shared) private (i, thread, H1pos, H2pos, ind, a, b, ind2, potential)
	{
	
	thread = omp_get_thread_num();
	if (thread == 0) {
		nthreads = omp_get_num_threads();
		cout << "Zero Potential Initialize Parallel Section Number of Threads: " << nthreads << endl;
	}
	
	//Calculate Energy for rotor in a given orientation
	
	//Calculate total H2 Energy
#pragma omp for schedule(guided) collapse(3)
	for (i=0; i<point_universe->sysSize; i++) {
		for (a=0; a<na; a++) {
			for (b=0; b<nb; b++) {
				
				ind2 = (i*na + a)*nb + b;
				
				H2potential[ind2] = 0.0;
				
				H2potentialPI[ind2] = 0.0;
		
			}
		}
	}	
	}

	cout << "Zero Potential Initialized." << endl;

}	

pointPotentialStorH2* preCalcPotential_ZeroPotential_ContGrid(int numDim, double *gridMax, int *gridPoints, string geometryFilename, interfaceStor *interface, int argc, char **argv) {
	
	//There are 3 spatial dimensions, so numDim = 3 and the gridMax and gridPoints arrays are of length 3
	
	//Get the grid for the potential using the same values as the system grid
	universeProp *potentialUniverse = new universeProp();
	fiveDGrid *grids;
	
	double *sysGridMax = new double [3];
	int *sysGridPoints = new int [3];
	
	grids = interface->grids;
	
	sysGridMax[0] = grids->x_max;
	sysGridPoints[0] = grids->nx;
	
	sysGridMax[1] = grids->y_max;
	sysGridPoints[1] = grids->ny;
	
	sysGridMax[2] = grids->z_max;
	sysGridPoints[2] = grids->nz;
	
	int ni = sysGridPoints[0];
	int nj = sysGridPoints[1];
	int nk = sysGridPoints[2];
	
	potentialUniverse->numDim = 3;
	potentialUniverse->grid_max = sysGridMax;
	potentialUniverse->grid_num = sysGridPoints;
	potentialUniverse->sysSize = ni*nj*nk;
	
	//Store dx, dy, dz
	potentialUniverse->d_i = new double [potentialUniverse->numDim];
	potentialUniverse->d_i[0] = sysGridMax[0] / double(sysGridPoints[0] - 1);
	potentialUniverse->d_i[1] = sysGridMax[1] / double(sysGridPoints[1] - 1);
	potentialUniverse->d_i[2] = sysGridMax[2] / double(sysGridPoints[2] - 1);
	
	//Build the grid
	try {
		potentialUniverse->grid = new VECT [potentialUniverse->sysSize];
    }
    catch( bad_alloc a) {
        const char * temp = a.what();
        cout << temp << endl;
        cout << "Threw a bad_alloc exception when allocating the potential vector grid." << endl;
    }
	
	int i, j, k;
	
	//Insert grid values
	int ind;
	
	//cout << "Beginning grid generation." << endl;
	
	#pragma omp parallel for default(shared) private (i, j, k, ind) schedule(guided) collapse(3)
	for (i=0; i<ni; i++) {
		for (j=0; j<nj; j++) {
			for (k=0; k<nk; k++) {
				
				ind = (i*nj + j)*nk + k;
				
				potentialUniverse->grid[ind].DIM(potentialUniverse->numDim);
				
				potentialUniverse->grid[ind].COOR(0) = grids->x_Grid[i];
				potentialUniverse->grid[ind].COOR(1) = grids->y_Grid[j];
				potentialUniverse->grid[ind].COOR(2) = grids->z_Grid[k];
			}
		}
	}
	
	cout << "Grid generated." << endl;
	
	//Get the final potentials
	int na, nb;
	
	na = interface->quadrature->GLnum;
	nb = interface->quadrature->GCnum;
	
	double *H2potential = new double [potentialUniverse->sysSize * na * nb];
	double *H2potentialPI = new double [potentialUniverse->sysSize * na * nb];
	
	ZeroPotential_Eng_contGrid(H2potential, H2potentialPI, potentialUniverse, interface);
	
	cout << "----------------------------------------------------------------------------" << endl;
	
	//Store everything for interface
	pointPotentialStorH2 *partialPotential = new pointPotentialStorH2();
	
	partialPotential->fullPotential2 = H2potential;
	partialPotential->fullPotential3 = H2potentialPI;
	
	partialPotential->potentialUniverse = potentialUniverse;
	
	return partialPotential;
}
