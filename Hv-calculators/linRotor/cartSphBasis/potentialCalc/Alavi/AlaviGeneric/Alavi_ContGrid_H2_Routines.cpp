#include "Alavi_ContGrid_H2_Routines.h"

using namespace std;

//The following function calculates the energies for the Centre of Mass and Hydrogen sites of the Alavi hydrogen model
void Alavi_SiteSite_Eng_contGrid(double *H2potential, double *H2potentialPI, universeProp *point_universe, sysAtoms *atomGeo, interfaceStor *interface) {
	
	int i, thread, nthreads;
	
	double (*SummedCoulombPotential)(VECT, double, char*, VECT*, int) = interface->fcnPointers->SummedCoulombPotential;
	double (*SummedLJPotentialFast)(VECT, char*, VECT*, int) = interface->fcnPointers->SummedLJPotentialFast;
	
	//Quadrature Variables
	double *cosThetaAbscissae, *cosPhiAbscissae, *phiAbscissae, *PIphiAbscissae, *thetaAbscissae;
	int na, nb;
	
	cosThetaAbscissae = interface->quadrature->GLabscissae;
	thetaAbscissae = interface->quadrature->GLacosAbscissae;
	na = interface->quadrature->GLnum;
	
	cosPhiAbscissae = interface->quadrature->GCabscissae;
	phiAbscissae = interface->quadrature->GCacosAbscissae;
	PIphiAbscissae = interface->quadrature->GCPIacosAbscissae;
	nb = interface->quadrature->GCnum;
	
	double *CMpotential = new double [point_universe->sysSize];
	
	VECT H1pos, H2pos;
	VECT *CM_H = new VECT [na*nb];
	VECT *CM_H_PI = new VECT [na*nb];
	
	int ind, a, b, ind2;
	double potential;
	double potCeil = interface->potentialCeiling;
	
#pragma omp parallel default(shared) private (i, thread, H1pos, H2pos, ind, a, b, ind2, potential)
	{
	
	thread = omp_get_thread_num();
	if (thread == 0) {
		nthreads = omp_get_num_threads();
		cout << "Alavi Point Potential Parallel Section Number of Threads: " << nthreads << endl;
	}
	
	//Calculate Energy for Centre of Mass
	
	//Calculate Coulomb energy
#pragma omp for schedule(guided)
	for (i=0; i<point_universe->sysSize; i++) {
		CMpotential[i] = (*SummedCoulombPotential)(point_universe->grid[i], Q_H2_CM, atomGeo->atomType, atomGeo->atomPos, atomGeo->nAtoms);
		
	}
	
	if (thread == 0) {
		cout << "Centre of Mass Coulomb Potential Calculated." << endl;
	}
	
	//Calculate Lennard Jones energy
#pragma omp for schedule(guided)
	for (i=0; i<point_universe->sysSize; i++) {
		CMpotential[i] += (*SummedLJPotentialFast)(point_universe->grid[i], atomGeo->atomType, atomGeo->atomPos, atomGeo->nAtoms);
		
	}
	
	if (thread == 0) {
		cout << "Centre of Mass Lennard-Jones Potential Calculated." << endl;
	}
	
	//Calculate Energy for H2 in a given orientation
	
	//Get unit vector pointing along H-H axis for both phi = [0,pi) and [pi, 2pi)
#pragma omp for schedule(guided) collapse(2)
	for (a=0; a<na; a++) {
		for (b=0; b<nb; b++) {
			
			ind = a*nb+b;
			
			CM_H[ind].DIM(3);
			CM_H[ind].COOR(0) = H2_H1_Z*sin(thetaAbscissae[a])*cosPhiAbscissae[b];
			CM_H[ind].COOR(1) = H2_H1_Z*sin(thetaAbscissae[a])*sin(phiAbscissae[b]);
			CM_H[ind].COOR(2) = H2_H1_Z*cosThetaAbscissae[a];
			
			CM_H_PI[ind].DIM(3);
			CM_H_PI[ind].COOR(0) = CM_H[ind].COOR(0);
			CM_H_PI[ind].COOR(1) = -1.0*(CM_H[ind].COOR(1));
			CM_H_PI[ind].COOR(2) = CM_H[ind].COOR(2);
		}
	}
			

	H1pos.DIM(3);
	H2pos.DIM(3);
	
	//Calculate total H2 Energy
#pragma omp for schedule(guided) collapse(3)
	for (i=0; i<point_universe->sysSize; i++) {
		for (a=0; a<na; a++) {
			for (b=0; b<nb; b++) {
				
				ind = a*nb+b;
				ind2 = (i*na + a)*nb + b;
				
				// Get positions for phi = [0, pi)
				H1pos = point_universe->grid[i] + CM_H[ind];
				H2pos = point_universe->grid[i] - CM_H[ind];
		
				potential = (*SummedCoulombPotential)(H1pos, Q_H2_H, atomGeo->atomType, atomGeo->atomPos, atomGeo->nAtoms);
				potential += (*SummedCoulombPotential)(H2pos, Q_H2_H, atomGeo->atomType, atomGeo->atomPos, atomGeo->nAtoms);
				potential += CMpotential[i];
				
				if (potential >= potCeil) {
					potential = potCeil;
				}
				
				H2potential[ind2] = potential;
				
				// Get positions for phi = [pi, 2pi)
				H1pos = point_universe->grid[i] + CM_H_PI[ind];
				H2pos = point_universe->grid[i] - CM_H_PI[ind];
				
				potential = (*SummedCoulombPotential)(H1pos, Q_H2_H, atomGeo->atomType, atomGeo->atomPos, atomGeo->nAtoms);
				potential += (*SummedCoulombPotential)(H2pos, Q_H2_H, atomGeo->atomType, atomGeo->atomPos, atomGeo->nAtoms);
				potential += CMpotential[i];
				
				if (potential >= potCeil) {
					potential = potCeil;
				}
				
				H2potentialPI[ind2] = potential;
		
			}
		}
	}	
	}

	cout << "H2 Potential Calculated." << endl;
	
	delete [] CMpotential;
	delete [] CM_H;
	delete [] CM_H_PI;
}	

//The following function precalculates the Alavi hydrogen site potentials by setting up the PES grid, getting the system geometry, 
//   and calculating the CM and H-atom potentials
pointPotentialStorH2* preCalcPotential_Alavi_ContGrid(int numDim, double *gridMax, int *gridPoints, string geometryFilename, interfaceStor *interface, int argc, char **argv) {
	
	//There are 3 spatial dimensions, so numDim = 3 and the gridMax and gridPoints arrays are of length 3
	
	//Get the grid for the potential using the same values as the system grid
	universeProp *potentialUniverse;
	
	double sysGridMax[3];
	int sysGridPoints[3];
	
	grids = interface->grids
	
	sysGridMax[0] = grids->x_max;
	sysGridPoints[0] = grids->nx;
	
	sysGridMax[1] = grids->y_max;
	sysGridPoints[1] = grids->ny;
	
	sysGridMax[2] = grids->z_max;
	sysGridPoints[2] = grids->nz;
	
	potentialUniverse->numDim = 3;
	potentialUniverse->grid_max = sysGridMax;
	potentialUniverse->grid_num = sysGridPoints;
	potentialUniverse->sysSize = sysGridMax[0]*sysGridMax[1]*sysGridMax[2];
	
	//Store dx, dy, dz
	potentialUniverse->d_i = new double [potentialUniverse->numDim];
	potentialUniverse->d_i[0] = sysGridMax[0] / double(sysGridPoints[0] - 1);
	potentialUniverse->d_i[1] = sysGridMax[1] / double(sysGridPoints[1] - 1);
	potentialUniverse->d_i[2] = sysGridMax[2] / double(sysGridPoints[2] - 1);
	
	//Build the grid
	potentialUniverse->grid = new VECT [potentialUniverse->sysSize];
	
	int i, j, k;
	
	//Initialize vectors
	for (i=0; i<potentialUniverse->sysSize; i++) {
		potentialUniverse->grid[i].DIM(potentialUniverse->numDim)
	}
	
	//Insert grid values
	int ind;
	
	int ni = sysGridPoints[0];
	int nj = sysGridPoints[1];
	int nk = sysGridPoints[2];
	
	for (i=0; i<ni; i++) {
		for (j=0; j<nj; j++) {
			for (k=0; k<nk; k++) {
				
				ind = (i*nj + j)*nk + k;
				
				potentialUniverse->grid[ind].COOR(0) = grids->x_Grid[i];
				potentialUniverse->grid[ind].COOR(1) = grids->y_Grid[j];
				potentialUniverse->grid[ind].COOR(2) = grids->z_Grid[k];
			}
		}
	}
	
	cout << "Grid generated." << endl;
	
	potentialUniverse = generateGrid(3, sysGridMax, sysGridPoints);
	
	//Check if the grids are the same
//	int i, j, k;
//	int index[3];
//	for (i=0; i<sysGridPoints[0] ; i++) {
//		for (j=0; j<sysGridPoints[1]; j++) {
//			for (k=0; k<sysGridPoints[2]; k++) {
//						
//				cout << "_______________________________________________________" << endl;
//				cout << "Point " << i << "," << j << "," << k << endl;
//				
//				index[0] = i;
//				index[1] = j;
//				index[2] = k;
//				
//				cout << potentialUniverse->grid[disp(index, 3, sysGridPoints)].COOR(0) - grids->x_Grid[i] << endl;
//				cout << potentialUniverse->grid[disp(index, 3, sysGridPoints)].COOR(1) - grids->y_Grid[j] << endl;
//				cout << potentialUniverse->grid[disp(index, 3, sysGridPoints)].COOR(2) - grids->z_Grid[k] << endl;
//				
//			}
//		}
//		
//	}
	
	//Get the system geometry for the water molecules
	sysAtoms *atomGeo = new sysAtoms();
	
	void (*getAtomGeometry)(char**, VECT**, int*, string) = interface->fcnPointers->getAtomGeometry;
	
	(*getAtomGeometry)(&(atomGeo->atomType), &(atomGeo->atomPos), &(atomGeo->nAtoms), geometryFilename);
	
	//Get the final potentials
	int na, nb;
	
	na = interface->quadrature->GLnum;
	nb = interface->quadrature->GCnum;
	
	double *H2potential = new double [potentialUniverse->sysSize * na * nb];
	double *H2potentialPI = new double [potentialUniverse->sysSize * na * nb];
	
	Alavi_SiteSite_Eng_contGrid(H2potential, H2potentialPI, potentialUniverse, atomGeo, interface);
	
	cout << "----------------------------------------------------------------------------" << endl;
	
	//Store everything for interface
	pointPotentialStorH2 *partialPotential = new pointPotentialStorH2();
	
	partialPotential->fullPotential2 = H2potential;
	partialPotential->fullPotential3 = H2potentialPI;
	
	partialPotential->potentialUniverse = potentialUniverse;
	
	//interface->atomGeo = atomGeo;
	
	delete [] atomGeo;
	
	return partialPotential;
}
