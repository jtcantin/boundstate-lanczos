#include "Alavi_H2_Routines.h"

using namespace std;

int AlaviMapH2GridToHGrid(H2_orient *H2_mol, universeProp *H_universe, double H_CM_dist) {
	int j, *indices, displacement;
	double c_i, c_i_min, d_c, e_i, numerator, theta, phi;
	
	phi = H2_mol->phi;
	theta = H2_mol->theta;
	
	indices = new int [H_universe->numDim];
	
	for (j=0; j<H_universe->numDim; j++) {
		
		c_i = H2_mol->CM->COOR(j);
		
		//Get unit vector pointing along H-H axis
		switch (j) {
			case 0:
				e_i = sin(theta)*cos(phi);
				break;
			case 1:
				e_i = sin(theta)*sin(phi);
				break;
			case 2:
				e_i = cos(theta);
				break;
		}
		
		//Calculate minimum grid value (for a [min max] interval centred about 0)
		c_i_min = (-1.0/2.0)*H_universe->grid_max[j];
		
		//Delta x
		d_c = H_universe->d_i[j];
		
		//Calculate nearest index
		numerator = c_i + e_i*H_CM_dist - c_i_min;			
		indices[j] = round(numerator/d_c); 
		
		//Keep index within confines of universe
		if (indices[j]<=0) {
			indices[j] = 0;
		}
		else if (indices[j]>=H_universe->grid_num[j]) {
			indices[j] = H_universe->grid_num[j]-1;
		}
	}
	//Get displacement in a 1D storage array corresponding to an (x, y, z, etc.) index in a numDim system
	displacement = disp(indices, H_universe->numDim, H_universe->grid_num);
	
	delete [] indices;
	
	return displacement;
}

void Alavi_H2_Eng(double *H2potential, double *CMpotential, double *H_potential, H2_orient *H2_mol, universeProp *point_universe) {
	
	int i, thread, nthreads, index, chunk;
	H2_orient H2_mol_local;
	
	H2_mol_local.phi = H2_mol->phi;
	H2_mol_local.theta = H2_mol->theta;
	
	if (point_universe->sysSize > CHUNK_RATIO) {
		chunk = (int) (point_universe->sysSize/CHUNK_RATIO);
	}
	else if (point_universe->sysSize > (CHUNK_RATIO/10)){
		chunk = (int) (point_universe->sysSize/(CHUNK_RATIO/10));
	}
	else {
		chunk = 1;
	}
	
#pragma omp parallel default(shared) private (i,thread,index,H2_mol_local)
	{
	thread = omp_get_thread_num();
	if (thread == 0) {
		nthreads = omp_get_num_threads();
		cout << "Alavi H2 Parallel Section Number of Threads: " << nthreads << endl;
	}
	//Calculate Energy for H2 molecule
	
#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<point_universe->sysSize; i++) {
		H2potential[i] = 0.0;
	}
	
	//Add Centre of Mass energy
#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<point_universe->sysSize; i++) {
		H2potential[i] += CMpotential[i];
	}
	
	if (thread == 0) {
		cout << "H2 Centre of Mass Potential Calculated." << endl;
	}
	
	//Add H1 energy
#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<point_universe->sysSize; i++) {
		
		H2_mol_local.CM = &(point_universe->grid[i]);
		
		//Get which point in the point potential corresponds to H1 of the H2 molecule at universe_grid[i] and theta and phi
		index = AlaviMapH2GridToHGrid(&H2_mol_local, point_universe, H2_H1_Z);
		
		//Add point potential to molecule potential
		H2potential[i] += H_potential[index];
	}
	
	if (thread == 0) {
		cout << "H2 H1 Potential Calculated." << endl;
	}
	
	//Add H2 energy
#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<point_universe->sysSize; i++) {
		
		H2_mol_local.CM = &(point_universe->grid[i]);
		
		//Get which point in the point potential corresponds to H2 of the H2 molecule at universe_grid[i] and theta and phi
		index = AlaviMapH2GridToHGrid(&H2_mol_local, point_universe, H2_H2_Z);
		
		//Add point potential to molecule potential
		H2potential[i] += H_potential[index];
	}
	
	if (thread == 0) {
		cout << "H2 H2 Potential Calculated." << endl;
	}
	}
}

void Alavi_point_Eng(double *CMpotential, double *Hpotential, universeProp *point_universe, sysAtoms *atomGeo) {
	
	int i, thread, nthreads, chunk;
	
	if (point_universe->sysSize > CHUNK_RATIO) {
		chunk = (int) (point_universe->sysSize/CHUNK_RATIO);
	}
	else if (point_universe->sysSize > (CHUNK_RATIO/10)){
		chunk = (int) (point_universe->sysSize/(CHUNK_RATIO/10));
	}
	else {
		chunk = 1;
	}
	
#pragma omp parallel default(shared) private (i,thread)
	{
	
	thread = omp_get_thread_num();
	if (thread == 0) {
		nthreads = omp_get_num_threads();
		cout << "Alavi Point Potential Parallel Section Number of Threads: " << nthreads << endl;
	}
	
	//Set potential to 0.0
#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<point_universe->sysSize; i++) {
		CMpotential[i] = 0.0;
		Hpotential[i] = 0.0;
	}
	
	
	if (thread == 0) {
		cout << "Potential Initialized." << endl;
	}
	
	//Calculate Energy for Centre of Mass
	
	//Calculate Coulomb energy
#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<point_universe->sysSize; i++) {
		CMpotential[i] += Q_TIP4P_Eng(point_universe->grid[i], Q_H2_CM, atomGeo->atomType, atomGeo->atomPos, atomGeo->nAtoms);
		
	}
	
	if (thread == 0) {
		cout << "Centre of Mass Coulomb Potential Calculated." << endl;
	}
	
	//Calculate Lennard Jones energy
#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<point_universe->sysSize; i++) {
		CMpotential[i] += LJ_TIP4P_Eng_Fast(point_universe->grid[i], atomGeo->atomType, atomGeo->atomPos, atomGeo->nAtoms);
		
	}
	
	if (thread == 0) {
		cout << "Centre of Mass Lennard-Jones Potential Calculated." << endl;
	}
	
	//Calculate Energy for H atom
	
	//Calculate Coulomb energy
#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<point_universe->sysSize; i++) {
		Hpotential[i] += Q_TIP4P_Eng(point_universe->grid[i], Q_H2_H, atomGeo->atomType, atomGeo->atomPos, atomGeo->nAtoms);
		
	}
	
	if (thread == 0) {
		cout << "H-atom Coulomb Potential Calculated." << endl;
	}
	
	}
}	


double Alavi_H2_Eng_Point(double *CMpotential, double *H_potential, H2_orient *H2_mol, universeProp *point_universe) {
	
	int i, index, *indices;
	double H2potential, grid_min;
	
	indices = new int [point_universe->numDim];
	
	
	
	//Calculate Energy for H2 molecule
	H2potential = 0.0;
	
	//
	//Add Centre of Mass energy
	//
	for (i=0; i<point_universe->numDim; i++) {
		grid_min = (-1.0/2.0)*point_universe->grid_max[i]; //Calculate minimum grid value (for a [min max] interval centred about 0)
		indices[i] = round((H2_mol->CM->COOR(i) - grid_min)/point_universe->d_i[i]); //Get index corresponding to spatial value
	}
	
	//Determine displacement within CMpotential
	index = disp(indices, point_universe->numDim, point_universe->grid_num);
	
	H2potential += CMpotential[index];
	
	//
	//Add H1 energy
	//
	
	//Get which point in the point potential corresponds to H1 of the H2 molecule at CM, theta, and phi
	index = AlaviMapH2GridToHGrid(H2_mol, point_universe, H2_H1_Z);
	
	//Add point potential to molecule potential
	H2potential += H_potential[index];
	
	//
	//Add H2 energy
	//
	
	//Get which point in the point potential corresponds to H2 of the H2 molecule at CM, theta, and phi
	index = AlaviMapH2GridToHGrid(H2_mol, point_universe, H2_H2_Z);
	
	//Add point potential to molecule potential
	H2potential += H_potential[index];
	
	delete [] indices;
	
	return H2potential;
}



