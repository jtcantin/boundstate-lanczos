#include "Alavi_H2_Routines.h"

using namespace std;

void Alavi_TIP4P_point_Eng(double *CMpotential, double *Hpotential, universeProp *point_universe, sysAtoms *atomGeo) {
	
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

pointPotentialStorH2* preCalcPotential_Alavi_TIP4P(int numDim, double *gridMax, int *gridPoints, string geometryFilename){
	
	//There are 3 spatial dimensions, so numDim = 3 and the gridMax and gridPoints arrays are of length 3
	
	//Get the grid for the potential
	universeProp *potentialUniverse;
	
	potentialUniverse = generateGrid(numDim, gridMax, gridPoints);
	
	//Get the system geometry for TIP4P molecules
	sysAtoms *atomGeo = new sysAtoms();
	
	getTIP4Patoms(&(atomGeo->atomType), &(atomGeo->atomPos), &(atomGeo->nAtoms), geometryFilename);
	
	//Get the partial potentials
	double *CMpotential, *Hpotential;
	
	CMpotential = new double [potentialUniverse->sysSize];
	Hpotential = new double [potentialUniverse->sysSize];
	
	Alavi_TIP4P_point_Eng(CMpotential, Hpotential, potentialUniverse, atomGeo);
	
	cout << "----------------------------------------------------------------------------" << endl;
	
	//Store everything for interface
	pointPotentialStorH2 *partialPotential = new pointPotentialStorH2();
	
	partialPotential->CMpotential = CMpotential;
	partialPotential->H_potential = Hpotential;
	partialPotential->potentialUniverse = potentialUniverse;
	
	//delete [] atomGeo->atomType;
	//	delete [] atomGeo->atomPos;
	delete atomGeo;
	
	return partialPotential;
}

