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

using namespace std;

#define CHUNK_RATIO 100 //How large the chunks should be relative to the total data size; 10 or 100 seem to be reasonable values.

#define NM_PER_BOHR 0.052918
#define ANG_PER_BOHR 0.52918

/* To run:
 ./PointPotentials_TIP4P_AH2 [dim] [max(1) num(1) max(2) num(2) ... max(dim) num(dim) inputFilename OutputFilename]
 */
int main (int argc, char** argv)
{	
	
	int *box_i_num, dim_num;
	int systemSize;
	double *box_i_max, *d_i;
	double *potential, *coor_temp;
	VECT *universe_grid;
	int *universe_grid_index;
	
	double *CMpotential, *Hpotential, *H2potential;
	
	VECT *dimGenerator, origin;
	int *atomNum, unitFlag;
	double *atomCharges;
	string comments[2];
	string filename, filenameCM, filenameH, filenameH2;
	
	//TIP4P potential variables
	char *atomType;
	VECT *atomPos, pos;
	int nAtoms;
	double *Eng;
	
	int i, j;
	
	VECT sourceCharge, testCharge;
	double sourceQ, testQ;
	
	double phi, theta;
	
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
	
	if (argc != ((dim_num*2)+2+1+1+2)) { //dim_max dim_size, dim_num+runCommand, filename_in, filename_out, angles
		cout << "Invalid number of arguments, need " << (dim_num*2)+1+1+1; 
		cout << " arguments for " << dim_num << " dimensions." << endl;
		exit (1);
	}
	
	//Get dimension sizes
	box_i_max = new double [dim_num];
	box_i_num = new int [dim_num];
	d_i = new double [dim_num]; //ODDLY, this seems to have given a value of 0 once, though that was before the memory problems with origin and atompos were fixed.
	systemSize = 1;
	for (i=0; i<dim_num; i++) {
		box_i_max[i] = atof(argv[2*i+2]);
		box_i_num[i] = atoi(argv[2*i+3]);
		d_i[i] = (box_i_max[i]) / ((double) box_i_num[i]);
		systemSize *= box_i_num[i];
	}
	universe_grid = new VECT [systemSize];
	universe_grid_index = new int [systemSize];
	CMpotential = new double [systemSize];
	Hpotential = new double [systemSize];
	H2potential = new double [systemSize];
	//coor_temp = new double [dim_num];
	
	phi = atof(argv[2*dim_num+2]);
	theta = atof(argv[2*dim_num+2+1]);
	
	//cout << "phi: " << phi << ", theta: " << theta << endl;
	
	//Fill up the universe_grid
	/*for (i=0; i<systemSize; i++) {
		universe_grid[i].DIM(dim_num);
		sizes = 1;
		for (j=0; j<dim_num; j++) { //SHOULD try to replace this with indexDim function, but it still needs to be tested.
			if (j>0) {
				sizes *= box_i_num[j-1];
			}
			quo = (int) (i/sizes);
			rem = quo % box_i_num[j];
			universe_grid[i].COOR(j) = rem;
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
		universe_grid[i].DIM(dim_num);
		for (j=0; j<dim_num; j++) {
			//universe_grid[i].COOR(j) = indexDim(j, i, dim_num, box_i_num)*d_i[j]; //Seems to work, but not fully tested.
			universe_grid[i].COOR(j) = (indexDim(j, i, dim_num, box_i_num)*d_i[j]) - (box_i_max[j]/2.0); //Seems to work, but not fully tested.; subtract 1/2*box_i_max to centre origin
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
			cout << universe_grid[i].COOR(j) << "  ";
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
	 cout << universe_grid[disp(index_test, dim_num, box_i_num)].COOR(j) << "  ";
	 }
	 cout << endl;
	 */
	
	cout << "Grid generated." << endl;
	
	
	//Acquire system geometry
	filename = argv[(dim_num*2)+4];
	
	getTIP4Patoms(&atomType, &atomPos, &nAtoms, filename);
	
	
	//Set up Coulomb Potential Calculation
	
	// Set Charges
	sourceQ = 1; //e
	testQ = 1; //e
	
	coor_temp = new double [dim_num];
	
	for (i=0; i<dim_num; i++) {
		coor_temp[i] = 0.0;
	}
	
	sourceCharge.DIM(dim_num);
	sourceCharge.SET_COOR(coor_temp, dim_num);
	
	

#pragma omp parallel default(shared) private (i,j,thread)
	{
	
	thread = omp_get_thread_num();
	if (thread == 0) {
		nthreads = omp_get_num_threads();
		cout << "2nd Parallel Section Number of Threads: " << nthreads << endl;
	}
	
	//Set potential to 0.0
	#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<systemSize; i++) {
		CMpotential[i] = 0.0;
		Hpotential[i] = 0.0;
		H2potential[i] = 0.0;
		
	}
	
	
	if (thread == 0) {
		cout << "Potential Initialized." << endl;
	}
	
	//Calculate Energy for Centre of Mass
	
	//Calculate Coulomb energy
	#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<systemSize; i++) {
		CMpotential[i] += Q_TIP4P_Eng(universe_grid[i], Q_H2_CM, atomType, atomPos, nAtoms);
		
	}
	
	if (thread == 0) {
		cout << "Centre of Mass Coulomb Potential Calculated." << endl;
	}
	
	//Calculate Lennard Jones energy
	#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<systemSize; i++) {
		CMpotential[i] += LJ_TIP4P_Eng_Fast(universe_grid[i], atomType, atomPos, nAtoms);
		
	}
	
	if (thread == 0) {
		cout << "Centre of Mass Lennard-Jones Potential Calculated." << endl;
	}
	
	//Calculate Energy for H atom
	
	//Calculate Coulomb energy
	#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<systemSize; i++) {
		Hpotential[i] += Q_TIP4P_Eng(universe_grid[i], Q_H2_H, atomType, atomPos, nAtoms);
		
	}
	
	if (thread == 0) {
		cout << "H-atom Coulomb Potential Calculated." << endl;
	}
	
	//Calculate Energy for H2 molecule
	
	
	
	//Add Centre of Mass energy
#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<systemSize; i++) {
		H2potential[i] += CMpotential[i];
	}
	
	if (thread == 0) {
		cout << "H2 Centre of Mass Potential Calculated." << endl;
	}
	
	int index, *indices;
	
	indices = new int [dim_num];
	double x_i, x_i_min, d_x, e_i, numerator;
	
	//Add H1 energy
#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<systemSize; i++) {
		for (j=0; j<dim_num; j++) {
			x_i = universe_grid[i].COOR(j);
			switch (j) {
				case 0:
					e_i = sin(theta)*cos(phi);
					//cout << j << ", e_i: " << e_i << endl;
					break;
				case 1:
					e_i = sin(theta)*sin(phi);
					//cout << j << ", e_i: " << e_i << endl;
					break;
				case 2:
					e_i = cos(theta);
					//cout << j << ", e_i: " << e_i << endl;
					break;
			}
			x_i_min = (-1.0/2.0)*box_i_max[j];
			
			//cout << box_i_max[j] << endl;
			d_x = d_i[j];
			numerator = x_i + e_i*H2_H1_Z - x_i_min;			
			indices[j] = round(numerator/d_x); //add (box_i_max[j]/2.0) to displace coodinate
			
			//Keep index within confines of universe
			if (indices[j]<=0) {
				indices[j] = 0;
			}
			else if (indices[j]>=box_i_num[j]) {
				indices[j] = box_i_num[j]-1;
			}
			
			//if (thread == 0) {
//				cout << "Index " << j << ":" << endl;
//				cout << "x_i" << ", " << "e_i" << ", " << "x_i_min" << ", " << "d_x" << ", " << "indices[j]" << ", " << "numerator" << endl;
//				cout << x_i << ", " << e_i << ", " << x_i_min << ", " << d_x << ", " << indices[j] << ", " << numerator << endl;
//			}
		}
		index = disp(indices, dim_num, box_i_num);
		//if (thread == 0) {
//			cout << index << endl;
//			cout << universe_grid[index].COOR(0) << ", " << universe_grid[index].COOR(1) << ", " << universe_grid[index].COOR(2) << endl;
//		}
		
		H2potential[i] += Hpotential[index];
	}
	
	if (thread == 0) {
		cout << "H2 H1 Potential Calculated." << endl;
	}
	
	//Add H2 energy
#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<systemSize; i++) {
		for (j=0; j<dim_num; j++) {
			x_i = universe_grid[i].COOR(j);
			switch (j) {
				case 0:
					e_i = sin(theta)*cos(phi);
					//cout << j << ", e_i: " << e_i << endl;
					break;
				case 1:
					e_i = sin(theta)*sin(phi);
					//cout << j << ", e_i: " << e_i << endl;
					break;
				case 2:
					e_i = cos(theta);
					//cout << j << ", e_i: " << e_i << endl;
					break;
			}
			x_i_min = (-1.0/2.0)*box_i_max[j];
			
			//cout << box_i_max[j] << endl;
			d_x = d_i[j];
			numerator = x_i + e_i*H2_H2_Z - x_i_min;			
			indices[j] = round(numerator/d_x); //add (box_i_max[j]/2.0) to displace coodinate
			
			//Keep index within confines of universe
			if (indices[j]<=0) {
				indices[j] = 0;
			}
			else if (indices[j]>=box_i_num[j]) {
				indices[j] = box_i_num[j]-1;
			}
			
			//if (thread == 0) {
			//				cout << "Index " << j << ":" << endl;
			//				cout << "x_i" << ", " << "e_i" << ", " << "x_i_min" << ", " << "d_x" << ", " << "indices[j]" << ", " << "numerator" << endl;
			//				cout << x_i << ", " << e_i << ", " << x_i_min << ", " << d_x << ", " << indices[j] << ", " << numerator << endl;
			//			}
		}
		index = disp(indices, dim_num, box_i_num);
		//if (thread == 0) {
		//			cout << index << endl;
		//			cout << universe_grid[index].COOR(0) << ", " << universe_grid[index].COOR(1) << ", " << universe_grid[index].COOR(2) << endl;
		//		}
		
		H2potential[i] += Hpotential[index];
	}
	
	if (thread == 0) {
		cout << "H2 H2 Potential Calculated." << endl;
	}
	
	}
	
	//potential[0] = 100000.0;
	/*
	for (i=0; i<systemSize; i += 100) {
			cout << "System Energy: " << potential[i] << " kJ/mol" << endl;
	}
	
	
	
	cout << "Potential Calculated." << endl;
	
	/*
	//Write the data to a cube file
	*/
	//Make generators for each dimension
	dimGenerator = new VECT [dim_num];
	for (i=0; i<dim_num; i++) {
		dimGenerator[i].DIM(dim_num);
		
		for (j=0; j<dim_num; j++) {
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
	origin.DIM(dim_num);
	for (j=0; j<dim_num; j++) {
		origin.COOR(j) = -1.0*box_i_max[j]/2/NM_PER_BOHR; //Convert to Bohr for Cube file
		//cout << origin.COOR(j) << endl;
	}
	
	for (i=0; i<nAtoms; i++) {
		for (j=0; j<dim_num; j++) {
			atomPos[i].COOR(j) /= NM_PER_BOHR; //Convert to Bohr for Cube file
		}
	}
	
		
	unitFlag = -1;
	
	comments[0] = "Comment 1";
	comments[1] = "Comment 2";
	
	//strcpy(filename, "testGaussian.cube");
	filename = argv[(dim_num*2)+5];	
	
	filenameCM = filename + "_CM.cube";
	
	/*
	for (j=0; j<dim_num; j++) {
		cout << origin.COOR(j) << endl;
	}*/
	
	writeCubeFile(dim_num, box_i_num, dimGenerator, CMpotential, systemSize, 
				  nAtoms, atomNum, atomCharges, atomPos, origin, comments, 
				  unitFlag, filenameCM);
	
	cout << "Centre of Mass data written to file '" << filenameCM << "'." << endl;
	
	filenameH = filename + "_H.cube";
	
	writeCubeFile(dim_num, box_i_num, dimGenerator, Hpotential, systemSize, 
				  nAtoms, atomNum, atomCharges, atomPos, origin, comments, 
				  unitFlag, filenameH);
	
	cout << "H-atom data written to file '" << filenameH << "'." << endl;
	
	
	filenameH2 = filename + "_H2.cube";
	
	writeCubeFile(dim_num, box_i_num, dimGenerator, H2potential, systemSize, 
				  nAtoms, atomNum, atomCharges, atomPos, origin, comments, 
				  unitFlag, filenameH2);
	
	cout << "H2 data written to file '" << filenameH2 << "'." << endl;
	
	
	
	//Free memory
	delete [] atomNum;
	delete [] atomCharges;
	delete [] atomPos;	
	delete [] dimGenerator;
	delete [] coor_temp;
	delete [] box_i_max;
	delete [] box_i_num;
	delete [] d_i;
	delete [] universe_grid;
	delete [] universe_grid_index;
	delete [] CMpotential;
	delete [] Hpotential;
	delete [] H2potential;

	return 0;
}
