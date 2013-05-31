#include "TIP4P_AH2_EngRoutines.h"
#include <omp.h>
#include "gridFcns.h"

#define CHUNK_RATIO 100 //How large the chunks should be relative to the total data size; 10 or 100 seem to be reasonable values.

//FOR TESTING PURPOSES
#define	PERIODIC
#define	SYS_SIZE 1.731 //nm; sII unitcell parameter from Takeuchi_2013
#define	C_CUT_OFF 50 //nm; cutoff used for Coulomb potential by Takeuchi_2003

using namespace std;

double Q_TIP4P_EngTest(VECT pos, double q, char *atomType, VECT *atomPos, int nAtoms) {
	int i;
	double Eng;
	
	VECT diffVec;
	double distance;
	
	Eng = 0;
	
	for (i=0; i<nAtoms; i++) {
		switch (atomType[i]) {
			case 'H':
				Eng += CoulombEng(pos, q, atomPos[i], Q_H);
				break;
			case 'M':
				Eng += CoulombEng(pos, q, atomPos[i], Q_M);
				break;
			case 'O':
				diffVec.DIM(pos.NUM_DIM);
				diffVec = atomPos[i] - pos;
				distance = diffVec.eNorm();
				if (distance>C_CUT_OFF) {
					i += 3;
				}
				break;
			default:
				cerr << "Atom type not recognized while executing Q_TIP4P_Eng()." << endl;
				exit(1);
				break;
		}
	}
	
	return Eng;
}

// To run: ./TIP4P_AH2_EngRoutinesTest_Periodic [water_molecule.xyz file] [testPos.xyz file] [sys_size] [numReplicates]
int main(int argc, char **argv) {
	string filename, posFile, *cageType;
	string line;
	char *atomType;
	VECT *atomPos, pos, *cagePos;
	int nAtoms;
	double Eng, C_Eng, LJ_Eng, *C_engArray, *cageEng, sys_size;
	int i, j, k, l, m, replicates, replicateSize, repSize[3], nCages;
	
	int chunk, nthreads, thread;
	
	//Get command line arguments
	if (argc != (4+1)) {
		cerr <<	"Error, incorrect number of arguments" << endl;
		exit(1);
	}
	filename = argv[1];
	posFile = argv[2];
	sys_size = atof(argv[3]); //nm
	replicates = atoi(argv[4]);
	
	replicateSize = replicates*replicates*replicates; //Total number of replicates (includes original cell)
	C_engArray = new double [replicateSize];
	
	
	cout << "Program Initialized." << endl;
	
	// Read in water molecule locations
	getTIP4Patoms(&atomType, &atomPos, &nAtoms, filename);
	
	cout << "Atom Positions read from file '" << filename << "'." << endl;
	
	pos.DIM(3);
	
	
	// Read in test positions
	ifstream datafile;
	
	datafile.open(posFile.c_str(), ios::in);
	if (datafile.is_open()) {
		datafile>>nCages;
		getline (datafile, line); //Eat up end of line character from first line
		
		getline (datafile, line); //Get rid of comment line
		
		cagePos = new VECT[nCages]; //Values in Angstroms
		cageType = new string [nCages];
		
		for (i=0;i<nCages;i++) {
			
			cagePos[i].DIM(3);
			
			//Read in position label
			datafile>>cageType[i];
			
			//Read in positions
			datafile>>cagePos[i].COOR(0); 		 		
			datafile>>cagePos[i].COOR(1); 		
			datafile>>cagePos[i].COOR(2); 		
			
			//Change to nm from Ang
			cagePos[i].COOR(0) *= NM_PER_ANG; 		 		
			cagePos[i].COOR(1) *= NM_PER_ANG; 		
			cagePos[i].COOR(2) *= NM_PER_ANG;
		}
		
		cout << "Test locations read from '" << posFile << "'." << endl;
	}
	else {
		cerr << "File " << filename << " could not be opened." << endl;
		exit(1);
	}
	
	
	datafile.close();
	
	cageEng = new double [nCages];
	
	cout << "Beginning potential calculation." << endl;
	
	// Calculate potential for each test location (number of cages)
	for (m=0; m<nCages; m++) {
		
		Eng = 0.0; //kJ/mol
		C_Eng = 0.0;
		LJ_Eng = 0.0;
		
		cout << "Cage " << cageType[m] << ": ";
		
		
		//Set and shift test location to be in the middle of the replicate unit cells
		for (l=0; l<3; l++) {
			pos.COOR(l) = cagePos[m].COOR(l);
			pos.COOR(l) += ((int) (replicates/2))*sys_size;
		}
		
		
		chunk = 1;
		
		repSize[0] = replicates;
		repSize[1] = replicates;
		repSize[2] = replicates;
		
		//Parallelize the calculation across each replicate unit cell
		#pragma omp parallel default(shared) private (i,j,k,l,m,thread)
		{
		thread = omp_get_thread_num();
		if (thread == 0) {
			nthreads = omp_get_num_threads();
			//cout << "Parallel Section Number of Threads: " << nthreads << endl;
		}
		
		// Make a replicates by replicates by replicates supercell of the original unit cell and use this to calculate the energy
		// Each cell is not actually made and stored, but the positions of each atom are shifted and the energy calculated for each new location
		// The test charge is kept in a unit cell that is at the centre of the super cell
		// There is a cutoff as specified by C_CUT_OFF above.
		
		#pragma omp for schedule(dynamic, chunk) collapse (3) // Can't compile "collapse (3)" on the Mac as you need g++ v4.4 or later http://stackoverflow.com/questions/5882095/collapse-clause-being-ignored-in-pragma-omp-for
		for (i=0; i<replicates; i++) {
			for (j=0; j<replicates; j++) {
				for (k=0; k<replicates; k++) {
					
					//Declare these in the loop to keep them private (should try to move these out of the loop, it may make it faster; pointers are an issue, though)
					int index[3], displacement;
					VECT *atomPos2, displ;
					
					atomPos2 = new VECT [nAtoms];
					
					displ.DIM(3);
					
					// Calculate displacement vectors
					displ.COOR(0) = i*sys_size;
					displ.COOR(1) = j*sys_size;
					displ.COOR(2) = k*sys_size;
					
					//Shift each atom appropriately
					for (l=0; l<nAtoms; l++) {
						atomPos2[l].DIM(3);
						atomPos2[l] = atomPos[l] + displ;
					}
					
					//Set up energy index for disp() function
					index[0] = i;
					index[1] = j;
					index[2] = k;
					
					// Determine what displacement there is in a 1D vector storing this 3D system
					displacement = disp(index, 3, repSize);
					
					// Calculate the energy due to the replicate unit cell
					C_engArray[displacement] = Q_TIP4P_EngTest(pos, 1.0, atomType, atomPos2, nAtoms);
					
					delete [] atomPos2; // Prevent a memory leak	
				} 
			}
		}
		}
		
		//Sum up energy from each replicate
		for (i=0; i<replicateSize; i++) {
			C_Eng += C_engArray[i];
		}
		cout << setprecision(6) << C_Eng << " kJ/mol" << endl;
		
	}
	
	delete [] atomPos;
	delete [] atomType;
}

//	//Output data to test in VMD
//	ofstream outdatafile;
//	outdatafile.open("unitcell.xyz", ios::out);
//	if (outdatafile.is_open()) {
//		outdatafile << (nAtoms+1) << endl;
//		outdatafile << "# Comment" << endl;
//		for (i=0; i<(nAtoms); i++) {
//			outdatafile << atomType[i] << " ";
//			for (j=0; j<3; j++) {
//				outdatafile << atomPos[i].COOR(j)/NM_PER_ANG << " ";
//			}
//			outdatafile << endl;
//			
//		}
//		outdatafile << "Ar" << " "; 
//		for (j=0; j<3; j++) {
//			outdatafile << pos.COOR(j)/NM_PER_ANG << " ";
//		}
//		outdatafile << endl;
//	}
//		
//	outdatafile.close();
