#include "TIP4P_AH2_EngRoutines.h"
#include <omp.h>
#include "gridFcns.h"

#define CHUNK_RATIO 100 //How large the chunks should be relative to the total data size; 10 or 100 seem to be reasonable values.

using namespace std;

void getTIP4Patoms(char **atomType, VECT **atomPos, int *nAtoms, string filename) {
	ifstream datafile;
	string line;
	VECT *EA, *CM;
	int nMol;
	
	int i, j;
	
	
	//Read in the Centres of Mass and Euler Angles
	datafile.open(filename.c_str(), ios::in);
	if (datafile.is_open()) {
		datafile>>nMol;
		//cout<<nMol<<endl;		
		getline (datafile, line);//Eat up end of line character from first line
		
		// get rid of comment line
		getline (datafile, line);
		//cout<<line<<endl;
		
		EA = new VECT[nMol]; //Values in radians
		CM = new VECT[nMol]; //Values in Angstroms
		
		char name[12];
		for (i=0;i<nMol;i++) {
			EA[i].DIM(3);
			CM[i].DIM(3);
			
			datafile>>name;
			//cout<<name<<"  ";
			datafile>>CM[i].COOR(0); 		
			datafile>>EA[i].COOR(0); 		
			datafile>>CM[i].COOR(1); 		
			datafile>>EA[i].COOR(1); 		
			datafile>>CM[i].COOR(2); 		
			datafile>>EA[i].COOR(2); 
			/*
			 cout << scientific;
			 for (k=0;k<3;k++) {
			 cout.width(12);
			 cout.precision(5);
			 cout<<CM[i].COOR(k)<<"  ";
			 cout.width(12);
			 cout.precision(5);
			 cout<<EA[i].COOR(k)<<"  ";
			 }
			 cout<<endl;*/
		}
	}
	else {
		cerr << "File " << filename << " could not be opened." << endl;
	}

	
	datafile.close();
	
	*nAtoms = nMol*N_MOL_ATOMS;
	*atomPos = new VECT[*nAtoms];
	*atomType = new char[*nAtoms];
	
	//Set up a generic water molecule
	
	VECT H1(3), H2(3), O(3), M(3);
	VECT H1_mol(3), H2_mol(3), O_mol(3), M_mol(3);
	
	//Set intitial water molecule atom locations assuming the WFF is aligned with the SFF
	H1.COOR(0) = H1_X;
	H1.COOR(1) = H1_Y;
	H1.COOR(2) = H1_Z;
	
	H2.COOR(0) = H2_X;
	H2.COOR(1) = H2_Y;
	H2.COOR(2) = H2_Z;
	
	O.COOR(0) = O_X;
	O.COOR(1) = O_Y;
	O.COOR(2) = O_Z;
	
	M.COOR(0) = M_X;
	M.COOR(1) = M_Y;
	M.COOR(2) = M_Z;
	
	for (i=0; i<nMol; i++) {
		//Rotate water molecules using the Euler angles
		MAT rotMat(3,3);
		double EA_array[3];
		
		for (j=0; j<3; j++) {
			EA_array[j] = EA[i].COOR(j);
		}
		
		rotMat = rotMat3D(EA_array);
		
		O_mol = rotMat*O;
		H1_mol = rotMat*H1;
		H2_mol = rotMat*H2;
		M_mol = rotMat*M;
		
		/*
		 cout << "O: ";
		 for (j=0; j<3; j++) {
		 cout << O_mol.COOR(j) << " ";
		 }
		 cout << endl;
		 
		 cout << "H1: ";
		 for (j=0; j<3; j++) {
		 cout << H1_mol.COOR(j) << " ";
		 }
		 cout << endl;
		 
		 cout << "H2: ";
		 for (j=0; j<3; j++) {
		 cout << H2_mol.COOR(j) << " ";
		 }
		 cout << endl;
		 
		 cout << "M: ";
		 for (j=0; j<3; j++) {
		 cout << M_mol.COOR(j) << " ";
		 }
		 cout << endl;
		 */
		
		
		//Add centre of mass to each atom to translate molecule
		
		O_mol = O_mol + NM_PER_ANG*CM[i];
		H1_mol = H1_mol + NM_PER_ANG*CM[i];
		H2_mol = H2_mol + NM_PER_ANG*CM[i];
		M_mol = M_mol + NM_PER_ANG*CM[i];
		
		/*
		 cout << "O: ";
		 for (j=0; j<3; j++) {
		 cout << O_mol.COOR(j) << " ";
		 }
		 cout << endl;
		 
		 cout << "H1: ";
		 for (j=0; j<3; j++) {
		 cout << H1_mol.COOR(j) << " ";
		 }
		 cout << endl;
		 
		 cout << "H2: ";
		 for (j=0; j<3; j++) {
		 cout << H2_mol.COOR(j) << " ";
		 }
		 cout << endl;
		 
		 cout << "M: ";
		 for (j=0; j<3; j++) {
		 cout << M_mol.COOR(j) << " ";
		 }
		 cout << endl;
		 */
		
		
		(*atomPos)[N_MOL_ATOMS*i + 0].DIM(3);
		(*atomPos)[N_MOL_ATOMS*i + 1].DIM(3);
		(*atomPos)[N_MOL_ATOMS*i + 2].DIM(3);
		(*atomPos)[N_MOL_ATOMS*i + 3].DIM(3);
		
		(*atomPos)[N_MOL_ATOMS*i + 0] = O_mol;
		(*atomPos)[N_MOL_ATOMS*i + 1] = H1_mol;
		(*atomPos)[N_MOL_ATOMS*i + 2] = H2_mol;
		(*atomPos)[N_MOL_ATOMS*i + 3] = M_mol;
		
		
		(*atomType)[N_MOL_ATOMS*i + 0] = 'O';
		(*atomType)[N_MOL_ATOMS*i + 1] = 'H';
		(*atomType)[N_MOL_ATOMS*i + 2] = 'H';
		(*atomType)[N_MOL_ATOMS*i + 3] = 'M';
		
		/*
		cout << "O: ";
		for (j=0; j<3; j++) {
			cout << atomPos[N_MOL_ATOMS*i + 0].COOR(j) << " ";
		}
		cout << endl;
		
		cout << "H1: ";
		for (j=0; j<3; j++) {
			cout << atomPos[N_MOL_ATOMS*i + 1].COOR(j) << " ";
		}
		cout << endl;
		
		cout << "H2: ";
		for (j=0; j<3; j++) {
			cout << atomPos[N_MOL_ATOMS*i + 2].COOR(j) << " ";
		}
		cout << endl;
		
		cout << "M: ";
		for (j=0; j<3; j++) {
			cout << atomPos[N_MOL_ATOMS*i + 3].COOR(j) << " ";
		}
		cout << endl;
		 */
	}
	
	/*
	//Output data to test in VMD
	ofstream outdatafile;
	outdatafile.open("unitcell.xyz", ios::out);
	if (outdatafile.is_open()) {
		outdatafile << (*nAtoms) << endl;
		outdatafile << "# Comment" << endl;
		for (i=0; i<(*nAtoms); i++) {
			outdatafile << atomType[i] << " ";
			for (j=0; j<3; j++) {
				outdatafile << atomPos[i].COOR(j)/NM_PER_ANG << " ";
			}
			outdatafile << endl;
			
		}
	}
	
	outdatafile.close();
	*/
	
	delete [] EA;
	delete [] CM;
}

double Q_TIP4P_Eng(VECT pos, double q, char *atomType, VECT *atomPos, int nAtoms) {
	int i;
	double Eng;
	
	VECT diffVec;
	double distance;
	
	Eng = 0;
	
	for (i=0; i<nAtoms; i++) {
		switch (atomType[i]) {
			case 'H':
				Eng += CoulombEng(pos, q, atomPos[i], Q_H);
				//cout << "H" << endl;
				break;
			case 'M':
				Eng += CoulombEng(pos, q, atomPos[i], Q_M);
				//cout << "M" << endl;
				break;
			case 'O':
				diffVec.DIM(pos.NUM_DIM);
				diffVec = atomPos[i] - pos;
				distance = diffVec.eNorm();
				if (distance>C_CUT_OFF) {
					i += 3;
					//cout << "Skipping... distance:" << distance << endl;
				}
				
				//cout << "O" << endl;
				break;
			default:
				cerr << "Atom type not recognized while executing Q_TIP4P_Eng()." << endl;
				exit(1);
				break;
		}
		//cout << Eng << " kJ/mol" << endl;
	}
	
	return Eng;
}

double LJ_TIP4P_Eng(VECT pos, char *atomType, VECT *atomPos, int nAtoms) {
	int i;
	double Eng;
	
	Eng = 0;
	
	for (i=0; i<nAtoms; i++) {
		switch (atomType[i]) {
			case 'H':
				break;
			case 'M':
				break;
			case 'O':
				Eng += LJEng(pos, atomPos[i], EPS_OH2, SIGMA_OH2);
				break;
			default:
				cerr << "Atom type not recognized while executing LJ_TIP4P_Eng()." << endl;
				exit(1);
				break;
		}
		//cout << Eng << " kJ/mol" << endl;
	}
	
	return Eng;
}

double LJ_TIP4P_Eng_Fast(VECT pos, char *atomType, VECT *atomPos, int nAtoms) {
	int i;
	double Eng;
	
	Eng = 0;
	
	for (i=0; i<nAtoms; i++) {
		switch (atomType[i]) {
			case 'H':
				break;
			case 'M':
				break;
			case 'O':
				Eng += LJEngFast(pos, atomPos[i], A_LJ_OH, B_LJ_OH);
				break;
			default:
				cerr << "Atom type not recognized while executing LJ_TIP4P_Eng()." << endl;
				exit(1);
				break;
		}
		//cout << Eng << " kJ/mol" << endl;
	}
	
	return Eng;
}

// To run: ./TIP4P_AH2_EngRoutinesTest_Periodic [water_molecule.xyz file]
int main(int argc, char **argv) {
	string filename, posFile, *cageType;
	string line;
	char *atomType;
	VECT *atomPos, pos, *cagePos;
	int nAtoms;
	double Eng, C_Eng, LJ_Eng, *C_engArray, *cageEng;
	int i, j, k, l, m, n, o, replicates, replicateSize, repSize[3], nCages;
	
	int chunk, nthreads, thread;
	
	filename = argv[1];
	posFile = argv[2];
	replicates = atoi(argv[3]);
	replicateSize = replicates*replicates*replicates; //Total number of replicates (includes original cell)
	C_engArray = new double [replicateSize];
	
	getTIP4Patoms(&atomType, &atomPos, &nAtoms, filename);
	
	pos.DIM(3);
	
	
	/*
	// Small Cage 1
	pos.COOR(0) = 0.6491; //nm
	pos.COOR(1) = 0.2165;
	pos.COOR(2) = 0.6492;
	*/
	
	/*
	// Small Cage 4
	pos.COOR(0) = 1.0824; //nm
	pos.COOR(1) = 0.2164;
	pos.COOR(2) = 1.0819;
	 */
	
	// Large Cage 1
	pos.COOR(0) = 0.8651; //nm
	pos.COOR(1) = 0.8655;
	pos.COOR(2) = 0.8654;
	
	ifstream datafile;
	
	datafile.open(posFile.c_str(), ios::in);
	if (datafile.is_open()) {
		datafile>>nCages;
		//cout<<nCages<<endl;		
		getline (datafile, line);//Eat up end of line character from first line
		
		// get rid of comment line
		getline (datafile, line);
		//cout<<line<<endl;
		
		cagePos = new VECT[nCages]; //Values in Angstroms
		cageType = new string [nCages];
		
		char name[12];
		for (i=0;i<nCages;i++) {

			cagePos[i].DIM(3);
			
			datafile>>cageType[i];
			//cout<<name<<"  ";
			datafile>>cagePos[i].COOR(0); 		 		
			datafile>>cagePos[i].COOR(1); 		
			datafile>>cagePos[i].COOR(2); 		

			cagePos[i].COOR(0) *= NM_PER_ANG; 		 		
			cagePos[i].COOR(1) *= NM_PER_ANG; 		
			cagePos[i].COOR(2) *= NM_PER_ANG;
			/*
			 cout << scientific;
			 for (k=0;k<3;k++) {
			 cout.width(12);
			 cout.precision(5);
			 cout<<CM[i].COOR(k)<<"  ";
			 cout.width(12);
			 cout.precision(5);
			 cout<<EA[i].COOR(k)<<"  ";
			 }
			 cout<<endl;*/
		}
	}
	else {
		cerr << "File " << filename << " could not be opened." << endl;
	}
	
	
	datafile.close();
	
	cageEng = new double [nCages];
	
for (m=0; m<nCages; m++) {

	Eng = 0.0; //kJ/mol
	C_Eng = 0.0;
	LJ_Eng = 0.0;
	
	cout << "Cage " << cageType[m] << ": ";
	
	
	
	for (l=0; l<3; l++) {
		pos.COOR(l) = cagePos[m].COOR(l);
		pos.COOR(l) += ((int) (replicates/2))*SYS_SIZE;
	}
	
	//for (m=0; m<2; m++) {
//		for (n=0; n<2; n++) {
//			for (o=0; 0<2; o++) {
	
	if (replicates > CHUNK_RATIO) {
		chunk = (int) (replicates/CHUNK_RATIO);
	}
	else if (replicates > (CHUNK_RATIO/10)){
		chunk = (int) (replicates/(CHUNK_RATIO/10));
	}
	else {
		chunk = 1;
	}
	repSize[0] = replicates;
	repSize[1] = replicates;
	repSize[2] = replicates;
	
#pragma omp parallel default(shared) private (i,j,k,l,m,n,o,thread)
	{
	thread = omp_get_thread_num();
	if (thread == 0) {
		nthreads = omp_get_num_threads();
		//cout << "Parallel Section Number of Threads: " << nthreads << endl;
	}
	#pragma omp for schedule(dynamic, chunk)  //Possibly try: http://software.intel.com/en-us/articles/openmp-loop-collapse-directive
				for (i=0; i<replicates; i++) {
					for (j=0; j<replicates; j++) {
						for (k=0; k<replicates; k++) {
							
							int index[3], displacement;
							VECT *atomPos2, displ;
							
							atomPos2 = new VECT [nAtoms];
							
							
							displ.DIM(3);
							
							//displ.COOR(0) = ((-1)^m)*i*SYS_SIZE;
//							displ.COOR(1) = ((-1)^n)*j*SYS_SIZE;
//							displ.COOR(2) = ((-1)^o)*k*SYS_SIZE;
							
							displ.COOR(0) = i*SYS_SIZE;
							displ.COOR(1) = j*SYS_SIZE;
							displ.COOR(2) = k*SYS_SIZE;
							
							
							
							
							
							for (l=0; l<nAtoms; l++) {
								atomPos2[l].DIM(3);
								atomPos2[l] = atomPos[l] + displ;
							}
							
											index[0] = i;
											index[1] = j;
											index[2] = k;
							

							
											displacement = disp(index, 3, repSize);
											C_engArray[displacement] = Q_TIP4P_Eng(pos, 1.0, atomType, atomPos2, nAtoms);
							
							
							//C_Eng += Q_TIP4P_Eng(pos, 1.0, atomType, atomPos2, nAtoms);
							
							delete [] atomPos2;
							
							
						} 
					}
				}
//			}
//		}
//	}
	}

	
	for (i=0; i<replicateSize; i++) {
		C_Eng += C_engArray[i];
	}
	cout << setprecision(6) << C_Eng << " kJ/mol" << endl;
	
}
	//cout <<	displ.COOR(0) << endl;
	//				cout <<	displ.COOR(1) << endl;
	//				cout <<	displ.COOR(2) << endl;
	//cout <<	atomPos2[l].COOR(0) << endl;
	//					cout <<	atomPos2[l].COOR(1) << endl;
	//					cout <<	atomPos2[l].COOR(2) << endl;
	//				index[0] = i;
	//				index[1] = j;
	//				index[2] = k;
	//				
	//				repSize[0] = replicates;
	//				repSize[1] = replicates;
	//				repSize[2] = replicates;
	//				
	//				
	//				displacement = disp(index, 3, repSize);
	//				C_engArray[displacement] = Q_TIP4P_Eng(pos, -1.0, atomType, atomPos2, nAtoms);
	
	//cout << "Coulomb Potential: " << C_Eng << " kJ/mol" << endl; //Should be -30.2331 kJ/mol for a 1.0e charge at the centre of the sII small cage 4 (1.0824nm, 0.2164nm, 1.0819nm) in the sII unitcell (sII_unitCelltest_CM_EA_ZYZ.xyz)
	
	//LJ_Eng += LJ_TIP4P_Eng_Fast(pos, atomType, atomPos, nAtoms);
//	
//	cout << "LJ: " << LJ_Eng << " kJ/mol" << endl; 
//	
//	Eng += C_Eng;
//	Eng += LJ_Eng;
//	
//	cout << "Coulomb+LJ: " << Eng << " kJ/mol" << endl; 
	
	

	
	
	
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
	 
	 
	
	delete [] atomPos;
	delete [] atomType;
}
