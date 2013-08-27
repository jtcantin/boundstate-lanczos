#include "TIP4P_AH2_EngRoutines.h"

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
	
		getline (datafile, line);//Eat up end of line character from first line
		
		getline (datafile, line); // Get rid of comment line
		
		EA = new VECT[nMol]; //Values in radians
		CM = new VECT[nMol]; //Values in Angstroms
		
		char name[12];
		for (i=0;i<nMol;i++) {
			EA[i].DIM(3);
			CM[i].DIM(3);
			
			datafile>>name;

			datafile>>CM[i].COOR(0); 		
			datafile>>EA[i].COOR(0); 		
			datafile>>CM[i].COOR(1); 		
			datafile>>EA[i].COOR(1); 		
			datafile>>CM[i].COOR(2); 		
			datafile>>EA[i].COOR(2); 
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
				
		//Add centre of mass to each atom to translate molecule
		
		O_mol = O_mol + NM_PER_ANG*CM[i];
		H1_mol = H1_mol + NM_PER_ANG*CM[i];
		H2_mol = H2_mol + NM_PER_ANG*CM[i];
		M_mol = M_mol + NM_PER_ANG*CM[i];
		
		//Store atom positions
		(*atomPos)[N_MOL_ATOMS*i + 0].DIM(3);
		(*atomPos)[N_MOL_ATOMS*i + 1].DIM(3);
		(*atomPos)[N_MOL_ATOMS*i + 2].DIM(3);
		(*atomPos)[N_MOL_ATOMS*i + 3].DIM(3);
		
		(*atomPos)[N_MOL_ATOMS*i + 0] = O_mol;
		(*atomPos)[N_MOL_ATOMS*i + 1] = H1_mol;
		(*atomPos)[N_MOL_ATOMS*i + 2] = H2_mol;
		(*atomPos)[N_MOL_ATOMS*i + 3] = M_mol;
		
		//Store atom types
		(*atomType)[N_MOL_ATOMS*i + 0] = 'O';
		(*atomType)[N_MOL_ATOMS*i + 1] = 'H';
		(*atomType)[N_MOL_ATOMS*i + 2] = 'H';
		(*atomType)[N_MOL_ATOMS*i + 3] = 'M';
	}
	
	delete [] EA;
	delete [] CM;
}

double Q_TIP4P_Eng(VECT pos, double q, char *atomType, VECT *atomPos, int nAtoms) {
	int i;
	double Eng;
	
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
				break;
			default:
				cerr << "Atom type not recognized while executing Q_TIP4P_Eng()." << endl;
				exit(1);
				break;
		}
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
	}
	
	return Eng;
}


//int main(int argc, char **argv) {
//	string filename;
//	char *atomType;
//	VECT *atomPos, pos;
//	int nAtoms;
//	double Eng, C_Eng, LJ_Eng;
//	int i, j;
//	
//	//int i, j;
//	
//	filename = argv[1];
//	
//	getTIP4Patoms(&atomType, &atomPos, &nAtoms, filename);
//	
//	Eng = 0.0; //kJ/mol
//	C_Eng = 0.0;
//	LJ_Eng = 0.0;
//	
//	pos.DIM(3);
//	
//	/*
//	// Small Cage 1
//	pos.COOR(0) = 0.6491; //nm
//	pos.COOR(1) = 0.2165;
//	pos.COOR(2) = 0.6492;
//	*/
//	
//	/*
//	// Small Cage 4
//	pos.COOR(0) = 1.0824; //nm
//	pos.COOR(1) = 0.2164;
//	pos.COOR(2) = 1.0819;
//	 */
//	
//	// Large Cage 1
//	pos.COOR(0) = 0.8651; //nm
//	pos.COOR(1) = 0.8655;
//	pos.COOR(2) = 0.8654;
//	
//	
//	C_Eng += Q_TIP4P_Eng(pos, 1.0, atomType, atomPos, nAtoms);
//	
//	cout << "Coulomb: " << C_Eng << " kJ/mol" << endl; //Should be -30.2331 kJ/mol for a 1.0e charge at the centre of the sII small cage 4 (1.0824nm, 0.2164nm, 1.0819nm) in the sII unitcell (sII_unitCelltest_CM_EA_ZYZ.xyz)
//	
//	LJ_Eng += LJ_TIP4P_Eng_Fast(pos, atomType, atomPos, nAtoms);
//	
//	cout << "LJ: " << LJ_Eng << " kJ/mol" << endl; 
//	
//	Eng += C_Eng;
//	Eng += LJ_Eng;
//	
//	cout << "Coulomb+LJ: " << Eng << " kJ/mol" << endl; 
//	
//	
//
//	
//	
//	
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
//	 
//	 
//	
//	delete [] atomPos;
//	delete [] atomType;
//}
