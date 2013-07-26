#ifndef SPCE_AH2_ENGROUTINES_H
#define SPCE_AH2_ENGROUTINES_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
//#include <cstdio>
#include "vectClass.h" //This is to select which linear algebra class should be used.
#include "rotmat.h"
#include "energyRoutines.h"
#include "lanczosUnits.h"
using namespace std;

//SPCE parameters; Values from a calculation using the values from: Water models http://www.lsbu.ac.uk/water/models.html (accessed May 17, 2013).
#define Q_H_SPCE 0.4238 //e from Alavi_2005
#define Q_O_SPCE -0.8476 //e from Alavi_2005

#define O_X_SPCE 0.
#define O_Y_SPCE 0.
#define O_Z_SPCE 0.006461408220 //nm

#define H1_X_SPCE -8.164965810E-2 //nm
#define H1_Y_SPCE 0.
#define H1_Z_SPCE -0.05127361870 //nm

#define H2_X_SPCE 8.164965810E-2 //nm
#define H2_Y_SPCE 0.
#define H2_Z_SPCE -0.05127361870 //nm

#define N_MOL_ATOMS_SPCE 3

//Alavi H2 parameters; It should be noted that Alavi actually used the SPC/E model, not TIP4P
#define EPS_OH2_SPCE 0.4306 //kJ/mol; Lorentz-Bertholet Combination of EPS_O-O = 0.6502 kJ/mol and EPS_H2-H2 = 0.2852 kJ/mol; EPS_OH2 = sqrt(EPS_O-O*EPS_H2-H2) = 0.4306 kJ/mol; Alavi, S.; Ripmeester, J. A.; Klug, D. D. Molecular-dynamics study of structure II hydrogen clathrates. J. Chem. Phys. 2005, 123, 024507.
#define SIGMA_OH2_SPCE 0.3101 //nm; Lorentz-Bertholet Combination of SIGMA_O-O = 3.166 Ang and SIGMA_H2-H2 = 3.038 Ang; SIGMA_OH2 = sqrt(SIGMA_O-O*SIGMA_H2-H2) = 3.101 Ang; Alavi, S.; Ripmeester, J. A.; Klug, D. D. Molecular-dynamics study of structure II hydrogen clathrates. J. Chem. Phys. 2005, 123, 024507.
#define A_LJ_OH_SPCE 1.364e-06 //kJ.nm^12/mol; 4*EPS_OH2*SIGMA_OH2^12
#define B_LJ_OH_SPCE 1.533e-3 //kJ.nm^6/mol; 4*EPS_OH2*SIGMA_OH2^6


//#define NM_PER_ANG 0.1

void getSPCEatoms(char **atomType, VECT **atomPos, int *nAtoms, string filename);

double Q_SPCE_Eng(VECT pos, double q, char *atomType, VECT *atomPos, int nAtoms);

double LJ_SPCE_Eng(VECT pos, char *atomType, VECT *atomPos, int nAtoms);

double LJ_SPCE_Eng_Fast(VECT pos, char *atomType, VECT *atomPos, int nAtoms);


#endif