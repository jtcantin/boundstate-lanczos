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

#ifndef TIP4P_AH2_ENGROUTINES_H
#define TIP4P_AH2_ENGROUTINES_H

//TIP4P parameters; Values from a calculation using the values from: Water models http://www.lsbu.ac.uk/water/models.html (accessed May 17, 2013).
#define Q_H 0.520 //e 
#define Q_M -1.040 //e 

#define O_X 0.
#define O_Y 0.
#define O_Z 0.00655538 //nm

#define H1_X -0.075695 //nm
#define H1_Y 0.
#define H1_Z -0.0520326 //nm

#define H2_X 0.075695 //nm
#define H2_Y 0.
#define H2_Z -0.0520326 //nm

#define M_X 0.
#define M_Y 0.
#define M_Z -0.00844462 //nm

#define N_MOL_ATOMS 4

//Alavi H2 parameters; It should be noted that Alavi actually used the SPC/E model, not TIP4P
#define EPS_OH2 0.4306 //kJ/mol; Lorentz-Bertholet Combination of EPS_O-O = 0.6502 kJ/mol and EPS_H2-H2 = 0.2852 kJ/mol; EPS_OH2 = sqrt(EPS_O-O*EPS_H2-H2) = 0.4306 kJ/mol; Alavi, S.; Ripmeester, J. A.; Klug, D. D. Molecular-dynamics study of structure II hydrogen clathrates. J. Chem. Phys. 2005, 123, 024507.
#define SIGMA_OH2 0.3101 //nm; Lorentz-Bertholet Combination of SIGMA_O-O = 3.166 Ang and SIGMA_H2-H2 = 3.038 Ang; SIGMA_OH2 = sqrt(SIGMA_O-O*SIGMA_H2-H2) = 3.101 Ang; Alavi, S.; Ripmeester, J. A.; Klug, D. D. Molecular-dynamics study of structure II hydrogen clathrates. J. Chem. Phys. 2005, 123, 024507.
#define A_LJ_OH 1.362e-06 //kJ.nm^12/mol; 4*EPS_OH2*SIGMA_OH2^12
#define B_LJ_OH 1.532e-3 //kJ.nm^6/mol; 4*EPS_OH2*SIGMA_OH2^6

#define Q_H2_H 0.4932 //e, from Alavi, S.; Ripmeester, J. A.; Klug, D. D. Molecular-dynamics study of structure II hydrogen clathrates. J. Chem. Phys. 2005, 123, 024507.
#define Q_H2_CM -0.9864 //e, from Alavi, S.; Ripmeester, J. A.; Klug, D. D. Molecular-dynamics study of structure II hydrogen clathrates. J. Chem. Phys. 2005, 123, 024507.

#define H2_H1_X 0.
#define H2_H1_Y 0.
#define H2_H1_Z 0.03707 //nm; Bondlength is 0.7414 Ang from Alavi_2005

#define H2_H2_X 0.
#define H2_H2_Y 0.
#define H2_H2_Z -0.03707 //nm; Bondlength is 0.7414 Ang from Alavi_2005

//#define NM_PER_ANG 0.1

void getTIP4Patoms(char **atomType, VECT **atomPos, int *nAtoms, string filename);

double AH2_TIP4P_Eng(VECT H2cm, double *H2angles, char *atomType, VECT *atomPos);

double Q_TIP4P_Eng(VECT pos, double q, char *atomType, VECT *atomPos, int nAtoms);

double LJ_TIP4P_Eng(VECT pos, char *atomType, VECT *atomPos, int nAtoms);

double LJ_TIP4P_Eng_Fast(VECT pos, char *atomType, VECT *atomPos, int nAtoms);


#endif