#ifndef ALAVI_H2_ROUTINES_H
#define	ALAVI_H2_ROUTINES_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <string.h>
#include "gridFcns.h"
#include "vectClass.h" //This is to select which linear algebra class should be used.
#include "boundStateContainers.h"
#include "lanczosUnits.h"

//#define VCEIL 100 //kJ/mol

//Alavi H2 parameters (NOTE: Lennard-Jones Parameters not included here as only mixed species are currently used and the LJ parameters are dependent on these other species)
#define Q_H2_H 0.4932 //e, from Alavi, S.; Ripmeester, J. A.; Klug, D. D. Molecular-dynamics study of structure II hydrogen clathrates. J. Chem. Phys. 2005, 123, 024507.
#define Q_H2_CM -0.9864 //e, from Alavi, S.; Ripmeester, J. A.; Klug, D. D. Molecular-dynamics study of structure II hydrogen clathrates. J. Chem. Phys. 2005, 123, 024507.

#define H2_H1_X 0.
#define H2_H1_Y 0.
#define H2_H1_Z 0.03707 //nm; Bondlength is 0.7414 Ang from Alavi_2005

#define H2_H2_X 0.
#define H2_H2_Y 0.
#define H2_H2_Z -0.03707 //nm; Bondlength is 0.7414 Ang from Alavi_2005

int AlaviMapH2GridToHGrid(H2_orient *H2_mol, universeProp *H_universe, double H_CM_dist);

void Alavi_H2_Eng(double *H2potential, double *CMpotential, double *H_potential, H2_orient *H2_mol, universeProp *point_universe);

double Alavi_H2_Eng_Point(interfaceStor *interface, H2_orient *H2_mol);

void Alavi_SiteSite_Eng(double *CMpotential, double *Hpotential, universeProp *point_universe, sysAtoms *atomGeo, interfaceStor *interface);

pointPotentialStorH2* preCalcPotential_Alavi(int numDim, double *gridMax, int *gridPoints, string geometryFilename, interfaceStor *interface);

#endif
