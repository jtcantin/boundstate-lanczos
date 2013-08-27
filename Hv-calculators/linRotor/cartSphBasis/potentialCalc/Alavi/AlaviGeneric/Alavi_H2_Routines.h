#ifndef ALAVI_H2_ROUTINES_H
#define	ALAVI_H2_ROUTINES_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <string.h>
#include "gridFcns.h"
#include "vectClass.h" //This is to select which linear algebra class should be used.
#include "interfaceContainers.h"
#include "linRotCartSphContainers.h"
#include "lanczosUnits.h"
#include "Alavi_Parameters.h"

int AlaviMapH2GridToHGrid(H2_orient *H2_mol, universeProp *H_universe, double H_CM_dist);

void Alavi_H2_Eng(double *H2potential, double *CMpotential, double *H_potential, H2_orient *H2_mol, universeProp *point_universe);

double Alavi_H2_Eng_Point(interfaceStor *interface, H2_orient *H2_mol);

void Alavi_SiteSite_Eng(double *CMpotential, double *Hpotential, universeProp *point_universe, sysAtoms *atomGeo, interfaceStor *interface);

pointPotentialStorH2* preCalcPotential_Alavi(int numDim, double *gridMax, int *gridPoints, string geometryFilename, interfaceStor *interface, int argc, char** argv);

#endif
