#ifndef ALAVI_CONTGRID_H2_ROUTINES_H
#define	ALAVI_CONTGRID_H2_ROUTINES_H

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

void Alavi_SiteSite_Eng_contGrid(double *H2potential, double *H2potentialPI, universeProp *point_universe, sysAtoms *atomGeo, interfaceStor *interface);

pointPotentialStorH2* preCalcPotential_Alavi_ContGrid(int numDim, double *gridMax, int *gridPoints, string geometryFilename, interfaceStor *interface, int argc, char** argv);

#endif
