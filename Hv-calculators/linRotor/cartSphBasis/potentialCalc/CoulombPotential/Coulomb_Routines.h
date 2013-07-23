#ifndef COULOMB_ROUTINES_H
#define	COULOMB_ROUTINES_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <string.h>
#include "vectClass.h" //This is to select which linear algebra class should be used.
#include "interfaceContainers.h"
#include "linRotCartSphContainers.h"
#include "lanczosUnits.h"
#include "energyRoutines.h"

double CoulombPotential(interfaceStor *interface, H2_orient *lin_mol);

pointPotentialStorH2* preCalcPotential_Coulomb(int numDim, double *gridMax, int *gridPoints, string geometryFilename, interfaceStor *interface, int argc, char **argv);


#endif
