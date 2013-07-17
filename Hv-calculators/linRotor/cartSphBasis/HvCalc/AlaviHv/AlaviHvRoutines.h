#ifndef ALAVIHVROUTINES_H
#define	ALAVIHVROUTINES_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <string.h>
#include <iomanip>
#include <cfloat>
#include <omp.h>
#include "lanczosUnits.h"
#include "vectClass.h"
#include "boundStateContainers.h"
#include "Alavi_H2_Routines.h"
#include "Alavi_H2_TIP4P_Routines.h"
#include "generalHvRoutines.h"

double* calc_ulm(double x, double y, double z, double *v_lpmp, interfaceStor *interface, int rangeFlag);

void HvPrep_Internal(int argc, char **argv, interfaceStor *interface, lanczosStor *lanczos);

#endif