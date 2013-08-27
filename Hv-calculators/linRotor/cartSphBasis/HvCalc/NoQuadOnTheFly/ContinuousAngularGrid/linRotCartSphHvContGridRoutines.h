#ifndef	LINROTCARTSPHHVCONTGRIDROUTINES_H
#define	LINROTCARTSPHHVCONTGRIDROUTINES_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <string.h>
#include <iomanip>
#include <cfloat>
#include <omp.h>
#include "vectClass.h" //This is to select which linear algebra class should be used.
#include "interfaceContainers.h"
#include "linRotCartSphContainers.h"
#include "lanczosUnits.h"
#include "linRotCartSphHvRoutines.h"

double* calc_Vlmlpmp_NoQuad_ContGrid(interfaceStor *interface);

#endif