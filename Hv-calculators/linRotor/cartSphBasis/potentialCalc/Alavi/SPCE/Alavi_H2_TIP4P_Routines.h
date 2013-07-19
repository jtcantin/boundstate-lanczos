#ifndef ALAVI_H2_TIP4P_ROUTINES_H
#define	ALAVI_H2_TIP4P_ROUTINES_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <string.h>
#include "gridFcns.h"
#include "energyRoutines.h"
#include "TIP4P_AH2_EngRoutines.h"
#include "rotmat.h"
#include "vectClass.h" //This is to select which linear algebra class should be used.
#include "boundStateContainers.h"

//#define VCEIL 100 //kJ/mol

void Alavi_TIP4P_point_Eng(double *CMpotential, double *Hpotential, universeProp *point_universe, sysAtoms *atomGeo);

#endif
