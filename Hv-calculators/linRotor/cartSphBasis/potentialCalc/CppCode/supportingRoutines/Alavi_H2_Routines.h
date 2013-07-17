#ifndef ALAVI_H2_ROUTINES_H
#define	ALAVI_H2_ROUTINES_H

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
#include "potentialUnits.h"
#include "boundStateContainers.h"

#define VCEIL 100 //kJ/mol

int AlaviMapH2GridToHGrid(H2_orient *H2_mol, universeProp *H_universe, double H_CM_dist);
void Alavi_H2_Eng(double *H2potential, double *CMpotential, double *H_potential, H2_orient *H2_mol, universeProp *point_universe);
void Alavi_point_Eng(double *CMpotential, double *Hpotential, universeProp *point_universe, sysAtoms *atomGeo);
double Alavi_H2_Eng_Point(double *CMpotential, double *H_potential, H2_orient *H2_mol, universeProp *point_universe);

#endif
