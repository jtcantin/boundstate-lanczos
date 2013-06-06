#include "boundStateContainers.h"

#ifndef POTENTIALUNITS_H
#define	POTENTIALUNITS_H

#define CHUNK_RATIO 100 //How large the chunks should be relative to the total data size; 10 or 100 seem to be reasonable values.
#define NM_PER_BOHR 0.052918
#define ANG_PER_BOHR 0.52918

int AlaviMapH2GridToHGrid(H2_orient *H2_mol, universeProp *H_universe, double H_CM_dist);
void Alavi_H2_Eng(double *H2potential, double *CMpotential, double *H_potential, H2_orient *H2_mol, universeProp *point_universe);
void Alavi_point_Eng(double *CMpotential, double *Hpotential, double *H2potential, universeProp *point_universe, sysAtoms *atomGeo);



#endif
