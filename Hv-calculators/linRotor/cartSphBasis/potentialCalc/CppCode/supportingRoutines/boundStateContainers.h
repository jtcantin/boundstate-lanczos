#include "vectClass.h" //This is to select which linear algebra class should be used.

#ifndef BOUNDSTATECONTAINERS_H
#define	BOUNDSTATECONTAINERS_H

#define CHUNK_RATIO 100 //How large the chunks should be relative to the total data size; 10 or 100 seem to be reasonable values.
#define NM_PER_BOHR 0.052918
#define ANG_PER_BOHR 0.52918

struct H2_orient {
	VECT *CM;
	double theta;
	double phi;
};

struct universeProp {
	VECT *grid;
	double *grid_max;
	double *d_i;
	int *grid_num;
	int numDim;
	int sysSize;
};

struct sysAtoms {
	char *atomType;
	VECT *atomPos;
	int nAtoms;
};

#endif
