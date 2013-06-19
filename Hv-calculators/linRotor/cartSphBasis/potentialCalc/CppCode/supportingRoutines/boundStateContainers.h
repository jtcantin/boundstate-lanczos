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

struct quadStor {
	double *GCabscissae;
	double *GCweights;
	int GCnum;
	double *GLabscissae;
	double *GLweights;
	int GLnum;
};

struct fiveDGrid {
	double *x_points;
	int nx;
	int x_max;
	double *y_points;
	int ny;
	int y_max;
	double *z_points;
	int nz;
	int z_max;
	double *theta_points;
	int ntheta;
	int theta_max;
	double *phi_points;
	int nphi;
	int phi_max;
};

struct pointPotentialStor {
	double *CMpotential;
	double *H_potential;
	universeProp *point_universe;
};

struct tesseralStor {
	double **L_lpmp; //[a][n]
	double **S_mp; //[b][m]
	double **L_lm; //[n][a]
	double **S_m; //[m][b]
};

#endif
