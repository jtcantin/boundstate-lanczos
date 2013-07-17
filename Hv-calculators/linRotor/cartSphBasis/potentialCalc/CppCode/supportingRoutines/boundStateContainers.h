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
	
	~universeProp() {
		delete [] grid;
		delete [] grid_max;
		delete [] d_i;
		delete [] grid_num;
	};
};

struct sysAtoms {
	char *atomType;
	VECT *atomPos;
	int nAtoms;
	
	~sysAtoms() {
		delete [] atomPos;
		delete [] atomType;
	};
};

struct quadStor {
	double *GCabscissae;
	double *GCweights;
	int GCnum;
	
	double *GLabscissae;
	double *GLweights;
	int GLnum;
	
	~quadStor() {
		delete [] GCabscissae;
		delete [] GCweights;
		
		delete [] GLabscissae;
		delete [] GLweights;
	};
};

struct fiveDGrid {
	double *x_Grid;
	int nx;
	int x_max;
	double *xKinMat;
	
	double *y_Grid;
	int ny;
	int y_max;
	double *yKinMat;
	
	double *z_Grid;
	int nz;
	int z_max;
	double *zKinMat;
	
	//double *theta_Grid;
	//int ntheta;
	
	//double *phi_Grid;
	//int nphi;
	
	~fiveDGrid() { //DO A BUNCH OF IF STATEMENTS ON THE POINTER VALUE
		if ((x_Grid==y_Grid) && (x_Grid==z_Grid)) {
			;
		}
		else if ((x_Grid==y_Grid) && (x_Grid!=z_Grid)) {
			delete [] z_Grid;
			delete [] zKinMat;
		}
		else if ((x_Grid!=y_Grid) && (y_Grid==z_Grid)) {
			delete [] y_Grid;
			delete [] yKinMat;
		}
		else if ((x_Grid!=y_Grid) && (y_Grid!=z_Grid) && (x_Grid==z_Grid)){
			delete [] y_Grid;
			delete [] yKinMat;
		}
		else if ((x_Grid!=y_Grid) && (y_Grid!=z_Grid) && (x_Grid!=z_Grid)){
			delete [] y_Grid;
			delete [] yKinMat;
			
			delete [] z_Grid;
			delete [] zKinMat;
		}
		
		delete [] x_Grid;
		delete [] xKinMat;
	};
};

struct lmFBR {
	int lmax;
	
	int **qNum;
	int length;
	
	int **index;
	int *dims;
	
	double *rotKinMat;
	
	~lmFBR() {
		int i;
		
		for (i=0; i<length; i++) {
			delete [] qNum[i];
		}
		delete [] qNum;
		
		for (i=0; i<=lmax; i++) {
			delete [] index[i];
		}
		delete [] index;
		
		delete [] dims;
		delete [] rotKinMat;
	};
};

struct pointPotentialStorH2 {
	double *CMpotential;
	double *H_potential;
	universeProp *potentialUniverse;
	double potentialCeiling;
	
	~pointPotentialStorH2() {
		delete [] CMpotential;
		delete [] H_potential;
		delete potentialUniverse;
	};
};

struct tesseralStor {
	int na;
	int nb;
	int lmax; //These three integers are here only to know the array sizes for the destructor
	
	double **L_lpmp; //[a][n]
	double **S_mp; //[b][m]
	double **L_lm; //[n][a]
	double **S_m; //[m][b]
	
	~tesseralStor() {
		int i;
		
		for (i=0; i<na; i++) {
			delete [] L_lpmp[i];
		}
		delete [] L_lpmp;
		
		for (i=0; i<nb; i++) {
			delete [] S_mp[i];
		}
		delete [] S_mp;
		
		for (i=0; i<((lmax+1)*(lmax+1)); i++) {
			delete [] L_lm[i];
		}
		delete [] L_lm;
		
		for (i=-lmax; i<=lmax; i++) {
			delete [] S_m[i + lmax];
		}
		delete [] S_m;
	};
};

struct interfaceStor {
	quadStor *quadrature;
	fiveDGrid *grids;
	pointPotentialStorH2 *potential;
	tesseralStor *tesseral;
	tesseralStor *tesseral2PI;
	lmFBR *lmBasis;
	
	/*
	interfaceStor() {
		quadrature = new *quadStor [1];
		grids = new *fiveDGrid [1];
		potential = new *pointPotentialStorH2 [1];
		tesseral = new *tesseralStor [1];
		tesseral2PI = new *tesseral2PI [1];
		lmBasis = new *lmFBR [1];
	};*/
	
	~interfaceStor() {
		delete quadrature;
		delete grids;
		delete potential;
		delete tesseral;
		delete tesseral2PI;
		delete lmBasis;
	};
};

struct lanczosStor {
	int total_basis_size;
	string sim_descr;
	string sim_descr_short;
};

struct generalStor {
	double dummy;
};

#endif
