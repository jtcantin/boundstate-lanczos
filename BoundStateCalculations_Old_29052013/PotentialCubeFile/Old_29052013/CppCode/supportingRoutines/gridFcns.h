#include "vectClass.h" //This is to select which linear algebra class should be used.
#include <string.h>
using namespace std;

#ifndef GRIDFCNS_H
#define GRIDFCNS_H

int disp(int *index, int dim, int *dim_size);

int indexDim(int dim, int disp, int dim_num, int *dim_size);

void writeCubeFile(int nDim, int *dimSize, VECT *dimGenerator, double *data, int dataSize, int nAtoms, int *atomNum, double *atomCharge, VECT *atomPos, const VECT& origin, string *comments, int unitFlag, string filename);

#endif