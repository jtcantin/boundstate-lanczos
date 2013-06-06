//This is to select which linear algebra class should be used.
#define VECT_CLASS 1

#if VECT_CLASS == 0
#include "Cvectors.h" //Written by Joshua Cantin
#define VECT CNvect
#define MAT mat2
#define COOR(x) coor[x]
#define DIM(x) setDim(x)
#define SET_COOR(x, y) setCoor(x ,y)
#define	NUM_DIM num_dim
#elif VECT_CLASS == 1
#include "matvec2.h" //Written by Dr. Pierre-Nicholas Roy
#define VECT vector
#define MAT matrix
#define COOR(x) TheVector[x]
#define DIM(x) setSz(x)
#define SET_COOR(x, y) setCoor(x)
#define	NUM_DIM Sz
#endif

/*
 Note: The folowing properties are the same in both classes:
 -eNorm
 */