#ifndef GENERALHVROUTINES_H
#define	GENERALHVROUTINES_H

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

//#define VCEIL 100 //kJ/mol

double factorial(int num);

void genIndices_lm(int l_max, int ***qNum, int *length, int ***index, int *dims);

void legendrePoly_zerothOrder(int l_max, double **legendre, double **legendreDeriv, double x);

double* legendreRoots(int l);

double* legendreWeights(int l, double *roots);

void gauleg(double x1,double x2,double *x,double *w,int n);

double plgndr(int l,int m,double x);

double* normAssocLegendrePoly(int **qNum, int length, double x);

double* tesseralTrigTerm(int **qNum, int length, double phi);

void tesseralHarmonicsTerms(int **qNum, int length, double **legendre, double **trig, double cosTheta, double phi);

void gaussChebyshev(int numPoints, double **abscissae, double **weights);

void gaussLegendre(int numPoints, double **abscissae, double **weights);

void cartKinGrid(double x_max, int nPoints, double totalMass, double **kinMat, double **grid);

double* rotKinEng(int **qNum, int length, double momentOfInertia);

void tesseralTest(int l_max, int thetaPoints, int phiPoints);

double* calc_ulm(double x, double y, double z, double *v_lpmp, interfaceStor *interface, int rangeFlag);

void HvPrep_Internal(int argc, char **argv, interfaceStor *interface, lanczosStor *lanczos);

double* Mv_5D_oneCompositeIndex(double *v_ipjkn, double *mat_iip, int ni, int nj, int nk, int nn);

double* diagMv_5D_oneCompositeIndex(double *v_npijk, double *mat_n, int nn, int ni, int nj, int nk);

double* reshuffleIndices_5D_oneCompositeIndex(double *v_ijkn, int ni, int nj, int nk, int nn);

double* Tv_5D_oneCompositeIndex(interfaceStor *interface, double *v_ipjkn);

double* Vv_5D_oneCompositeIndex(interfaceStor *interface, double *v_ipjkn);

double* Hv_5D_oneCompositeIndex(interfaceStor *interface, double *v_ipjkn);

#endif