#ifndef	LINROTCARTSPHHVFULLVPREPROUTINES_H
#define	LINROTCARTSPHHVFULLVPREPROUTINES_H

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

double* calc_Vlmlpmp_NoQuad(interfaceStor *interface);

void HvPrep_Internal_NoQuad(int argc, char **argv, interfaceStor *interface, lanczosStor *lanczos);

void quadratureConvergenceStudy_NoQuad(interfaceStor *interface, lanczosStor *lanczos);

double* Vv_5D_oneCompositeIndex_NoQuad(interfaceStor *interface, double *v_ijknp);

double* Hv_5D_oneCompositeIndex_NoQuad(interfaceStor *interface, double *v_ipjkn);

void Hv_Prep_linRotCartSph_NoQuad(int argc, char **argv, generalStor *general_data, lanczosStor *lanczos_data);

void Hv_linRotCartSph_NoQuad(int argc, char **argv, generalStor *general_data, lanczosStor *lanczos_data, double *vec, double *uec);

double* calcSym(interfaceStor *interface, string symFlag);

#endif