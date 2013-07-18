#ifndef INTERFACE_H
#define	INTERFACE_H

#include <iostream>
#include "lanczos.h"
#include "boundStateContainers.h"
#include "generalHvRoutines.h"
#include "Alavi_H2_Routines.h"

void HvInterfaceSetup(string HvCalculatorSwitch, generalStor **general_data, void (**HvPrepPtr)(int, char**, generalStor*, lanczosStor*), void (**HvPtr)(int, char**, generalStor*, lanczosStor*, double*, double*));

void Hv_Prep_linRotCartSph(int argc, char **argv, generalStor *general_data, lanczosStor *lanczos_data);

void Hv_linRotCartSph(int argc, char **argv, generalStor *general_data, lanczosStor *lanczos_data, double *vec, double *uec);


#endif