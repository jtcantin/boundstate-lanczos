#ifndef INTERFACE_H
#define	INTERFACE_H

#include <iostream>
#include "interfaceContainers.h"
#include "lanczosUnits.h"

//Vector and Matrix Classes
#include "vectClass.h" //This is to select which linear algebra class should be used.

//Linear Rotor Cartesian Position and Spherical Harmonics Rotational Basis Hv
#include "linRotCartSphHvRoutines.h"
#include "linRotCartSphContainers.h"

//Linear Rotor Cartesian Position and Spherical Harmonics Rotational Basis Hv with the Potential Fully Pre-calculated
#include "linRotCartSphHvFullVPrepRoutines.h"

//Linear Rotor Cartesian Position and Spherical Harmonics Rotational Basis Hv with the Potential Fully Pre-calculated and an angular continuous grid
#include "linRotCartSphHvContGridRoutines.h"

//Alavi Hydrogen Model
#include "Alavi_H2_Routines.h"

//Alavi Hydrogen Model - Continuous Grid
#include "Alavi_ContGrid_H2_Routines.h"

//TIP4P Model
#include "TIP4P_AH2_EngRoutines.h"

//SPC/E Model
#include "SPCE_AH2_EngRoutines.h"

//Coulomb Potential
#include "Coulomb_Routines.h"

//Isotropic Harmonic Oscillator Potential
#include "IsoHarmOsc_Routines.h"


void HvInterfaceSetup(string HvCalculatorSwitch, generalStor **general_data, void (**HvPrepPtr)(int, char**, generalStor*, lanczosStor*), void (**HvPtr)(int, char**, generalStor*, lanczosStor*, double*, double*));

#endif