#include <iostream>
#include <cmath>
#include <cstdlib>
#include "vectClass.h" //This is to select which linear algebra class should be used.
#include "lanczosUnits.h"

#ifndef ENERGYROUTINES_H
#define	ENERGYROUTINES_H

double CoulombEng(const VECT& q1, double q1_q, const VECT& q2, double q2_q);

double LJEng(const VECT& a1, const VECT& a2, double epsilon, double sigma);

double LJEngFast(const VECT& a1, const VECT& a2, double A, double B);

#endif
