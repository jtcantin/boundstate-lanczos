#include "vectClass.h" //This is to select which linear algebra class should be used.

#ifndef ENERGYROUTINES_H
#define	ENERGYROUTINES_H

#define PI 3.14159265359 //From MMTK
#define EPS0 0.000572765638445 //From MMTK
#define K_CONST 13.89354846 // 1/(4*PI*EPS0)

double CoulombEng(const VECT& q1, double q1_q, const VECT& q2, double q2_q);
double LJEng(const VECT& a1, const VECT& a2, double epsilon, double sigma);
double LJEngFast(const VECT& a1, const VECT& a2, double A, double B);


#endif
