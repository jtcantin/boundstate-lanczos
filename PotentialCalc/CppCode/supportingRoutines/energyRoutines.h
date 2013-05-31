#include "vectClass.h" //This is to select which linear algebra class should be used.

#ifndef ENERGYROUTINES_H
#define	ENERGYROUTINES_H

//Uses the same unit system as MMTK (length - nm, time - ps, energy - kJ/mol)
#define PI 3.141592653589793 //From NumPy (16 sig figs)
#define EPS0 0.0005727656384448188 //From MMTK (16 sig figs)
//#define K_CONST 138.93548461110962 // 1/(4*PI*EPS0) (17 sig figs)
#define K_CONST 138.9354846111096 // 1/(4*PI*EPS0) (16 sig figs)

double CoulombEng(const VECT& q1, double q1_q, const VECT& q2, double q2_q);
double LJEng(const VECT& a1, const VECT& a2, double epsilon, double sigma);
double LJEngFast(const VECT& a1, const VECT& a2, double A, double B);


#endif
