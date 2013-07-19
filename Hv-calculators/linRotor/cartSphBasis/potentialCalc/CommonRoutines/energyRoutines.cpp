#include "energyRoutines.h"

using namespace std;

double CoulombEng(const VECT& q1, double q1_q, const VECT& q2, double q2_q){
	double distance, k, Eng;
	VECT diffVec;
	
	//k = 1/(4*PI*EPS0);
	
	diffVec.DIM(q1.NUM_DIM);
	
	//cout << q1.NUM_DIM << endl;
	//cout << q2.NUM_DIM << endl;
	
	diffVec = q2 - q1;
	
	distance = diffVec.eNorm();
	
	//Eng = k*q1_q*q2_q/distance;
	
	Eng = K_CONST*q1_q*q2_q/distance;
	
	return Eng;
}

double LJEng(const VECT& a1, const VECT& a2, double epsilon, double sigma){
	double distance, Eng;
	VECT diffVec;

	diffVec.DIM(a1.NUM_DIM);

	diffVec = a2 - a1;

	distance = diffVec.eNorm();

	Eng = 4*epsilon*(pow((sigma/distance),12) - pow((sigma/distance),6));
	
	return Eng;
}

double LJEngFast(const VECT& a1, const VECT& a2, double A, double B){
	double distance, Eng;
	VECT diffVec;
	
	diffVec.DIM(a1.NUM_DIM);
	
	diffVec = a2 - a1;
	
	distance = diffVec.eNorm();
	
	Eng = (A/pow(distance,12)) - (B/pow(distance,6));

	return Eng;
}