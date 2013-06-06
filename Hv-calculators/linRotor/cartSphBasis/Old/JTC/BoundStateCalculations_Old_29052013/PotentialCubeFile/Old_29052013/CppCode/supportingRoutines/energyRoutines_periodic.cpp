#include <iostream>
#include <cmath>
#include <cstdlib>
#include "energyRoutines.h"

#include "TIP4P_AH2_EngRoutines.h"

using namespace std;

double CoulombEng(const VECT& q1, double q1_q, const VECT& q2, double q2_q){
	double distance, k, Eng;
	VECT diffVec;
	int j;
	
	k = 1/(4*PI*EPS0);
	
	diffVec.DIM(q1.NUM_DIM);
	
	//cout << q1.NUM_DIM << endl;
	//cout << q2.NUM_DIM << endl;
	
	diffVec = q2 - q1;
	
//#ifdef PERIODIC
//	for (j=0; j<q1.NUM_DIM; j++) {
//		if (diffVec.COOR(j) > SYS_SIZE*0.5) {
//			diffVec.COOR(j) -= SYS_SIZE;
//		}
//		else if (diffVec.COOR(j) < SYS_SIZE*0.5) {
//			diffVec.COOR(j) += SYS_SIZE;
//		}
//	}	
//#endif
	//Periodic BCs
//	for (j=0; j<q1.NUM_DIM; j++) {
//		diffVec.COOR(j) -= SYS_SIZE*(round(diffVec.COOR(j)/SYS_SIZE));
//	}
	
	distance = diffVec.eNorm();
	
	

//	if (distance<C_CUT_OFF) {
		Eng = k*q1_q*q2_q/distance;
//	}
//	else {
//		Eng = 0.0;
//		//cout << "Too far: Distance: " << distance << endl;
//	}

	
	return Eng;
}

double LJEng(const VECT& a1, const VECT& a2, double epsilon, double sigma){
	double distance, Eng;
	VECT diffVec;
	int j;

	diffVec.DIM(a1.NUM_DIM);

	diffVec = a2 - a1;
	
//	for (j=0; j<a1.NUM_DIM; j++) {
//		if (diffVec.COOR(j) > SYS_SIZE*0.5) {
//			diffVec.COOR(j) -= SYS_SIZE;
//		}
//		else if (diffVec.COOR(j) < SYS_SIZE*0.5) {
//			diffVec.COOR(j) += SYS_SIZE;
//		}
//	}

	//Periodic BCs
//	for (j=0; j<a1.NUM_DIM; j++) {
//		diffVec.COOR(j) -= SYS_SIZE*(round(diffVec.COOR(j)/SYS_SIZE));
//	}
		
	distance = diffVec.eNorm();

	Eng = 4*epsilon*(pow((sigma/distance),12) - pow((sigma/distance),6));
	
	return Eng;
}

double LJEngFast(const VECT& a1, const VECT& a2, double A, double B){
	double distance, Eng;
	VECT diffVec;
	int j;
	
	diffVec.DIM(a1.NUM_DIM);
	
	diffVec = a2 - a1;
	
	//	for (j=0; j<a1.NUM_DIM; j++) {
	//		if (diffVec.COOR(j) > SYS_SIZE*0.5) {
	//			diffVec.COOR(j) -= SYS_SIZE;
	//		}
	//		else if (diffVec.COOR(j) < SYS_SIZE*0.5) {
	//			diffVec.COOR(j) += SYS_SIZE;
	//		}
	//	}
	
	//Periodic BCs
//	for (j=0; j<a1.NUM_DIM; j++) {
//		diffVec.COOR(j) -= SYS_SIZE*(round(diffVec.COOR(j)/SYS_SIZE));
//	}
	
	distance = diffVec.eNorm();
	
	Eng = (A/pow(distance,12)) - (B/pow(distance,6));

	return Eng;
}