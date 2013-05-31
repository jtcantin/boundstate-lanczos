#include <iostream>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include "Cvectors.h"

using namespace std;

#define PI 3.14159265359 //From MMTK
#define EPS0 0.000572765638445 //From MMTK

double CoulombEng(const C2vect& q1, double q1_q, const C2vect& q2, double q2_q)
{
	double distance, k, Eng;
	C2vect diffVec;
	
	k = 1/(4*PI*EPS0);
	
	diffVec = q2 - q1;
	
	distance = diffVec.eNorm();
	
	Eng = k*q1_q*q2_q/distance;
	
	return Eng;
}
	
	

int main (int argc, char** argv)
{
	double box_x_max, box_y_max;
	int box_x_num, box_y_num;
	C2vect sourceCharge, testCharge;
	double sourceQ, testQ;
	
	double Eng;
	
	//MPI Variables
	int rank, size;
	MPI::Status status;
	
	if (argc != (4+1)) {
		cout << "Invalid number of arguments, need 4" << endl;
		exit (1);
	}
	
	//Initialize MPI
	MPI::Init(argc, argv);
	rank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();
	
	// Get universe size
	box_x_max = atof(argv[1]); //nm
	box_x_num = atoi(argv[2]);
	
	box_y_max = atof(argv[3]); //nm
	box_y_num = atoi(argv[4]);
	
	// Set Charges
	sourceQ = 1; //e
	testQ = 1; //e
	
	sourceCharge.setCor(0.0, 0.0);
	
	testCharge.setCor(1.0,1.0);
	
	Eng = CoulombEng (sourceCharge, sourceQ, testCharge, testQ);
	
	cout << "System Energy: " << Eng << " kJ/mol" << endl;
	
	
	
	MPI::Finalize();	
	return 0;
}
