#include <iostream>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <string.h>

using namespace std;

//Algorithm from http://www.livephysics.com/computational-physics/fortran/fortran-subroutines-functions/
int factorial(int num) {
	int fact = 1;
	int i;
	
	for (i=2; i<=num; i++) {
		fact *= i;
	}
	
	return fact;
}

void genIndices_lm(int l_max, int ***qNum, int *length, int ***index, int *dims){
	int l, m, i;
	//cout <<	((l_max +1)^2) << endl;
	*length = (int) pow(double(l_max+1),2);
	
	*qNum = new int* [*length];
	*index = new int* [l_max];
	
	for (l=0; l<=l_max; l++) {
		(*index)[l] = new int [(2*l_max + 1)];
		for (m=0; m<=(2*l_max + 1); m++) {
			(*index)[l][m] = -1;
		}
	}
	
	i = 0;
	for (m=(-1*l_max); m<=l_max; m++) {
		
		for (l=abs(m); l<=l_max; l++) {
			//cout << "(l, m): (" << l << ", " << m << ")" << endl;
			
			//cout << i << " " <<	l << " " << m << endl;
			
			(*qNum)[i] = new int [2];
			(*qNum)[i][0] = l;
			(*qNum)[i][1] = m;
			cout << int (m+round((2*l_max + 1)/2)) << endl;
			cout << (*index)[l][int (m+round((2*l_max + 1)/2))] << endl;
			(*index)[l][int (m+round((2*l_max + 1)/2))] = i;
			
			i++;
		}
	}
	
	int j;
	for (l=0; l<=l_max; l++) {
		for (m=(-1*l); m<=l; m++) {
			cout << (*index)[l][int (m+round((2*l_max + 1)/2))] << " ";
		}
		cout << endl;
	}
	
	for (l=0; l<=l_max; l++) {
		for (m=0; m<=(2*l_max + 1); m++) {
			cout << (*index)[l][m] << " ";
		}
		cout << endl;
	}
	
}

double* assocLegendrePoly(int **qNum, int length, double x){
	double *legendre;
	int l_max;
	
	l_max = qNum[length][0];
	
	return legendre;
}


void tesseralHarmonics(int **qNum, int length, double *legendre, double *trig, double theta, double phi) {
	int m, l;
}

int main(int argc, char** argv) {
	int l_max, length;
	int **qNum, **index, dims[2];
	int i;
	
	l_max = 4;
	
	genIndices_lm(l_max, &qNum, &length, &index, dims);
	
	cout << "length: " << length << endl;
	
	for (i=0 ; i<length; i++) {
		//cout << i << " " <<	qNum[i][0] << " " << qNum[i][1] << endl;
		
		delete [] qNum[i];
	}
	delete [] qNum;
	
	for (i=0; i<l_max; i++) {
		delete [] index[i];
	}
	delete [] index;
	
	
	return 0;
}


