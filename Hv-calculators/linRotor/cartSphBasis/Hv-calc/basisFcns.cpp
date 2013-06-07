#include <iostream>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <string.h>

using namespace std;

#define m_shift(X) (X+l_max)

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
	//qNum is a map from n to l and m; ie. qNum[n][0] = l and qNum[n][1] = m.
	//index is a map of l and m to an n value; ie. index[l][m + l_max] = n. Note: l_max is added to m as the array index starts at 0, not -l_max.
	//l_max is the maximum l value to be used
	//length will store the number of n values (ie. the length of the first dimension of qNum, such that max_n = length-1)
	//dims stores the length of each dimension of index; ie. l_max+1 = dims[0] and (2*l_max+1) = dims[1]
	
	int l, m, i;
	//cout <<	((l_max +1)^2) << endl;
	*length = (int) pow(double(l_max+1),2);
	
	*qNum = new int* [*length];
	
	dims[0] = l_max + 1;
	dims[1] = 2*l_max + 1;
	*index = new int* [dims[0]];
	
	//Initialize index to -1 so that non-existent lm pairs return a -1 index.
	for (l=0; l<dims[0]; l++) {
		(*index)[l] = new int [dims[1]];
		for (m=0; m<dims[1]; m++) {
			(*index)[l][m] = -1;
		}
	}
	
	i = 0;
	for (m=(-1*l_max); m<=l_max; m++) {
		
		for (l=abs(m); l<=l_max; l++) {
			cout << i << " " <<	l << " " << m << endl;
			
			(*qNum)[i] = new int [2];
			(*qNum)[i][0] = l;
			(*qNum)[i][1] = m;

			(*index)[l][m_shift(m)] = i;
			
			i++;
		}
	}
	
	cout << endl;
	cout << "Run through m-values:" << endl;
	for (l=0; l<=l_max; l++) {
		for (m=(-1*l); m<=l; m++) {
			cout << (*index)[l][m + l_max] << " ";
		}
		cout << endl;
	}
	
	cout << endl;
	cout << "Run through index:" << endl;
	
	for (l=0; l<=l_max; l++) {
		for (m=0; m<dims[1]; m++) {
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
	
	l_max = atoi(argv[1]);
	
	cout <<	"l_max: " << l_max << endl;
	
	genIndices_lm(l_max, &qNum, &length, &index, dims);
	
	cout << "length: " << length << endl;
	
	cout << "qNum: 2" << endl;
	cout << "(l,m): (" << qNum[2][0] << "," << qNum[2][1] << ")" << endl;
	
	cout << "qNum: 0" << endl;
	cout << "(l,m): (" << qNum[0][0] << "," << qNum[0][1] << ")" << endl;
	
	cout << "qNum: 3" << endl;
	cout << "(l,m): (" << qNum[3][0] << "," << qNum[3][1] << ")" << endl;
	
	cout << "l m: 0 -1" << endl;
	cout << "n: " << index[0][m_shift(-1)] << endl;
	
	cout << "l m: 1 -1" << endl;
	cout << "n: " << index[1][m_shift(-1)] << endl;
	
	cout << "l m: 1 0" << endl;
	cout << "n: " << index[1][m_shift(0)] << endl;
	
	
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


