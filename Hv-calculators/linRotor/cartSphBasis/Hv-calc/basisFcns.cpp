#include <iostream>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <string.h>
#include <iomanip>
#include <cfloat>
#include "lanczosUnits.h"

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

//This function generates two maps for the lm basis. One, qNum, goes from n to l,m and the other, index, goes from l,m to n. 
//   The mapping order is based on method of partial summation to be used, which sums from m=-l_max to l_max as the outer loop and l = abs(m) to l_max
//  as the inner loop. i.e. first l is summed over, then m is.
void genIndices_lm(int l_max, int ***qNum, int *length, int ***index, int *dims){
	//l_max is the maximum l value to be used
	//qNum is a map from n to l and m; ie. qNum[n][0] = l and qNum[n][1] = m.
	//length will store the number of n values (ie. the length of the first dimension of qNum, such that max_n = length-1)
	//index is a map of l and m to an n value; ie. index[l][m + l_max] = n.  
	//		Note: l_max is added to m as the array index starts at 0, not -l_max. The m_shift(x) macro can do this shift automatically.
	//dims stores the length of each dimension of index; ie. dims[0] = l_max+1 and dims[1] = 2*l_max+1
	
	int l, m, i;

	//Get size of basis; length = (l_max + 1) ^ 2
	*length = (int) pow(double(l_max+1),2);	
	*qNum = new int* [*length];
	
	//Get dimensions of l,m to n map
	dims[0] = l_max + 1;
	dims[1] = 2*l_max + 1;
	*index = new int* [dims[0]];
	
	//Initialize index to -1 so that non-existent lm pairs (within the dimensions of index) return a -1 index.
	for (l=0; l<dims[0]; l++) {
		(*index)[l] = new int [dims[1]];
		for (m=0; m<dims[1]; m++) {
			(*index)[l][m] = -1;
		}
	}
	
	//Populate both maps
	i = 0;
	for (m=(-1*l_max); m<=l_max; m++) {
		for (l=abs(m); l<=l_max; l++) {
			
			//Add in l and m to qNum for a given n; n to l,m map
			(*qNum)[i] = new int [2];
			(*qNum)[i][0] = l;
			(*qNum)[i][1] = m;

			//Add in n for a given l,m; l,m to n map
			(*index)[l][m_shift(m)] = i;
			
			i++;
		}
	}	
}

//Calculation of the Legendre polynomials of order 0 and their derivatives
//Bonnet's Recursion formula is used (Abramowitz and Stegun pg 334, Eq. 8.5.3 and 8.5.4 where mu = 0)
//Valid for -1<x<1.
void legendrePoly_zerothOrder(int l_max, double **legendre, double **legendreDeriv, double x){
	int l;
	double ld;
	
	//Initialize containers
	(*legendre) = new double [l_max+1];
	(*legendreDeriv) = new double [l_max+1];
	
	(*legendre)[0] = 1.0;
	(*legendre)[1] = x;
	
	(*legendreDeriv)[0] = 0.0;
	
	//Calculate the polynomials
	for (l=2; l<=l_max; l++) {
		ld = double (l) - 1.0;
		(*legendre)[l] = (1.0/(ld + 1.0)) * ( (2.0*ld+1.0) * x * (*legendre)[l-1] - ld*(*legendre)[l-2]);
	}
	
	//Calculate the derivatives
	for (l=1; l<=l_max; l++) {
		ld = double (l);
		(*legendreDeriv)[l] = (1/(x*x-1)) * (ld * x * (*legendre)[l] - ld * (*legendre)[l-1]);
	}
}

//Calculation of the zeros of the Legendre polynomials of order 0 using Newton's method (not secant as the derivatives are known)
double* legendreRoots(int l, int *numRootsFound) {
	double *legendre, *legendreDeriv, *roots;
	double oldRoot, newRoot;
	double absError, absError2, max_absError, d_x;
	int i, j, numRoots, n, rootCount, found;
	
	int iteration, max_iterations;
	
	numRoots = l; //There are l roots for a Legendre polynomial of degree l.
	
	max_absError = DBL_EPSILON;
	absError = max_absError + 1000000;
	
	max_iterations = 10000;
	
	roots = new double [l];
	
	//Step through domain (-1 to 1) to find all unique roots
	n = l*100;
	d_x = 2.0/double(n);
	rootCount = 0;
	found = 0;
	iteration = 0;
	for (i=0; i<n; i++) {
		//Incrementally increase the initial guess to catch all roots; it is brute force, but it seems to work
		oldRoot = (i+1)*d_x-1.0;
		
		//Use the relation x_i ~= cos(pi*(i-1/4)/(l+1/2)) to get a more efficient guess
		//oldRoot = cos(PI * (double(i)-0.25) / (double(l)+1/2));
		
		iteration = 0; //Reset the iteration count
		
		//Use Newton's method to find the root (ie. x_(i+1) = x_i - f(x_i)/f'(x_i))
		while ((absError>max_absError) && (iteration<max_iterations)) {
			
			//Get f(x) and f'(x)
			legendrePoly_zerothOrder(l, &legendre, &legendreDeriv, oldRoot);
			
			//x_(i+1) = x_i - f(x_i)/f'(x_i)
			newRoot = oldRoot - legendre[l]/legendreDeriv[l];
			
			//Get the error
			absError = fabs(newRoot-oldRoot);
			
			//Update the root
			oldRoot = newRoot;
			iteration++;
		}
		
		if (iteration>=max_iterations) { //Root is not converging, so ignore it.
			continue;
		}
		
		absError = max_absError + 1000000; //Reset the error magnitude
		
		
		//If the root is nan or is outside of (-1,1), ignore it.
		if ((fabs(newRoot) >= 1.0)||isnan(newRoot)) { 
			continue;
		} 
		else {
			//Cycle through all current roots, checking if the newRoot has already been found.
			for (j=0; j<rootCount; j++) {
				absError2 = fabs(newRoot-roots[j]);
				if (absError2<DBL_EPSILON) { //Root has been already found if newRoot and roots[j] are within the machine epsilon
					found = 1;
					break;
				}
			}
			if (found == 0) { //It is a new root, so add it to the list
				roots[rootCount] = newRoot;
				rootCount++;
			}
			found = 0;
		}
	}
	//Sort the roots
	double storage;
	int swapped;
	
	swapped = 1;
	//Sort roots
	while (swapped == 1) {
		swapped = 0;
		for (i=0; i<(rootCount-1); i++) {
			if (roots[i] > roots[i+1]){
				storage = roots[i+1];
				roots[i+1] = roots[i];
				roots[i] = storage;
				swapped = 1;
			}
		}
	}
	
	*numRootsFound = rootCount;
	return roots;
}

// Dr. PN Roy's code for the calculation of the Gauss-Legendre abscissae and weights
void gauleg(double x1,double x2,double *x,double *w,int n)
{
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;
	
	m=(n+1)/2;
	xm=0.5*(x2+x1); //Mean x value
	xl=0.5*(x2-x1); //Length of domain
	
	//Cycle through each positive root
	for (i=1;i<=m;i++)  {
		
		z=cos(M_PI*((double)i-0.25)/((double)n+0.5)); //Make an initial guess that is efficient for Newton's Method
		
		//Perform Newton's method
		do {
			//Calculate the Legendre polynomial
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*(double)j-1.0)*z*p2-((double)j-1.0)*p3)/(double)j;
			}
			
			pp=(double)n*(z*p1-p2)/(z*z-1.0); //Calculate the derivative of the Legendre Polynomial
			
			z1=z;
			z=z1-p1/pp; //Newton's Method equation: x_(i+1) = x_i - f(x_i)/f'(x_i)
		} while (fabs(z-z1) > DBL_EPSILON);
		
		x[i-1]=xm-xl*z; //Store the negative root;
		x[n-i]=xm+xl*z; //Store the positive root;
		
		w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n-i]=w[i-1];
	}
}

//Calculation of the set of "spherical Legendre polynomials" for l = 0 to l_max and m = -l_max to l_max 
// that is the spherical harmonics excluding only the trigonometric term (i.e. cos(m*phi) and sin(|m|*phi) )
// This function returns a pointer to an array of the spherical Legendre polynomials
double* sphLegendrePoly(int **qNum, int length, double x){
	//Storage arrays
	double *legendre, **legenArr;
	
	//Calculation variables
	double sqrtTwoPiInv, recFactor, phaseFactor, legenRecFactor, partialNormFactor;
	int l_max;
	
	//Index variables
	int l, m, n;
	double ld, md;
	
	l_max = qNum[length - 1][0];
	
	//Calculation of the set of associated Legendre polynomials P_lm(x)
	// for l = 0, 1, ..., l_max and m = 0, 1, ..., l.
	// First the Legendre polynomials are calculated:
	//    P^m_l(x) = 1/(l!.2^l) * (1-x^2)^(m/2) * d^(l+m)/dx^(l+m) [(x^2-1)^l]
	// Then each is multiplied by:
	//    N_lm = (-1)^m * sqrt{ [(2l+1) * (l-m)!] / [4.PI*(l+m)!] }  =  (-1)^m * sqrt[1/(2.PI)]  * sqrt[ (l-m)! * (2l+1) / (2(l+m)!) ]
	// This returns P_lm(x).
	// Note: This code (and comment) is heavily based on the code by Toby Zeng in the file ylm_py2.f
	
	sqrtTwoPiInv = 1.0 / sqrt(2.0*PI);
	
	legenArr = new double* [l_max+1];
	
	for (l=0; l<=l_max; l++) {
		legenArr[l] = new double [l_max+1]; //l_max+1 and not 2l_max+1 as only m>=0 is used.
	}
	
	//First the Legendre Polynomials are going to be calculated
	
	//A recursion factor
	recFactor = sqrt(1.0-x*x);
	
	//Initial value; P^0_0 = 1.0
	legenArr[0][0] = 1.0;
	
	//Calculate the values when l = m (the diagonal elements of the 2D array) using a recursion formula	
	// Recursion formula: P^l_l = (2l-1)*P^(l-1)_(l-1)*sqrt(1-x^2)
	for (l=1; l<=l_max; l++) {
		ld = double(l);
		legenArr[l][l] = (2.0*ld-1.0) * legenArr[l-1][l-1] * recFactor; 
	}
	
	//Calculate the values for when m < l; calculate down each column (i.e. all l for a given m)
	for (m=0; m<=(l_max-1); m++) {
		for (l=m; l<=(l_max-1); l++) {
			ld = double(l);
			md = double(m);
			
			if ( (l-1) < m ) {
				legenRecFactor = 0.0; //If the previous l is above the diagonal (ie. l-1 < m), there is no contribution from this point.
			}
			else {
				legenRecFactor = legenArr[l-1][m]; //If the previous l has a value, use in the recursion formula.
			}
			
			//Recursion relation: P^(l+1)_m = [(2l+1)*x*P^l_m - (l+m)*P^(l-1)_m] / [l-m+1], where P^(l-1)_m = 0 if l-1 < m
			legenArr[l+1][m] = ( ( (2.0*ld+1.0) * x * legenArr[l][m] ) - ( (ld+md)*legenRecFactor ) ) / (ld-md+1.0);

		}
	}
	
	//At this point, the Legendre Polynomials have been calculated (except for a phase factor (-1)^m).
	//Now, the acquired polynomials are going to be renormalized, with:
	//    N_lm = (-1)^m * sqrt{ [(2l+1) * (l-m)!] / [4.PI*(l+m)!] }  =  (-1)^m * sqrt[1/(2.PI)]  * sqrt[ (l-m)! * (2l+1) / (2(l+m)!) ]
	
	for (l=0; l<=l_max; l++) {
		phaseFactor = 1.0;
		for (m=0; m<=l; m++) {
			ld = double(l);
			md = double(m);
			
			partialNormFactor = double(factorial(l-m)) * (2.0*ld+1.0) / (2.0*double(factorial(l+m)));
			legenArr[l][m] *= phaseFactor * sqrtTwoPiInv * sqrt(partialNormFactor);
			phaseFactor *= -1.0;
		}
	}
	
	
	//Now, the associated Legendre Polynomials have been calculated and will be converted into a form for all lm, not just m>=0, and stored in legendre.
	legendre = new double [length];
	for (n=0; n<length; n++) {
		l = qNum[n][0];
		m = qNum[n][1];
		
		//Build the "tesseral Harmonics", but exclude the trigonometric portion; these will be referred to as the spherical Legendre Polynomials and stored in legendre.
		if (m == 0) { 
			legendre[n] = legenArr[l][m];
		}
		else if (m > 0) {
			legendre[n] = sqrt(2.0) * legenArr[l][m]; //Excluding cos(m*phi)
		}
		else {
			m *= -1; //Make m positive to effectively take the absolute value of m.
			legendre[n] = sqrt(2.0) * legenArr[l][m]; //Excluding sin(|m|*phi)
		}

	}
	
	return legendre;
}


void tesseralHarmonicsTerms(int **qNum, int length, double **legendre, double **trig, double theta, double phi) {
	//Index variables
	int m, n;
	
	double x;
	
	//Calculate the trigonometric portion of the tesseral harmonics, which are a function of m and phi only.
	*trig = new double [length];
	
	for (n=0; n<length; n++) {
		m = qNum[n][1];
		
		if (m == 0) { 
			(*trig)[n] = 1.0;
		}
		else if (m > 0) {
			(*trig)[n] = cos(double(m)*phi); //cos(m*phi)
		}
		else {
			m *= -1; //Make m positive to effectively take the absolute value of m.
			(*trig)[n] = sin(double(m)*phi); //sin(|m|*phi)
		}
	}
	
	//Calculate the spherical Legendre polynomials (i.e. the associated Legendre Polynomials times the normalization and phase factors)
	// as a function of cos(theta) for all l and m in qNum.
	x = cos(theta);
	(*legendre) = sphLegendrePoly(qNum, length, x);
	
}

int main(int argc, char** argv) {
	int l_max, length;
	int **qNum, **index, dims[2];
	int i, l, m, n;
	
	double *trig, *legendre;
	double theta, phi, sphereHarm; 
	
	double max_absError, *roots;
	int rootCount;
	
	double *roots_pn, *weights_pn;
	
	l_max = atoi(argv[1]);
	theta = atof(argv[2]);
	phi = atof(argv[3]);
	//l = atoi(argv[4]);
	//m = atoi(argv[5]);
	//max_absError = atof(argv[4]);
	
	genIndices_lm(l_max, &qNum, &length, &index, dims);
	
	tesseralHarmonicsTerms(qNum, length, &legendre, &trig, theta, phi);
	
	cout << fixed << setprecision(15);
	
	//n = index[l][m_shift(m)];
//	sphereHarm = legendre[n] * trig[n];
	
	//cout << "Spherical Harmonic: " << sphereHarm << endl;
	
	//cout << "Trig: " << endl;
//	for (n=0; n<length; n++) {
//		l = qNum[n][0];
//		m = qNum[n][1];
//		
//		cout << l << " " << m << " " << trig[n] << endl;
//	}
//	
//	cout << "Legendre: " << endl;
//	for (n=0; n<length; n++) {
//		l = qNum[n][0];
//		m = qNum[n][1];
//		
//		cout << l << " " << m << " " << legendre[n] << endl;
//	}
	
//	cout << "SphereHarm: " << endl;
//	for (n=0; n<length; n++) {
//		l = qNum[n][0];
//		m = qNum[n][1];
//		
//		cout << l << " " << m << " " << legendre[n] * trig[n] << endl;
//	} 
	
	cout << "Finding Legendre polynomial roots:" << endl;
	
	roots = legendreRoots(l_max, &rootCount);
	
	cout << "Roots Found:" << endl;
	
	for (l=0; l<rootCount; l++) {
		cout << roots[l] << endl;
	}
	
	roots_pn = new double [l_max];
	weights_pn = new double [l_max];
	
	gauleg(-1.0, 1.0,roots_pn,weights_pn,l_max);
	
	cout << "Difference between JTC and PN's roots: " << endl;
	
	for (l=0; l<l_max; l++) {
		cout << roots[l]-roots_pn[l] << endl;
	}
	
	cout << "PN's weights: " << endl;
	
	for (l=0; l<l_max; l++) {
		cout << weights_pn[l] << endl;
	}
	
	
	
	
	for (i=0 ; i<length; i++) {		
		delete [] qNum[i];
	}
	delete [] qNum;
	
	for (i=0; i<l_max; i++) {
		delete [] index[i];
	}
	delete [] index;
	delete [] legendre;
	delete [] trig;
	
	
	return 0;
}


