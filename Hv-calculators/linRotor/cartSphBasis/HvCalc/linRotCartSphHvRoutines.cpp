#include "linRotCartSphHvRoutines.h"

using namespace std;

#define m_shift(X) (X+l_max)

//Algorithm from http://www.livephysics.com/computational-physics/fortran/fortran-subroutines-functions/
double factorial(int num) {
	double fact = 1.0;
	int i;
	
	for (i=2; i<=num; i++) {
		fact *= double (i);
	}
	
	return fact;
}

//This function generates two maps for the lm basis. One, qNum, goes from n to l,m and the other, index, goes from l,m to n. 
//   The mapping order is based on method of partial summation to be used, which sums from m=-l_max to l_max as the outer loop and l = abs(m) to l_max
//  as the inner loop. i.e. first l is summed over, then m is.
// !Make sure to delete (ie. dealloc) qNum and index!
// !!Make sure to access index with the m_shift() macro!!
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
// !Make sure to delete (ie. dealloc) legendre and legendreDeriv!
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
double* legendreRoots(int l) {
	double *legendre, *legendreDeriv, *roots;
	double oldRoot, newRoot;
	double absError, absError2, max_absError, d_x;
	int i, j, /*numRoots,*/ n, rootCount, found;
	
	int iteration, max_iterations;
	
	//numRoots = l; //There are l roots for a Legendre polynomial of degree l.
	
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
	
	if (rootCount != l) {
		cout << "Wrong number of Legendre Polynomial roots found!" << endl;
		exit(1);
	}
	return roots;
}

double* legendreWeights(int l, double *roots) {
	int i;
	
	double *weights, *legendre, *legendreDeriv;
	
	weights = new double [l];
	
	for (i=0; i<l; i++) {
		legendrePoly_zerothOrder(l, &legendre, &legendreDeriv, roots[i]);
		weights[i] = 2.0 / (1.0 - (roots[i] * roots[i])) / (legendreDeriv[l]*legendreDeriv[l]);
	}
	
	return weights;
}

// Dr. P.-N. Roy's code for the calculation of the Gauss-Legendre abscissae and weights
// !Make sure to delete (ie. dealloc) x and w!
void gauleg(double x1,double x2,double *x,double *w,int n) {
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

// Dr. P.-N. Roy's code for the calculation of the Legendre polynomials
double plgndr(int l,int m,double x) {
	double fact,pll,pmm,pmmp1,somx2;
	int i,ll;
	void nrerror();
	
	/*	if (m < 0 || m > l || fabs(x) > 1.0)
	 nrerror("Bad arguments in routine PLGNDR");*/
	pmm=1.0;
	if (m > 0) {
		somx2=sqrt((1.0-x)*(1.0+x));
		fact=1.0;
		for (i=1;i<=m;i++) {
			pmm *= -fact*somx2;
			fact += 2.0;
		}
	}
	if (l == m)
		return pmm;
	else {
		pmmp1=x*(2*m+1)*pmm;
		if (l == (m+1))
			return pmmp1;
		else {
			for (ll=(m+2);ll<=l;ll++) {
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
				pmm=pmmp1;
				pmmp1=pll;
			}
			return pll;
		}
	}
}

//Calculation of the set of normalized associated Legendre polynomials for l = 0 to l_max and m = -l_max to l_max 
// This function returns a pointer to an array of the normalized associated Legendre polynomials in the [n] composite basis
double* normAssocLegendrePoly(int **qNum, int length, double x){
	//Storage arrays
	double *legendre, **legenArr;
	
	//Calculation variables
	double recFactor, phaseFactor, legenRecFactor, partialNormFactor;
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
	//    N_lm = (-1)^m * sqrt[ (l-m)! * (2l+1) / (2(l+m)!) ]
	// This returns P_lm(x).
	// Note: This code (and comment) is heavily based on the code by Toby Zeng in the file ylm_py2.f
	
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
	
	legendre = new double [length];
	
	//At this point, the Legendre Polynomials have been calculated (except for a phase factor (-1)^m).
	//Now, the acquired polynomials are going to be renormalized, with:
	//    N_lm = (-1)^m * sqrt[ (l-m)! * (2l+1) / (2(l+m)!) ]
	//#pragma omp parallel default(shared) private (l,m,ld,md, phaseFactor, partialNormFactor, n) 
	{
	//#pragma omp for schedule(guided)
	for (l=0; l<=l_max; l++) {
		phaseFactor = 1.0;
		for (m=0; m<=l; m++) {
			ld = double(l);
			md = double(m);
			
			partialNormFactor = factorial(l-m) * (2.0*ld+1.0) / (2.0* factorial(l+m));
			legenArr[l][m] *= phaseFactor * sqrt(partialNormFactor);
			phaseFactor *= -1.0;
		}
	}
	
	//Now, the associated Legendre Polynomials have been calculated and will be converted into a form for all lm, not just m>=0, and stored in legendre.	
	//#pragma omp for schedule(guided)
	for (n=0; n<length; n++) {
		l = qNum[n][0];
		m = qNum[n][1];
		
		//Calculate the Legendre Polynomials for positive and negative m
		if (m == 0) { 
			legendre[n] = legenArr[l][m];
		}
		else if (m > 0) {
			legendre[n] = legenArr[l][m];
		}
		else {
			m *= -1; //Make m positive to effectively take the absolute value of m.
			legendre[n] = legenArr[l][m];
		}
		
	}
	}
	
	for (l=0; l<=l_max; l++) {
		delete [] legenArr[l];
	}
	delete [] legenArr;
	
	return legendre;
}

double* tesseralTrigTerm(int **qNum, int length, double phi) {
	//Index variables
	int m, n;
	
	double sqrtPiInv, *trig;
	
	sqrtPiInv = 1.0 / sqrt(PI);
	//sqrtPiInv = 1.0;
	
	//Calculate the trigonometric portion of the tesseral harmonics, which are a function of m and phi only and include the factors sqrt(2) * sqrt(1/2pi) = sqrt(1/pi);
	trig = new double [length];
	
	//#pragma omp parallel for default(shared) private (n,m) schedule(guided)
	for (n=0; n<length; n++) {
		m = qNum[n][1];
		
		if (m == 0) { 
			trig[n] = sqrtPiInv / sqrt(2.0); //Also have 1/sqrt(2) as if m!=0, there is a sqrt(2) term added to cancel this 1/sqrt(2)
		}
		else if (m > 0) {
			trig[n] = sqrtPiInv * cos(double(m)*phi); //cos(m*phi)
		}
		else {
			m *= -1; //Make m positive to effectively take the absolute value of m.
			trig[n] = sqrtPiInv * sin(double(m)*phi); //sin(|m|*phi)
		}
	}
	
	return trig;
}

// !Make sure to delete (ie. dealloc) legendre and trig!
void tesseralHarmonicsTerms(int **qNum, int length, double **legendre, double **trig, double cosTheta, double phi) {
	
	//Calculate the trigonometric portion of the tesseral harmonics, which are a function of m and phi only and include the factors sqrt(2) * sqrt(1/2pi)
	(*trig) = tesseralTrigTerm(qNum, length, phi);
	
	//Calculate the normalized associated Legendre polynomials (i.e. the associated Legendre Polynomials times the normalization and phase factors)
	// as a function of cos(theta) for all l and m in qNum.
	(*legendre) = normAssocLegendrePoly(qNum, length, cosTheta);
	
}

//Gauss-Chebyshev quadrature abscissae and weights
// !Make sure to delete (ie. dealloc) abscissae and weights!
void gaussChebyshev(int numPoints, double **abscissae, double **weights){
	int i, n;
	
	n = numPoints;
	//Calculate the abscissae x_i = cos[pi(2i-1)/2n]
	(*abscissae) = new double [n];
	(*weights) = new double	[n];
	
	//#pragma omp parallel for default(shared) private (i) schedule(guided)
	for (i=0; i<n; i++) {
		(*abscissae)[i] = cos(PI * (2.0*double(i+1) - 1.0) /(2.0 * double(n)));
		(*weights)[i] = PI / double(n);
	}
}

//Gauss-Legendre quadrature abscissae and weights
// !Make sure to delete (ie. dealloc) abscissae and weights!
void gaussLegendre(int numPoints, double **abscissae, double **weights){
	int n;
	
	n = numPoints;
	
	(*abscissae) = new double [n];
	(*weights) = new double [n];
	
	(*abscissae) = legendreRoots(n);
	
	(*weights) = legendreWeights(n, (*abscissae));
}

//Cartesian Kinetic Energy operator and grid
//The grid spans [-x_max/2, x_max/2]
// !Make sure to delete (ie. dealloc) kinMat and grid!
//Changed d_x = x_max/n to d_x = x_max/(n-1) to get a proper grid - JTC 2013/08/20
void cartKinGrid(double x_max, int nPoints, double totalMass, double **kinMat, double **grid) {
	int i, j;
	double d_x;
	
	d_x = x_max / double(nPoints - 1);
	
	(*grid) = new double [nPoints];
	(*kinMat) = new double [nPoints*nPoints];
	
	
	//#pragma omp parallel default(shared) private (i,j)
	{
	//#pragma omp for schedule(guided) collapse(2) //"guided" schedule means that the iterations are dynamically dispersed to each thread dependent on the number of threads and the number of iterations remaining
	for (i=0; i<nPoints; i++) {		
		for (j=0; j<nPoints; j++) {
			(*kinMat)[i*nPoints + j] = (H_BAR*H_BAR) / (2.0 * totalMass * d_x * d_x) * pow(-1.0, (i+1)-(j+1));
			if (i==j) {
				(*kinMat)[i*nPoints + j] *= (PI*PI)/3.0;
			}
			else {
				(*kinMat)[i*nPoints + j] *= 2.0/double( ((i+1)-(j+1)) * ((i+1)-(j+1)) );
			}
			
		}
		
	}
	
	//#pragma omp for schedule(guided)
	for (i=0; i<nPoints; i++) {	//This loop separated out of the above for better parallelization
		(*grid)[i] = double(i)*d_x - x_max/2.0;
	}
	
	}
}

//Rotational Kinetic Energy Operator in the |lm> basis
// !Make sure to delete (ie. dealloc) kinVec and grid!
// Note that kinVec = diag(kinMat), since kinMat is diagonal in the |lm> basis
double* rotKinEng(int **qNum, int length, double rotationalConstant) {
	double B, *rotEng, l;
	int i;
	
	rotEng = new double [length];
	
	//B = H_BAR*H_BAR / 2.0 / momentOfInertia;
	B = rotationalConstant;
	
	//#pragma omp parallel for default(shared) private (i,l) schedule(guided)
	for (i=0; i<length; i++) {
		l = double(qNum[i][0]);
		rotEng[i] = B * l * (l+1.0);
	}
	
	return rotEng;
}

//This tests the orthonormality of the Tesseral Harmonics terms
//To do this test, add the following to "main" and comment out all else:
//int l_max, thetaPoints, phiPoints;
//
//l_max = atoi(argv[1]);
//thetaPoints = atoi(argv[2]);
//phiPoints = atoi(argv[3]);
//
//tesseralTest(l_max, thetaPoints, phiPoints);
void tesseralTest(int l_max, int thetaPoints, int phiPoints) {
	int m, l, a, b, n, i, lp, mp;
	int **qNum, length, **index, dims[2];
	
	int width;
	
	double **legendreGrid, **trigGrid, **trigGrid2;
	
	double *phiAbscissae, *phiWeights;
	
	double *cosThetaAbscissae, *cosThetaWeights;
	
	double ***Lmla, ***resMatLegendre;
	
	double ***Llma, ***resMatLegendre_m;
	
	//double const_lFactor, const_mFactor;
	
	//double normConst;
	
	double ***trigMat, ***trigResMat, ***trigMat2;
	
	//Find Gauss-Legendre and Gauss-Chebyshev grids
	gaussChebyshev(phiPoints, &phiAbscissae, &phiWeights);
	
	//gaussLegendre(thetaPoints, &cosThetaAbscissae, &cosThetaWeights);
	
	cosThetaAbscissae = new double [thetaPoints];
	cosThetaWeights = new double [thetaPoints];
	
	gauleg(-1.0, 1.0, cosThetaAbscissae, cosThetaWeights, thetaPoints);
	
	
	
	//Get the lm basis indices
	genIndices_lm(l_max, &qNum, &length, &index, dims);
	
	//Calculate the Legendre polynomials for each of the Gauss-Legendre Abscissae (cosThetaAbscissae) 
	legendreGrid = new double* [thetaPoints];
	
	//#pragma omp parallel for default(shared) private (a) schedule(guided)
	for (a=0; a<thetaPoints; a++) {
		legendreGrid[a] = normAssocLegendrePoly(qNum, length, cosThetaAbscissae[a]);
	}
	
	//	for (a=0; a<thetaPoints; a++) {
	//		legendreGrid[a] = new double [length];
	//		for (n=0; n<length; n++) {
	//			l = qNum[n][0];
	//			m = qNum[n][1];
	//			
	//			normConst = sqrt(double (2*l+1) * factorial(l-m) / 2.0 / factorial(l+m));
	//		
	//			if (m<0) {
	//				m = -m;
	//				legendreGrid[a][n] = pow(-1.0, m) * factorial(l-m) / factorial(l+m) * normConst * plgndr(l, m, cosThetaAbscissae[a]);
	//			}
	//			else {
	//				legendreGrid[a][n] = normConst * plgndr(l, m, cosThetaAbscissae[a]);
	//			}
	//
	//		}
	//	}
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//Calculate the L_la^m matrix elements (ie. constant m)
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Allocate memory for the matrix elements
	Lmla = new double** [2*l_max+1]; //Outer Loop is m
	
	for (m=-l_max; m<=l_max; m++) {
		Lmla[m_shift(m)] = new double* [l_max+1];
		
		for (l=0; l<=l_max; l++) {
			Lmla[m_shift(m)][l] = new double [thetaPoints];
		}
	}
	
	//#pragma omp parallel for default(shared) private (m,l,a) schedule(guided) collapse(3)
	for (m=-l_max; m<=l_max; m++) {		
		for (l=0; l<=l_max; l++) {			
			for (a=0; a<thetaPoints; a++) {
				Lmla[m_shift(m)][l][a] = 0.0;
			}
			
		}
	}
	
	//Calculate elements, which is sqrt(gaussLegendreWeights(a)) * ~P_lm(x_a), where x_a = cos(theta_a), ~ means normalized
	//#pragma omp parallel for default(shared) private (n,a,l,m) schedule(guided) collapse(2)
	for (n=0; n<length; n++) {
		for (a=0; a<thetaPoints; a++) {
			l = qNum[n][0];
			m = qNum[n][1];
			Lmla[m_shift(m)][l][a] = sqrt(cosThetaWeights[a]) * legendreGrid[a][index[l][m_shift(m)]];
		}
	}
	
	//Calculate L^m(L^m)^T = sum(over a){L^m_la * L^m_l'a)
	resMatLegendre = new double** [2*l_max+1];
	
	for (m=-l_max; m<=l_max; m++) {
		resMatLegendre[m_shift(m)] = new double* [l_max+1];
		
		for (l=0; l<=l_max; l++) {
			resMatLegendre[m_shift(m)][l] = new double [l_max+1];
			
			for (lp=0; lp<=l_max; lp++) {
				resMatLegendre[m_shift(m)][l][lp] = 0.0; //Initialize matrix value
			}
		}
	}
	
	//#pragma omp parallel for default(shared) private (m,l,lp,a) schedule(guided) collapse(3) //Can only be at most collapse(3) and not collapse(4), because of resMatLegendre[m_shift(m)][l][lp] (ie. it is not dependent on a, so two processors may write to the same memory location)
	for (m=-l_max; m<=l_max; m++) {		
		for (l=0; l<=l_max; l++) {			
			for (lp=0; lp<=l_max; lp++) {				
				for (a=0; a<thetaPoints; a++) {
					resMatLegendre[m_shift(m)][l][lp] += Lmla[m_shift(m)][l][a] * Lmla[m_shift(m)][lp][a];
				}
			}
		}
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//Calculate the L_ma^l matrix elements (ie. constant l)
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	// Allocate memory for the matrix elements
	Llma = new double** [l_max+1]; //Outer Loop is l
	
	for (l=0; l<=l_max; l++) {
		Llma[l] = new double* [2*l_max+1];
		
		for (m=-l_max; m<=l_max; m++) {
			Llma[l][m_shift(m)] = new double [thetaPoints];
		}
	}
	
	//#pragma omp parallel for default(shared) private (l,m,a) schedule(guided) collapse(3)
	for (l=0; l<=l_max; l++) {		
		for (m=-l_max; m<=l_max; m++) {			
			for (a=0; a<thetaPoints; a++) {
				Llma[l][m_shift(m)][a] = 0.0;
			}
		}
	}
	
	//for (a=0; a<thetaPoints; a++) {
	//		legendreGrid[a] = normAssocLegendrePoly(qNum, length, cos(cosThetaAbscissae[a]));
	//	}
	
	//Calculate elements, which is sqrt(gaussLegendreWeights(a)) * ~P_lm(x_a), where x_a = cos(theta_a), ~ means normalized
	//#pragma omp parallel for default(shared) private (n,a,l,m,mp) schedule(guided) collapse(2)
	for (n=0; n<length; n++) {
		for (a=0; a<thetaPoints; a++) { //Add sqrt(1.0-(cosThetaAbscissae[a]*cosThetaAbscissae[a])) for the othogonality condition shown at http://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Orthogonality
			l = qNum[n][0];
			m = qNum[n][1];
			
			mp = m;
			
			if (m<0) {
				mp *= -1.0;
			}
			
			Llma[l][m_shift(m)][a] = sqrt(double (2*mp) / double (2*l+1)) * sqrt(cosThetaWeights[a]) * legendreGrid[a][index[l][m_shift(m)]] / sqrt(1.0-(cosThetaAbscissae[a]*cosThetaAbscissae[a]));
			
			
		}
	}
	
	//Calculate L^l(L^l)^T = sum(over a){L_ma^l * L_m'a^l}
	resMatLegendre_m = new double** [l_max+1];
	
	for (l=0; l<=l_max; l++) {
		resMatLegendre_m[l] = new double* [2*l_max+1];
		
		for (m=-l_max; m<=l_max; m++) {
			resMatLegendre_m[l][m_shift(m)] = new double [2*l_max+1];
			
			for (mp=-l_max; mp<=l_max; mp++) {
				resMatLegendre_m[l][m_shift(m)][m_shift(mp)] = 0.0; //Initialize matrix value
			}
		}
	}
	
	//#pragma omp parallel for default(shared) private (l,m,mp,a) schedule(guided) collapse(3) //Can only be at most collapse(3) and not collapse(4), because of resMatLegendre_m[l][m_shift(m)][m_shift(mp)] (ie. it is not dependent on a, so two processors may write to the same memory location)
	for (l=0; l<=l_max; l++) {		
		for (m=-l_max; m<=l_max; m++) {			
			for (mp=-l_max; mp<=l_max; mp++) {				
				for (a=0; a<thetaPoints; a++) {
					resMatLegendre_m[l][m_shift(m)][m_shift(mp)] += Llma[l][m_shift(m)][a] * Llma[l][m_shift(mp)][a];
				}
			}
		}
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//Check the orthogonality of the trigonometric term of the tesseral harmonics
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Calculate the trigonometric term of the tesseral harmonics
	trigGrid = new double* [phiPoints];
	trigGrid2 = new double* [phiPoints];
	for (b=0; b<phiPoints; b++) { //Don't parallelize as the internal functions already are parallelized
		trigGrid[b] = tesseralTrigTerm(qNum, length, acos(phiAbscissae[b]));
		
		trigGrid2[b] = tesseralTrigTerm(qNum, length, (2.0*PI - acos(phiAbscissae[b])));
	}
	
	trigMat = new double** [l_max+1];
	trigMat2 = new double** [l_max+1];
	
	//Allocate memory
	for (l=0; l<=l_max; l++) {
		trigMat[l] = new double* [2*l_max+1];
		trigMat2[l] = new double* [2*l_max+1];
		
		for (m=-l_max; m<=l_max; m++) {
			trigMat[l][m_shift(m)] = new double [phiPoints];
			trigMat2[l][m_shift(m)] = new double [phiPoints];
			
			for (b=0; b<phiPoints; b++) {
				
				trigMat[l][m_shift(m)][b] = 0.0;
				trigMat2[l][m_shift(m)][b] = 0.0;
			}
		}
	}
	
	//Calculate the matrix elements T^l_mb = sqrt(w^GC_b) * [sqrt(2.0*PI) / sqrt(2.0)] * trig(l,m,b); [sqrt(2.0*PI) / sqrt(2.0)] is to renormalize the functions
	//#pragma omp parallel for default(shared) private (n,b,l,m) schedule(guided) collapse(2)
	for (n=0; n<length; n++) {
		for (b=0; b<phiPoints; b++) { 
			l = qNum[n][0];
			m = qNum[n][1];
			
			trigMat[l][m_shift(m)][b] = trigGrid[b][index[l][m_shift(m)]];
			trigMat2[l][m_shift(m)][b] = trigGrid2[b][index[l][m_shift(m)]];
		}
		
	}
	
	
	
	//Calculate T^l_mb(T^l_m'b)T = delta_mm'
	trigResMat = new double** [l_max+1];
	for (l=0; l<=l_max; l++) {
		trigResMat[l] = new double* [2*l_max+1];
		
		for (m=-l_max; m<=l_max; m++) {
			trigResMat[l][m_shift(m)] = new double [2*l_max+1];
			
			for (mp=-l_max; mp<=l_max; mp++) {
				trigResMat[l][m_shift(m)][m_shift(mp)] = 0.0; //Initialize values
			}
		}
	}
	
	//#pragma omp parallel for default(shared) private (l,m,mp,b) schedule(guided) collapse(3) //Can't have collapse(4) as this allows two different CPUs to access the same memory location as trigResMat does not depend on b, only l, m, and mp.
	for (l=0; l<=l_max; l++) {
		for (m=-l_max; m<=l_max; m++) {
			for (mp=-l_max; mp<=l_max; mp++) {				
				for (b=0; b<phiPoints; b++) {
					trigResMat[l][m_shift(m)][m_shift(mp)] += phiWeights[b] * ((trigMat[l][m_shift(m)][b] * trigMat[l][m_shift(mp)][b]) + (trigMat2[l][m_shift(m)][b] * trigMat2[l][m_shift(mp)][b]));
				}
			}
		}
	}
	
	//#pragma omp parallel for default(shared) private (l,m,mp,b) schedule(guided) collapse(3)
	//	for (l=0; l<(l_max+1); l++) {
	//		for (m=0; m<(2*l_max+1); m++) {
	//			for (mp=0; mp<(2*l_max+1); mp++) {				
	//				for (b=0; b<phiPoints; b++) {
	//					trigResMat[l][m][mp] += phiWeights[b] * ((trigMat[l][m][b] * trigMat[l][mp][b]) + (trigMat2[l][m][b] * trigMat2[l][mp][b]));
	//				}
	//			}
	//		}
	//	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//Print out the results
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//width = 10;
	//		cout << scientific;
	//	cout << setprecision(15);
	
	//cout << factorial(170) << endl; //170 is the max possible l value, given limitation of double format (170! = 7.257416e+306), max double = 1.797693e+308
	//cout << DBL_MAX << endl;
	
	width = 6;
	cout << fixed;
	cout << setprecision(3);
	
	//Print out Legendre Polynomial matrices
	//	cout << fixed << setprecision(3);
	//	cout << "Legendre Polynomial matrices" << endl;
	//	for (m=-l_max; m<=l_max; m++) {
	//		cout << "m = " << m << endl;
	//		
	//		for (l=0; l<=l_max; l++) {
	//			
	//			for (a=0; a<thetaPoints; a++) {
	//				cout << setw(width) << Lmla[m_shift(m)][l][a] << " ";
	//			}
	//			cout << endl;
	//		}
	//	}
	
	//Print out Legendre Polynomial result matrices for fixed m
	
	cout << "Legendre Polynomial Fixed m - Result matrices" << endl;
	for (m=-l_max; m<=l_max; m++) {
		cout << "m = " << m << endl;
		
		for (l=0; l<=l_max; l++) {
			
			for (lp=0; lp<=l_max; lp++) {
				cout << setw(width) << resMatLegendre[m_shift(m)][l][lp] << " ";
			}
			cout << endl;
		}
	}
	
	//Print out Legendre Polynomial result matrices for fixed l
	
	cout << "Legendre Polynomial Fixed l - Result matrices" << endl;
	for (l=0; l<=l_max; l++) {
		cout << "l = " << l << endl;
		
		for (m=-l_max; m<=l_max; m++) {
			
			for (mp=-l_max; mp<=l_max; mp++) {
				cout << setw(width) << resMatLegendre_m[l][m_shift(m)][m_shift(mp)] << " ";
			}
			cout << endl;
		}
	}
	
	//Print out the normalization constants
	
	//	cout << "Legendre Polynomial Fixed l - Expected Values" << endl;
	//	for (l=0; l<=l_max; l++) {
	//		cout << "l = " << l << endl;
	//		
	//		for (m=-l_max; m<=l_max; m++) {
	//			
	//			for (mp=-l_max; mp<=l_max; mp++) {
	//				if ((m==mp)&&(m!=0)) {
	//					//const_lFactor = factorial(l+m) / double (m) / factorial(l-m);
	//					//const_mFactor = factorial(l+m) * 2.0 / (double (2*l+1)) / factorial(l-m));
	//					//cout << const_lFactor/const_mFactor << " ";
	//					
	//					cout << setw(width) << double (2*l+1) / double (2*m) << " ";
	//				}
	//				else if ((m==mp)&&(m==0)) {
	//					cout << setw(width) << "Inf" << " ";
	//				}
	//				else if ((m!=mp)) {
	//					cout << setw(width) << 0.0 << " ";
	//				}
	//				else {
	//					cerr << "Error!" << endl;
	//					exit(1);
	//				}
	//
	//			}
	//			cout << endl;
	//		}
	//	}
	
	//Print out Trig Term result matrices for fixed l
	
	cout << "Trig Term Fixed l - Result matrices" << endl;
	for (l=0; l<=l_max; l++) {
		cout << "l = " << l << endl;
		
		for (m=-l_max; m<=l_max; m++) {
			
			for (mp=-l_max; mp<=l_max; mp++) {
				cout << setw(width) << trigResMat[l][m_shift(m)][m_shift(mp)] << " ";
			}
			cout << endl;
		}
	}
	
	//Deallocate Quadrature Arrays
	delete [] phiAbscissae;
	delete [] phiWeights;
	
	delete [] cosThetaAbscissae;
	delete [] cosThetaWeights;
	
	//Deallocate legendre polynomials
	for (a=0; a<thetaPoints; a++) {
		delete [] legendreGrid[a];
	}
	
	delete [] legendreGrid;
	
	//Deallocate Legendre Polynomial matrices
	for (m=-l_max; m<=l_max; m++) {
		
		for (l=0; l<=l_max; l++) {
			delete [] Lmla[m_shift(m)][l];
		}
		delete [] Lmla[m_shift(m)];
	}
	delete [] Lmla;
	
	//Deallocate Legendre Polynomial result matrices
	for (m=-l_max; m<=l_max; m++) {
		
		for (l=0; l<=l_max; l++) {
			delete [] resMatLegendre[m_shift(m)][l];
		}
		delete [] resMatLegendre[m_shift(m)];
	}
	delete [] resMatLegendre;
	
	//Deallocate |lm> basis maps	
	for (i=0; i<length; i++) {
		delete [] qNum[i];
	}
	delete [] qNum;
	
	for (i=0; i<dims[0]; i++) {
		delete [] index[i];
	}
	delete	[] index;
	
	
}

//This calculates the vector multiplied by the potential, for a given xyz and either acos(cosPhiAbscissae) or 2PI - acos(cosPhiAbscissae)
//Test this with both a potential = 1 and potential = cos(theta) or cos(phi)
double* calc_ulm(double x, double y, double z, double *v_lpmp, interfaceStor *interface, int rangeFlag) {
	int m, mp, n, a, b;
	int anmp, nmp, nm, bna, anb, mna;
	double potentialCeiling;
	
	//Extract out desired variables from the interface storage structure
	quadStor *gaussQuad;
	pointPotentialStorH2 *atomPotentials;
	tesseralStor *tesseralHarmonics;
	lmFBR *lmBasis;
	
	int l_max, **qNum, length;
	double (*linearMoleculePotential)(interfaceStor*, H2_orient*) = NULL;
	
	gaussQuad = interface->quadrature;
	atomPotentials = interface->potential;
	lmBasis = interface->lmBasis;
	
	
	l_max = lmBasis->lmax;
	qNum = lmBasis->qNum;
	length = lmBasis->length;
	
	potentialCeiling = atomPotentials->potentialCeiling;
	linearMoleculePotential = interface->fcnPointers->linearMoleculePotential;
	
	
	//Determine whether I am calculating phi = acos(cosPhiAbscissae) or phi = 2PI - acos(cosPhiAbscissae)
	if (rangeFlag == 0) {
		tesseralHarmonics = interface->tesseral;
	}
	else {
		tesseralHarmonics = interface->tesseral2PI;
	}
	
	//Get number of m values
	nmp = 2*l_max + 1;
	nm = 2*l_max + 1;
	
	//Quadrature Variables
	double *cosThetaAbscissae, *cosPhiAbscissae, *wa, *wb;
	int na, nb;
	
	cosThetaAbscissae = gaussQuad->GLabscissae;
	wa = gaussQuad->GLweights;
	na = gaussQuad->GLnum;
	
	cosPhiAbscissae = gaussQuad->GCabscissae;
	wb = gaussQuad->GCweights;
	nb = gaussQuad->GCnum;
	//
	
	//Tesseral Harmonics Variables
	double **L_lpmp, **S_mp, **S_m, **L_lm;
	
	L_lpmp = tesseralHarmonics->L_lpmp;
	S_mp = tesseralHarmonics->S_mp;
	S_m = tesseralHarmonics->S_m;
	L_lm = tesseralHarmonics->L_lm;
	//
	
	//Memory allocation and assignment
	//Loop 1
	double *vt_mpa;
	vt_mpa = new double [na*nmp];
	
	//Loop 2
	double	*u_ab;
	u_ab = new double [na*nb];
	
	//Loop 3
	double *ut_ab;
	ut_ab = new double [na*nb];
	VECT CMpos(3);
	H2_orient linearMolecule;	
	double V_ab;
	
	//Loop 4
	double *ut_ma;
	ut_ma = new double [nm*na];
	
	//Loop 5
	double *u_lm;
	u_lm = new double [length];
	
//#pragma omp parallel default(shared) private (a,mp,anmp,n,b,bna,m,mna,anb,linearMolecule,CMpos,V_ab)
	{
	
	//Loop 1 - vt_mpa = L_lpmp(q_a) * v_lpmp; FLOPS = na * length = na * (l_max+1)^2
	//#pragma omp for schedule(guided) collapse(2)
	for (a = 0; a<na; a++) {
		for (mp=-l_max; mp<=l_max; mp++) {
			anmp = a * nmp;			
			vt_mpa[anmp + m_shift(mp)] = 0.0; //Initialize vector
			
			for (n=0; n<length; n++) {
				vt_mpa[anmp + m_shift(mp)] += L_lpmp[a][n] * v_lpmp[n];
			}
		}
	}
	
	//Loop 2 - u_ab = S_mp(Pb) * vt_mpa; FLOPS = na * nb * nm = na * nb * (2l_max+1)
	//#pragma omp for schedule(guided) collapse(2)
	for (b=0; b<nb; b++) {
		for  (a=0; a<na; a++){
			bna = b * na;
			anmp = a * nmp;			
			u_ab[bna + a] = 0.0;
			
			for (mp=-l_max; mp<=l_max; mp++) {
				u_ab[bna + a] += S_mp[b][m_shift(mp)] * vt_mpa[anmp + m_shift(mp)]; //Pb------------------------------
			}
		}
	}
	
	
	//Loop 3 - ut_ab = V_ab(xyz, Qa, Pb) * u_ab; FLOPS = na * nb
	CMpos.DIM(3);
	CMpos.COOR(0) = x;
	CMpos.COOR(1) = y;
	CMpos.COOR(2) = z;
	
	linearMolecule.CM = &CMpos;
	
//#pragma omp for schedule(guided) collapse(2)
	for (b=0; b<nb; b++) {
		for (a=0; a<na; a++) {
			bna = b * na;
			anb = a * nb;
			
			linearMolecule.theta = acos(cosThetaAbscissae[a]);
			//Determine whether I am calculating phi = acos(cosPhiAbscissae) or phi = 2PI - acos(cosPhiAbscissae) //Pb------------------------------
			if (rangeFlag == 0) { //This if statement may be SLOW!!!!!!!!!!!!!
				linearMolecule.phi = acos(cosPhiAbscissae[b]); 
			}
			else {
				linearMolecule.phi = 2*PI - acos(cosPhiAbscissae[b]);
			}
			
			//Calculate potential at x, y, z, theta, phi EDIT
			V_ab = (*linearMoleculePotential)(interface, &linearMolecule);
			
			//cout << "VvPot: " << V_ab << " " << potentialCeiling << endl;
			
			//V_ab = 0.0; //Set to zero for debugging purposes
			
			if (V_ab >= potentialCeiling) {
				V_ab = potentialCeiling;
			}
			
			ut_ab[anb + b] = V_ab * u_ab[bna + a];
		}
	}
	
	//Loop 4 - ut_ma = wb * S_m(Pb) * ut_ab; FLOPS = na * nb * nm = na * nb * (2l_max+1)
	//#pragma omp for schedule(guided) collapse(2)
	for (m=-l_max; m<=l_max; m++) {
		for (a=0; a<na; a++) {
			mna = m_shift(m) * na;
			ut_ma[mna + a] = 0.0;
			anb = a * nb;
			
			for (b=0; b<nb; b++) {
				ut_ma[mna + a] += wb[b] * S_m[m_shift(m)][b] * ut_ab[anb + b]; //Pb------------------------------
			}
		}
	}
	
	//Loop 5 - u_lm = wa * L_lm(Qa) * ut_ma; FLOPS = na * length = na * (l_max+1)^2
	//#pragma omp for schedule(guided)
	for (n=0; n<length; n++) {
		u_lm[n] = 0.0;
		
		m = qNum[n][1];
		mna = m_shift(m) * na;
		
		for (a=0; a<na; a++) {
			u_lm[n] = wa[a] * L_lm[n][a] * ut_ma[mna + a];
		}
	}
	}
	
	delete [] ut_ma;
	delete [] u_ab;
	delete [] ut_ab;
	delete [] vt_mpa;
	
	return u_lm;
	
}

void HvPrep_Internal(int argc, char **argv, interfaceStor *interface, lanczosStor *lanczos) {
	int na, nb;
	
	cout << endl;
	cout << "Hv Preparation STARTED." << endl;
	cout << "////////////////////////////////////////////////////////////////////////////" << endl;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Get required information from input file 
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	string inputFilename;
	ifstream inputFile;
	
	int nx, ny, nz, l_max, thetaPoints, phiPoints, pnx, pny, pnz;
	double x_max, y_max, z_max, px_max, py_max, pz_max;
	double rotationalConstant, totalMass, ceilingPotential;
	string geometryFilename, line, junk, simulationFilename, quadConvergeStudy;
	
	inputFilename = argv[2];
	
	inputFile.open(inputFilename.c_str(), ios::in);
	if (inputFile.is_open()) {
		//Get rid of comment lines;
		getline(inputFile, line);
		getline(inputFile, line);
		
		//Gather data for input; each line has the name of the value separate from the value with a space
		//Always ignore the name of the value.
		inputFile >> junk;
		inputFile >> nx;		
		
		inputFile >> junk;
		inputFile >> x_max;
		
		inputFile >> junk;
		inputFile >> ny;
		
		inputFile >> junk;
		inputFile >> y_max;
		
		inputFile >> junk;
		inputFile >> nz;
		
		inputFile >> junk;
		inputFile >> z_max;
		
		inputFile >> junk;
		inputFile >> l_max;
		
		inputFile >> junk;
		inputFile >> thetaPoints;
		
		inputFile >> junk;
		inputFile >> phiPoints;
		
		inputFile >> junk;
		inputFile >> pnx;
		
		inputFile >> junk;
		inputFile >> px_max;
		
		inputFile >> junk;
		inputFile >> pny;
		
		inputFile >> junk;
		inputFile >> py_max;
		
		inputFile >> junk;
		inputFile >> pnz;
		
		inputFile >> junk;
		inputFile >> pz_max;
		
		inputFile >> junk;
		inputFile >> totalMass;
		
		inputFile >> junk;
		inputFile >> rotationalConstant;
		
		inputFile >> junk;
		inputFile >> geometryFilename;
		
		inputFile >> junk;
		inputFile >> simulationFilename;
		
		inputFile >> junk;
		inputFile >> ceilingPotential;
		
		inputFile >> junk;
		inputFile >> quadConvergeStudy;
		
	}
	else {
		cerr << "Input file '" << inputFilename << "' could not be opened." << endl;
		exit(1);
	}
	
	inputFile.close();
	
	//Check to ensure the geometry file can be opened before doing the rest of the calculations
	//NULL is to allow for cases where the geometry file is not needed
	inputFile.open(geometryFilename.c_str(), ios::in);
	if (!(inputFile.is_open()) && (geometryFilename != "NULL")) {
		cerr << "Atom geometry file '" << geometryFilename << "' could not be opened." << endl;
		exit(1);
	}
	inputFile.close();
	
	//Check that the system dimensions are reasonable
	if (x_max > px_max || y_max > py_max || z_max > pz_max) {
		cerr << "Error, the system's physical dimensions must be the same size or smaller than the potential universe dimensions." << endl;
		exit(1);
	}
	
	if (x_max*px_max*y_max*py_max*z_max*pz_max < DBL_EPSILON || nx*ny*nz*pnx*pny*pnz == 0) {
		cerr << "Error, none of the system dimensions nor the number of points can be zero and none should be smaller than the machine epsilon: " << DBL_EPSILON << endl;
		exit(1);
	}
	
	cout << "Hv Input file parameters read." << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Generate the x, y, and z bases and Kinetic Energy Operators
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double *xKinmat, *xGrid;
	
	cartKinGrid(x_max, nx, totalMass, &xKinmat, &xGrid); //grid[i], kinMat[(i*nx) + j]
	
	double *yKinmat, *yGrid;
	
	if ((fabs(y_max-x_max)<DBL_EPSILON) && (ny==nx)) { //If the y-grid is the same as the x-grid, reuse the matrix and save memory and time
		yKinmat = xKinmat;
		yGrid = xGrid;
	}
	else {
		cartKinGrid(y_max, ny, totalMass, &yKinmat, &yGrid); //grid[i], kinMat[(i*ny) + j]
	}
	
	double *zKinmat, *zGrid;
	
	if ((fabs(z_max-x_max)<DBL_EPSILON) && (nz==nx)) { //If the z-grid is the same as the x-grid or the y-grid, reuse the matrix and save memory and time
		zKinmat = xKinmat;
		zGrid = xGrid;
	}
	else if ((fabs(z_max-y_max)<DBL_EPSILON) && (nz==ny)) {
		zKinmat = yKinmat;
		zGrid = yGrid;
	}
	else {
		cartKinGrid(z_max, nz, totalMass, &zKinmat, &zGrid); //grid[i], kinMat[(i*nz) + j]
	}
	
	//Store everything for interface
	fiveDGrid *gridStor = new fiveDGrid();
	
	gridStor->x_Grid = xGrid;
	gridStor->nx = nx;
	gridStor->x_max = x_max;
	gridStor->xKinMat = xKinmat;
	
	gridStor->y_Grid = yGrid;
	gridStor->ny = ny;
	gridStor->y_max = y_max;
	gridStor->yKinMat = yKinmat;
	
	gridStor->z_Grid = zGrid;
	gridStor->nz = nz;
	gridStor->z_max = z_max;
	gridStor->zKinMat = zKinmat;	
	
	cout << "Cartesian grid and Kinetic Energy Operators generated." << endl;
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Generate the |lm> basis and Rotational Kinetic Energy Operator
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Generate the |lm> basis
	int **qNum, length, **index, *dims;
	dims = new int [2];
	genIndices_lm(l_max, &qNum, &length, &index, dims);
	
	//Calculate the rotational Kinetic Energy operator
	double *rotKinEngOp;
	
	rotKinEngOp = rotKinEng(qNum, length, rotationalConstant); //This vector is based on the [n] composite basis
	
	//Store everything for interface
	lmFBR *lmBasis = new lmFBR();
	
	lmBasis->lmax = l_max;
	
	lmBasis->qNum = qNum;
	lmBasis->length = length;
	
	lmBasis->index = index;
	lmBasis->dims = dims;
	
	lmBasis->rotKinMat = rotKinEngOp;
	
	cout << "|lm> basis and Rotational Kinetic Energy Operator generated." << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Pre-compute everything required for the potential
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Calculate abscissae and weights for the quadratures
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	quadStor *quadrature;
	
	quadrature = QuadraturePrep(thetaPoints, phiPoints);
	
	cout << "Gauss-Legendre and Gauss-Chebyshev quadrature weights and abscissae calculated." << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Calculate Tesseral Harmonics
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	tesseralStor *tessHarmonics, *tessHarmonics2PI;
	
	na = thetaPoints;
	nb = phiPoints;
	
	TesseralPrep(na, nb, quadrature, lmBasis, &tessHarmonics, &tessHarmonics2PI);
		
	cout << "Tesseral Harmonics terms calculated." << endl;
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Pre-Calculate the Potential
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	cout << "Pre-calculating the potential:" << endl;
	cout << "----------------------------------------------------------------------------" << endl;
	
	int numDim = 3; //There are 3 spatial dimensions
	double *gridMax = new double [3];
	int *gridPoints = new int [3];
	
	gridMax[0] = px_max;
	gridPoints[0] = pnx;
	
	gridMax[1] = py_max;
	gridPoints[1] = pny;
	
	gridMax[2] = pz_max;
	gridPoints[2] = pnz;
	
	
	pointPotentialStorH2 *partialPotential;
	
	partialPotential = (*(interface->fcnPointers->preCalcPotential))(numDim, gridMax, gridPoints, geometryFilename, interface, argc, argv);
	
	partialPotential->potentialCeiling = ceilingPotential;
	
	cout << "Potentials Pre-calculated." << endl;
	
	
	
	
	//Put everything in the interface sructure
	interface->grids = gridStor;
	interface->lmBasis = lmBasis;
	interface->quadrature = quadrature;
	interface->tesseral = tessHarmonics;
	interface->tesseral2PI = tessHarmonics2PI;
	interface->potential = partialPotential;
	
	//Put what is needed for the Lanczos algorithm in lanczos	
	lanczos->total_basis_size = nx*ny*nz*length;
	lanczos->sim_descr = simulationFilename; //Keep it at this for now.
	lanczos->sim_descr_short = simulationFilename;
	
	if (quadConvergeStudy == "TRUE") {
		quadratureConvergenceStudy(interface, lanczos);
	}
	else {
		cout << "No convergence study performed, using values in : " << inputFilename << endl;
	}

	cout << "////////////////////////////////////////////////////////////////////////////" << endl;
	
	cout << "Hv Preparation FINISHED." << endl;
}

quadStor* QuadraturePrep(int thetaPoints, int phiPoints) {
	//Gauss-Legendre
	double minVal, maxVal, *cosThetaAbscissae, *cosThetaWeights;
	minVal = -1.0;
	maxVal = 1.0;
	
	cosThetaAbscissae = new double [thetaPoints];
	cosThetaWeights = new double [thetaPoints];
	//na = thetaPoints;
	//nb = phiPoints;
	
	gauleg(minVal, maxVal, cosThetaAbscissae, cosThetaWeights, thetaPoints);
	
	//Gauss-Chebyshev
	double *cosPhiAbscissae, *cosPhiWeights;
	
	gaussChebyshev(phiPoints, &cosPhiAbscissae, &cosPhiWeights);
	
	//Pre-compute values using acos for the sake of speed
	int a,b;
	//Theta [0,PI]
	double *thetaAbscissae = new double [thetaPoints];
	
	for (a; a<thetaPoints; a++) {
		thetaAbscissae[a] = acos(cosThetaAbscissae[a]);
	}
	
	//Phi [0,PI) and [PI,2PI)
	double *phiAbscissae = new double [phiPoints];
	double *PIphiAbscissae = new double [phiPoints];
	
	for (b; b<phiPoints; b++) {
		phiAbscissae[b] = acos(cosPhiAbscissae[b]);
		PIphiAbscissae[b] = 2*PI - acos(cosPhiAbscissae[b]);
	}
	
	
	//Store everything for interface
	quadStor *quadrature = new quadStor();
	
	quadrature->GCabscissae = cosPhiAbscissae;
	quadrature->GCacosAbscissae = phiAbscissae;
	quadrature->GCPIacosAbscissae = PIphiAbscissae;
	quadrature->GCweights = cosPhiWeights;
	quadrature->GCnum = phiPoints;
	
	quadrature->GLabscissae = cosThetaAbscissae;
	quadrature->GLacosAbscissae = thetaAbscissae;
	quadrature->GLweights = cosThetaWeights;
	quadrature->GLnum = thetaPoints;	
	
	return quadrature;
}

void TesseralPrep(int na, int nb, quadStor *quadrature, lmFBR *lmBasis, tesseralStor **tessHarmonicsStor, tesseralStor **tessHarmonics2PIStor) {
	
	//Calculate the Tesseral Harmonics terms and rearrange appropriately for phi = acos(cosPhiAbscissae)
	tesseralStor *tessHarmonics = new tesseralStor();
	double *stor, *cosPhiAbscissae, *cosThetaAbscissae;
	int n, a, b, m, l, nm;
	
	int l_max, length, **qNum, **index;
	
	qNum = lmBasis->qNum;
	length = lmBasis->length;
	l_max = lmBasis->lmax;
	index = lmBasis->index;
	nm = 2*l_max + 1;
	
	cosPhiAbscissae = quadrature->GCabscissae;
	cosThetaAbscissae = quadrature->GLabscissae;
	
	tessHarmonics->na = na;
	tessHarmonics->nb = nb;
	tessHarmonics->lmax = l_max;
	
	//Memory allocation
	
	//Sec 1	
	tessHarmonics->L_lpmp = new double* [na];
	
	//Sec 2
	tessHarmonics->S_mp = new double* [nb];
	
	for (b=0; b<nb; b++) {
		tessHarmonics->S_mp[b] = new double [nm];
	}
	
	//Sec3
	tessHarmonics->L_lm = new double* [length];
	
	for (n=0; n<length; n++) {
		tessHarmonics->L_lm[n] = new double [na];
	}
	
	//Sec 4
	tessHarmonics->S_m = new double* [nm];
	
	for (m=-l_max; m<=l_max; m++) {
		tessHarmonics->S_m[m_shift(m)] = new double [nb];
	}
	
	//Sec 5
	tesseralStor *tessHarmonics2PI = new tesseralStor();
	
	tessHarmonics2PI->na = na;
	tessHarmonics2PI->nb = nb;
	tessHarmonics2PI->lmax = l_max;
	
	tessHarmonics2PI->L_lpmp = new double* [na];
	
	//Sec 6
	tessHarmonics2PI->S_mp = new double* [nb];
	
	for (b=0; b<nb; b++) {
		tessHarmonics2PI->S_mp[b] = new double [nm];
	}
	
	//Sec 7
	tessHarmonics2PI->L_lm = new double* [length];
	
	for (n=0; n<length; n++) {
		tessHarmonics2PI->L_lm[n] = new double [na];
	}
	
	//Sec 8
	tessHarmonics2PI->S_m = new double* [nm];
	
	for (m=-l_max; m<=l_max; m++) {
		tessHarmonics2PI->S_m[m_shift(m)] = new double [nb];
	}
	
	
	//#pragma omp parallel default(shared) private(a,b,m,l,n,stor)
	{
	
	//Sec 1	
	//#pragma omp for schedule(guided)
	for (a=0; a<na; a++) {
		tessHarmonics->L_lpmp[a] = normAssocLegendrePoly(qNum, length, cosThetaAbscissae[a]);
	}
	
	//Sec 2
	//#pragma omp for schedule(guided)
	for (b=0; b<nb; b++) {		
		stor = tesseralTrigTerm(qNum, length, acos(cosPhiAbscissae[b])); //Phi --------------------------------------------------------
		
		for (m=-l_max; m<=l_max; m++) { //Reshift indices around
			l=l_max;
			n = index[l][m_shift(m)];
			
			if (n<0) {
				cerr << "Error in finding the index within HvPrep_Internal for S_mp" << endl;
				exit(1);
			}
			
			tessHarmonics->S_mp[b][m_shift(m)] = stor[n];			
		}
		
		delete [] stor;
	}
	
	//Sec3
	//#pragma omp for schedule(guided)
	for (a=0; a<na; a++) {
		stor = normAssocLegendrePoly(qNum, length, cosThetaAbscissae[a]);
		
		for (n=0; n<length; n++) {
			tessHarmonics->L_lm[n][a] = stor[n];
		}
		delete [] stor;
	}
	
	//Sec 4	
	//#pragma omp for schedule(guided)
	for (b=0; b<nb; b++) {
		stor = tesseralTrigTerm(qNum, length, acos(cosPhiAbscissae[b])); //Phi --------------------------------------------------------
		
		for (m=-l_max; m<=l_max; m++) {
			l=l_max;
			n = index[l][m_shift(m)];
			
			if (n<0) {
				cerr << "Error in finding the index within HvPrep_Internal for S_m" << endl;
				exit(1);
			}
			
			tessHarmonics->S_m[m_shift(m)][b] = stor[n];
		}
		delete [] stor;
	}
	
	//Calculate the Tesseral Harmonics terms and rearrange appropriately for phi = 2PI - acos(cosPhiAbscissae)
	
	//Sec 5	
	//#pragma omp for schedule(guided)
	for (a=0; a<na; a++) {
		tessHarmonics2PI->L_lpmp[a] = normAssocLegendrePoly(qNum, length, cosThetaAbscissae[a]);
	}
	
	//Sec 6	
	//#pragma omp for schedule(guided)
	for (b=0; b<nb; b++) {		
		stor = tesseralTrigTerm(qNum, length, 2.0*PI - acos(cosPhiAbscissae[b])); //Phi --------------------------------------------------------
		
		for (m=-l_max; m<=l_max; m++) { //Reshift indices around
			l=l_max;
			n = index[l][m_shift(m)];
			
			if (n<0) {
				cerr << "Error in finding the index within HvPrep_Internal for S_mp" << endl;
				exit(1);
			}
			
			tessHarmonics2PI->S_mp[b][m_shift(m)] = stor[n];			
		}
		delete [] stor;
	}
	
	//Sec 7
	//#pragma omp for schedule(guided)
	for (a=0; a<na; a++) {
		stor = normAssocLegendrePoly(qNum, length, cosThetaAbscissae[a]);
		
		for (n=0; n<length; n++) {
			tessHarmonics2PI->L_lm[n][a] = stor[n];
		}
		delete [] stor;
	}
	
	//Sec 8
	//#pragma omp for schedule(guided)
	for (b=0; b<nb; b++) {
		stor = tesseralTrigTerm(qNum, length, 2.0*PI - acos(cosPhiAbscissae[b])); //Phi --------------------------------------------------------
		
		for (m=-l_max; m<=l_max; m++) {
			l=l_max;
			n = index[l][m_shift(m)];
			
			if (n<0) {
				cerr << "Error in finding the index within HvPrep_Internal for S_m" << endl;
				exit(1);
			}
			
			tessHarmonics2PI->S_m[m_shift(m)][b] = stor[n];
		}
		delete [] stor;
		
	}
	}
	
	//Return storage containers	
	*tessHarmonicsStor = tessHarmonics;
	*tessHarmonics2PIStor = tessHarmonics2PI;
	
}				

void quadratureConvergenceStudy(interfaceStor *interface, lanczosStor *lanczos) {

	int i,a,b,j;
	int basisSize;
	int na, nb;
	
	basisSize = lanczos->total_basis_size;
	
	na = interface->quadrature->GLnum;
	nb = interface->quadrature->GCnum;
	
	double* vec=new double[basisSize];
	double* uec;
	double expectVal = 0.0;
	
	quadStor *quadrature;
	tesseralStor *tessHarmonics, *tessHarmonics2PI;
	
	string outputFilename = lanczos->sim_descr_short;
	
	ofstream outputFile(outputFilename.c_str());

	cout << "Quadrature Convergence Study beginning." << endl;
	
	for (i=0; i<basisSize; i++) {
		vec[i] = 1.0/sqrt((double)basisSize);
	}
	
	outputFile << "#ThetaPoints" << " " << "PhiPoints" << " " << "<v|H|v>" << endl;
	
	j = 0;
	
	for (a=1; a<=na; a++) {
		for (b=1; b<=nb; b++) {
			j++;
			
			quadrature = QuadraturePrep(a, b);
			
			TesseralPrep(a, b, quadrature, interface->lmBasis, &tessHarmonics, &tessHarmonics2PI);
			
			delete interface->quadrature;
			delete interface->tesseral;
			delete interface->tesseral2PI;
			
			interface->quadrature = quadrature;
			interface->tesseral = tessHarmonics;
			interface->tesseral2PI = tessHarmonics2PI;
			
			uec = Hv_5D_oneCompositeIndex(interface, vec);
			
			for (i=0; i<basisSize; i++) {
				expectVal += vec[i]*uec[i];
			}
			
			outputFile << a << " " << b << " " << expectVal << endl;
			//outputFile << a << " " << expectVal << endl;
			
			if (j % 10 == 0) {
				cout << "Working... at " << j << " out of " << na*nb << endl;
			}
			
			expectVal = 0.0;
			
		}
	}
	
	
	
}

double* Mv_5D_oneCompositeIndex(double *v_ipjkn, double *mat_iip, int ni, int nj, int nk, int nn) { 
	int i, ip, j, k, n;
	
	double *v_ijkn = new double [ni*nj*nk*nn];
	
//#pragma omp parallel for default(shared) private(n,k,j,i,ip) schedule(guided) collapse(4)
	for (n=0; n<nn; n++) {
		for (k=0; k<nk; k++) {			
			for (j=0; j<nj; j++) {			
				for (i=0; i<ni; i++) {
					
					v_ijkn[((n*nk + k)*nj + j)*ni + i] = 0.0;
					
					for (ip=0; ip<ni; ip++) {
						v_ijkn[((n*nk + k)*nj + j)*ni + i] += mat_iip[i*ni + ip] * v_ipjkn[((n*nk + k)*nj + j)*ni + ip];
					}
				}
			}
		}
	}
	
	return v_ijkn;
}

double* diagMv_5D_oneCompositeIndex(double *v_npijk, double *mat_n, int nn, int ni, int nj, int nk) { 
	int i, j, k, n;
	
	double *v_nijk = new double [ni*nj*nk*nn];
	
	//#pragma omp parallel for default(shared) private(k,j,i,n) schedule(guided) collapse(4)
	for (k=0; k<nk; k++) {
		for (j=0; j<nj; j++) {
			for (i=0; i<ni; i++) {
				for (n=0; n<nn; n++) {
					v_nijk[((k*nj + j)*ni + i)*nn + n] = mat_n[n] * v_npijk[((k*nj + j)*ni + i)*nn + n];
				}
			}
		}
	}
	
	return v_nijk;
}

//This function rotates the indices of a vector (1D array with multiple, nested quantum numbers) left one; ie. v_ijkn -> v_jkni
double* reshuffleIndices_5D_oneCompositeIndex(double *v_ijkn, int ni, int nj, int nk, int nn) {
	int i, j, k, n;
	
	double *v_jkni = new double [ni*nj*nk*nn];
	
//#pragma omp parallel for default(shared) private(n,k,j,i) schedule(guided) collapse(4)
	for (n=0; n<nn; n++) {
		for (k=0; k<nk; k++) {			
			for (j=0; j<nj; j++) {				
				for (i=0; i<ni; i++) {
					//The below equivalent to v_jkni[((i*nn + n)*nk + k)*nj + j] = v_ijkn[((n*nk + k)*nj + j)*ni + i];
					v_jkni[((i*nn + n)*nk + k)*nj + j] = v_ijkn[((n*nk + k)*nj + j)*ni + i]; 
				}
			}
		}
	}
	
	return v_jkni;
}

//This function calculates the Kinetic Energy terms of the Hamiltonian multiplied by the vector
double* Tv_5D_oneCompositeIndex(interfaceStor *interface, double *v_ipjkn) {
	int p;
	int ni, nj, nk, nn;
	
	double *Tx, *Ty, *Tz, *Trot;
	
	ni = interface->grids->nx;
	nj = interface->grids->ny;
	nk = interface->grids->nz;
	
	nn = interface->lmBasis->length;
	
	int basis_size = ni*nj*nk*nn;
	
	//////////////////////////////////////
	//Calculate T_x * v_ipjkn
	//////////////////////////////////////
	double *v_ijkn, *v_jkni_TxTerm, *v_jpkni;
	
	Tx = interface->grids->xKinMat;
	v_ijkn = Mv_5D_oneCompositeIndex(v_ipjkn, Tx, ni, nj, nk, nn);
	
	//Reshuffle indices
	v_jkni_TxTerm = reshuffleIndices_5D_oneCompositeIndex(v_ijkn, ni, nj, nk, nn);
	v_jpkni = reshuffleIndices_5D_oneCompositeIndex(v_ipjkn, ni, nj, nk, nn);
	
	delete [] v_ijkn; //Remove redundant memory (NOTE: can't delete v_ipjkn as it is needed elsewhere in Hv)
	
	
	//////////////////////////////////////
	//Calculate T_y * v_jpkni
	//////////////////////////////////////
	double *v_jkni, *v_knij_TyTerm, *v_kpnij;
	
	Ty = interface->grids->yKinMat;
	v_jkni = Mv_5D_oneCompositeIndex(v_jpkni, Ty, nj, nk, nn, ni);
	
	//Perform the sum Tx*v + Ty*v = v_jkni_TxTerm + v_jkni
	//#pragma omp parallel for default(shared) private(p) schedule(guided) //Parallelization justified as no processor will access the same memory location at the same time (i.e. p is different for each processor)
	for (p=0; p<basis_size; p++) {
		v_jkni[p] += v_jkni_TxTerm[p];
	}
	delete [] v_jkni_TxTerm; //Don't need this anymore
	
	//Reshuffle indices
	v_knij_TyTerm = reshuffleIndices_5D_oneCompositeIndex(v_jkni, nj, nk, nn, ni);
	v_kpnij = reshuffleIndices_5D_oneCompositeIndex(v_jpkni, nj, nk, nn, ni);
	
	delete [] v_jkni;
	delete [] v_jpkni; //Only need v_ipjkn, not its reshuffled forms elsewhere in Hv
	
	
	//////////////////////////////////////
	//Calculate T_z * v_kpnij
	//////////////////////////////////////
	double *v_knij, *v_nijk_TzTerm, *v_npijk;
	
	Tz = interface->grids->zKinMat;
	v_knij = Mv_5D_oneCompositeIndex(v_kpnij, Tz, nk, nn, ni, nj);
	
	//Perform the sum (Tx*v + Ty*v) + Tz*v = v_knij_TyTerm + v_knij
	//#pragma omp parallel for default(shared) private(p) schedule(guided) //Parallelization justified as no processor will access the same memory location at the same time (i.e. p is different for each processor)
	for (p=0; p<basis_size; p++) {
		v_knij[p] += v_knij_TyTerm[p];
	}
	delete [] v_knij_TyTerm;
	
	//Reshuffle Indices
	v_nijk_TzTerm = reshuffleIndices_5D_oneCompositeIndex(v_knij, nk, nn, ni, nj);
	v_npijk = reshuffleIndices_5D_oneCompositeIndex(v_kpnij, nk, nn, ni, nj);
	
	delete [] v_knij;
	delete [] v_kpnij;
	
	//////////////////////////////////////
	//Calculate T_rot * v_npijk
	//////////////////////////////////////
	double *v_nijk;
	
	Trot = interface->lmBasis->rotKinMat;
	v_nijk = diagMv_5D_oneCompositeIndex(v_npijk, Trot, nn, ni, nj, nk);
	
	//Perform the sum (Tx*v + Ty*v + Tz*v) + Trot*v = v_nijk_TzTerm + v_nijk
	//#pragma omp parallel for default(shared) private(p) schedule(guided) //Parallelization justified as no processor will access the same memory location at the same time (i.e. p is different for each processor)
	for (p=0; p<basis_size; p++) {
		v_nijk[p] += v_nijk_TzTerm[p];
		//v_nijk[p] += 0.0; //Only keep rotational term - for debugging purposes
	}
	delete [] v_nijk_TzTerm;
	
	//Reshuffle Indices back to original state
	v_ijkn = reshuffleIndices_5D_oneCompositeIndex(v_nijk, nn, ni, nj, nk);
	
	delete [] v_npijk;
	delete [] v_nijk;
	
	return v_ijkn;
}

//This function calculates the potential multiplied by the vector
double* Vv_5D_oneCompositeIndex(interfaceStor *interface, double *v_ipjkn) {
	int i, j, k, n;
	
	int ni, nj, nk, nn, basis_size;	
	ni = interface->grids->nx;
	nj = interface->grids->ny;
	nk = interface->grids->nz;
	
	nn = interface->lmBasis->length;
	
	basis_size = ni*nj*nk*nn;
	
	double *v_np;
	
	double *ulm1, *ulm2, *ulm;
	
	ulm = new double [basis_size];
	
	double *x_grid, *y_grid, *z_grid;
	
	x_grid = interface->grids->x_Grid;
	y_grid = interface->grids->y_Grid;
	z_grid = interface->grids->z_Grid;
	
#pragma omp parallel for default(shared) private(i,j,k,n,ulm1,ulm2,v_np) schedule(guided) collapse(3)
	for (i=0; i<ni; i++) {
		for (j=0; j<nj; j++) {
			for (k=0; k<nk; k++) {
				
				v_np = new double [nn];
				
				for (n=0; n<nn; n++) {
					v_np[n] = v_ipjkn[((n*nk + k)*nj + j)*ni + i];
				}
				
				ulm1 = calc_ulm(x_grid[i], y_grid[j], z_grid[k], v_np, interface, 0); //phi = acos(cos phi)
				ulm2 = calc_ulm(x_grid[i], y_grid[j], z_grid[k], v_np, interface, 1); //phi = 2PI - acos(cos phi)
				
				for (n=0; n<nn; n++) {
					ulm[((n*nk + k)*nj + j)*ni + i] = ulm1[n] + ulm2[n];
				}
				
				delete [] v_np;
				delete [] ulm1;
				delete [] ulm2;
			}
		}
	}
	
	
	
	return ulm;
}

double* Hv_5D_oneCompositeIndex(interfaceStor *interface, double *v_ipjkn) {
	int p;
	
	int ni, nj, nk, nn, basis_size;	
	ni = interface->grids->nx;
	nj = interface->grids->ny;
	nk = interface->grids->nz;
	
	nn = interface->lmBasis->length;
	
	basis_size = ni*nj*nk*nn;	
	
	double *Tv_ijkn;
	double *Vv_ijkn;
	double *Hv_ijkn = new double [basis_size];
	
	//Calculate the kinetic energy terms
	Tv_ijkn = Tv_5D_oneCompositeIndex(interface, v_ipjkn);
	
	//cout << "Tv finished" << endl;
	
	//Calculate the potential energy terms
	Vv_ijkn = Vv_5D_oneCompositeIndex(interface, v_ipjkn);
	
	//cout << "Vv finished" << endl;
	
	//Sum Tv_ijkn and Vv_ijkn to get Hv_ijkn
	//#pragma omp parallel for default(shared) private(p) schedule(guided) //Parallelization justified as no processor will access the same memory location at the same time (i.e. p is different for each processor)
	for (p=0; p<basis_size; p++) {
		Hv_ijkn[p] = Tv_ijkn[p] + Vv_ijkn[p]; 
	//	Hv_ijkn[p] = Tv_ijkn[p]; // Free system - only for debugging purposes
	}
	
	delete [] Tv_ijkn;
	delete [] Vv_ijkn;
	
	return Hv_ijkn;
}

void Hv_Prep_linRotCartSph(int argc, char **argv, generalStor *general_data, lanczosStor *lanczos_data) {
	
	interfaceStor *Hv_data;
	Hv_data = reinterpret_cast<interfaceStor*> (general_data);
	
	HvPrep_Internal(argc, argv, Hv_data, lanczos_data);
	
};

void Hv_linRotCartSph(int argc, char **argv, generalStor *general_data, lanczosStor *lanczos_data, double *vec, double *uec) {
	int i;
	double *uec1;
	
	interfaceStor *Hv_data;
	Hv_data = reinterpret_cast<interfaceStor*> (general_data);
	
	uec1 = Hv_5D_oneCompositeIndex(Hv_data, vec);
	
	for (i=0; i<lanczos_data->total_basis_size; i++) {
		uec[i] += uec1[i];
	}
	
	delete [] uec1;
	
};

/*
int main(int argc, char** argv) {
	//int l_max, thetaPoints, phiPoints;
	//	
	//	l_max = atoi(argv[1]);
	//	thetaPoints = atoi(argv[2]);
	//	phiPoints = atoi(argv[3]);
	//	
	//	tesseralTest(l_max, thetaPoints, phiPoints);
	
	interfaceStor *interface = new interfaceStor();
	lanczosStor *lanczos = new lanczosStor();
	
	HvPrep_Internal(argc, argv, interface);
	
	int n;
	
	//Hv
	int ni, nj, nk, nn;
	
	ni = interface->grids->nx;
	nj = interface->grids->ny;
	nk = interface->grids->nz;
	
	nn = interface->lmBasis->length;
	
	double *Hv_ijkn, *Hv_temp;
	double *v_ipjkn = new double [ni*nj*nk*nn];
	
	for (n=0; n<(ni*nj*nk*nn); n++) {
		v_ipjkn[n] = 1.0 / sqrt(double(ni*nj*nk*nn));
	}
	
	Hv_ijkn =  Hv_5D_oneCompositeIndex(interface, v_ipjkn);
	
	for (n=0; n<atoi(argv[2]); n++) {
		Hv_temp = Hv_ijkn;
		Hv_ijkn =  Hv_5D_oneCompositeIndex(interface, Hv_ijkn);
		
		delete [] Hv_temp;
	}
	
	cout << "Hv finished" << endl;
	
//	for (n=(ni*nj*nk*nn-10); n<(ni*nj*nk*nn); n++) {
//		cout <<	Hv_ijkn[n] << endl;
//	}
	
//	for (n=0; n<(ni*nj*nk*nn); n++) {
//		cout <<	Hv_ijkn[n] << endl;
//	}
	
	delete interface;
	
	delete [] Hv_ijkn;
	delete [] v_ipjkn;
	
	return 0;
}
*/

