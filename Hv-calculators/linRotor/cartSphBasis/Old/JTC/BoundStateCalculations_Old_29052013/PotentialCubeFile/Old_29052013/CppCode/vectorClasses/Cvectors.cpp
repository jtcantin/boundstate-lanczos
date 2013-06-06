#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Cvectors.h"

using namespace std;

// ND Cartesian vectors #######################################################

//Constructor
CNvect::CNvect () {
	coor = new double [1];
	num_dim = 1;
}

CNvect::CNvect (int dim) {
	int i;
	num_dim = dim;
	coor = new double [num_dim];
	
	
	for (i=0; i<num_dim; i++) {
		coor[i] = 0.;
	}
}

CNvect::CNvect (double* new_coor, int dim) {
	int i;
	num_dim = dim;
	coor = new double [num_dim];
	
	for (i=0; i<num_dim; i++) {
		coor[i] = new_coor[i];
	}
}

/*
//Destructor
CNvect::~CNvect(){
	delete[] coor;
}
*/
 
void CNvect::setCoor(double* new_coor, int dim)
{
	int i;
	if (dim != num_dim) {
		cerr << "Error, incorrect number of dimensions passed during CNvect.setCor() call." << endl;
		exit(1);
	}
	
	for (i=0; i<dim; i++) {
		coor[i] = new_coor[i];
	}
}

// NOTE: Resizing the number of dimensions clears out the vector.
void CNvect::setDim(int dim){
	int i;
	num_dim = dim;
	//delete [] coor; //Need to dealocate the old memory, or there is a memory leak
	coor = new double [num_dim];
	
	
	for (i=0; i<num_dim; i++) {
		coor[i] = 0.;
	}
}

//Addition
CNvect CNvect::operator + (const CNvect& bVec) const
{
	int i;
	//CNvect cVec(this->num_dim);
	CNvect cVec;
	cout << "here" << endl;
	cout << this->num_dim << endl;
	cout << bVec.num_dim << endl;
	
	if (this->num_dim != bVec.num_dim) {
		cerr << "Error, you can only add vectors of the same dimension." << endl;
		exit(1);
	}
	
	for (i=0; i<(this->num_dim); i++) {
		cVec.coor[i] = this->coor[i] + bVec.coor[i];
	}
	
	return cVec;
}

//Subtraction
CNvect CNvect::operator - (const CNvect& bVec) const
{
	int i;
	CNvect cVec(this->num_dim);
	
	if (this->num_dim != bVec.num_dim) {
		cerr << "Error, you can only subtract vectors of the same dimension." << endl;
		exit(1);
	}
	
	for (i=0; i<(this->num_dim); i++) {
		cVec.coor[i] = this->coor[i] - bVec.coor[i];
	}
	
	return cVec;
}

//Scalar Multiplication
CNvect operator * (int a, const CNvect& Vec)
{
	int i;
	CNvect newVec(Vec.num_dim);
	
	for (i=0; i<(Vec.num_dim); i++) {
		newVec.coor[i] = Vec.coor[i] * a;
	}
	
	return newVec;
}

CNvect operator * (const CNvect& Vec, int a)
{
	int i;
	CNvect newVec(Vec.num_dim);
	
	for (i=0; i<(Vec.num_dim); i++) {
		newVec.coor[i] = Vec.coor[i] * a;
	}
	
	return newVec;
}

//Dot product
double CNvect::operator * (const CNvect& bVec) const
{
	int i;
	double dotProd = 0.0;
	
	if (this->num_dim != bVec.num_dim) {
		cerr << "Error, you can only dot vectors of the same dimension." << endl;
		exit(1);
	}
	
	for (i=0; i<(this->num_dim); i++) {
		dotProd += this->coor[i] * bVec.coor[i];
	}
	
	return dotProd;
}

//Assignment NOTE: This will overwrite the size of the assigned to vector, without warning.
CNvect& CNvect::operator = (const CNvect& bVec)
{
	int i;
	
	num_dim = bVec.num_dim;
	coor = new double [num_dim];
	
	
	for (i=0; i<num_dim; i++) {
		coor[i] = bVec.coor[i];
	}
	
	return *this;
}

//Get the Euclidean norm
double CNvect::eNorm ()  const
{
	double sumSquare;
	sumSquare = (*this) * (*this);
	return sqrt(sumSquare);
}

//Print out the coordinates
void CNvect::printVec ()  const
{
	int i;
	for (i=0; i<(this->num_dim); i++) {
		cout << "Dim " << i << ": " << this->coor[i] << endl;
	}
}



// 3D Cartesian vectors #######################################################

//Constructor
C3vect::C3vect (double a, double b, double c)
{
	x = a;
	y = b;
	z = c;
}

void C3vect::setCor(double a, double b, double c)
{
	x = a;
	y = b;
	z = c;
}

//Addition
C3vect C3vect::operator + (const C3vect& bVec) const
{
	C3vect cVec;
	cVec.x = this->x + bVec.x;
	cVec.y = this->y + bVec.y;
	cVec.z = this->z + bVec.z;
	
	return (cVec);
}

//Subtraction
C3vect C3vect::operator - (const C3vect& bVec) const
{
	C3vect cVec;
	cVec.x = this->x - bVec.x;
	cVec.y = this->y - bVec.y;
	cVec.z = this->z - bVec.z;
	
	return (cVec);
}

//Scalar Multiplication
C3vect operator * (int a, const C3vect& Vec)
{
	C3vect newVec;
	newVec.x = Vec.x * a;
	newVec.y = Vec.y * a;
	newVec.z = Vec.z * a;
	
	return newVec;
}

C3vect operator * (const C3vect& Vec, int a)
{
	C3vect newVec;
	newVec.x = Vec.x * a;
	newVec.y = Vec.y * a;
	newVec.z = Vec.z * a;
	
	return newVec;
}

//Dot product
double C3vect::operator * (const C3vect& bVec) const
{
	double dotProd = 0.0;
	dotProd += this->x * bVec.x;
	dotProd += this->y * bVec.y;
	dotProd += this->z * bVec.z;
	
	return dotProd;
}

//Assignment
C3vect& C3vect::operator = (const C3vect& bVec)
{
	x = bVec.x;
	y = bVec.y;
	z = bVec.z;
	
	return *this;
}

//Get the Euclidean norm
double C3vect::eNorm ()  const
{
	double sumSquare;
	sumSquare = (*this) * (*this);
	return sqrt(sumSquare);
}

//Print out the coordinates
void C3vect::printVec ()  const
{
	cout << "x: " << this->x << endl;
	cout << "y: " << this->y << endl;
	cout << "z: " << this->z << endl;
}

/* Testing Code for the C3vect structure
 C3vect vecA, vecB, vecC, vecD;
 double length, dotted;
 
 vecA.setCor(atof(argv[1]), atof(argv[2]),atof(argv[3]));
 //vecA.setCor(1,1,1); //values in the comments are for (1,1,1); (2.5, 3.6, 4.2)
 cout << "Vector A: " << endl;
 vecA.printVec(); //x=1, y=1, z=1; x=2.5, y=3.6, z=4.2
 cout << endl;
 
 
 vecB = vecA;
 cout << "Vector B: " << endl;
 vecB.printVec(); //x=1, y=1, z=1; x=2.5, y=3.6, z=4.2
 cout << endl;
 
 vecB = vecB*2;
 cout << "Vector B*2: " << endl;
 vecB.printVec(); //x=2, y=2, z=2; x=5, y=7.2, z=8.4
 cout << endl;
 
 vecB = 2*vecB;
 cout << "Vector 2*B: " << endl;
 vecB.printVec(); //x=4, y=4, z=4; x=10, y=14.4, z=16.8
 cout << endl;
 
 length = vecA.eNorm();
 cout << "Vector A length: " << endl;
 cout << length << endl; //sqrt(3) ~= 1.732050808; 6.07042
 cout << endl; 
 
 dotted = vecA * vecB;
 cout << "Vector A*B: " << endl;
 cout << dotted << endl; //12; 147.4
 cout << endl;
 
 dotted = vecB * vecA;
 cout << "Vector B*A: " << endl;
 cout << dotted << endl; //12; 147.4
 cout << endl;
 
 vecC = vecB + vecA;
 cout << "Vector C=B+A: " << endl;
 vecC.printVec(); //x=5, y=5, z=5; x=12.5, y=18, z=21
 cout << endl;
 
 vecC = vecA + vecB;
 cout << "Vector C=A+B: " << endl;
 vecC.printVec(); //x=5, y=5, z=5; x=12.5, y=18, z=21
 cout << endl;
 
 vecD = vecC - vecA;
 cout << "Vector D=C-A: " << endl;
 vecD.printVec(); //x=4, y=4, z=4; x=10, y=14.4, z=16.8
 cout << endl;
 
 vecD = vecA - vecC;
 cout << "Vector D=A-C: " << endl;
 vecD.printVec(); //x=-4, y=-4, z=-4; x=-10, y=-14.4, z=-16.8
 cout << endl;
 */







// 2D Cartesian vectors #######################################################

//Constructor
C2vect::C2vect (double a, double b)
{
	x = a;
	y = b;
}

void C2vect::setCor(double a, double b)
{
	x = a;
	y = b;
}

//Addition
C2vect C2vect::operator + (const C2vect& bVec) const
{
	C2vect cVec;
	cVec.x = this->x + bVec.x;
	cVec.y = this->y + bVec.y;
	
	return (cVec);
}

//Subtraction
C2vect C2vect::operator - (const C2vect& bVec) const
{
	C2vect cVec;
	cVec.x = this->x - bVec.x;
	cVec.y = this->y - bVec.y;
	
	return (cVec);
}

//Scalar Multiplication
C2vect operator * (int a, const C2vect& Vec)
{
	C2vect newVec;
	newVec.x = Vec.x * a;
	newVec.y = Vec.y * a;
	
	return newVec;
}

C2vect operator * (const C2vect& Vec, int a)
{
	C2vect newVec;
	newVec.x = Vec.x * a;
	newVec.y = Vec.y * a;
	
	return newVec;
}

//Dot product
double C2vect::operator * (const C2vect& bVec) const
{
	double dotProd = 0.0;
	dotProd += this->x * bVec.x;
	dotProd += this->y * bVec.y;
	
	return dotProd;
}

//Assignment
C2vect& C2vect::operator = (const C2vect& bVec)
{
	x = bVec.x;
	y = bVec.y;
	return *this;
}

//Get the Euclidean norm
double C2vect::eNorm ()  const
{
	double sumSquare;
	sumSquare = (*this) * (*this);
	return sqrt(sumSquare);
}

//Print out the coordinates
void C2vect::printVec ()  const
{
	cout << "x: " << this->x << endl;
	cout << "y: " << this->y << endl;
}

/* Testing Code for the C2vect structure
 C2vect vecA, vecB, vecC, vecD;
 double length, dotted;
 
 vecA.setCor(atof(argv[1]), atof(argv[2]));
 //vecA.setCor(1,1); //values in the comments are for (1,1)
 cout << "Vector A: " << endl;
 vecA.printVec(); //x=1, y=1
 cout << endl;
 
 
 vecB = vecA;
 cout << "Vector B: " << endl;
 vecB.printVec(); //x=1, y=1
 cout << endl;
 
 vecB = vecB*2;
 cout << "Vector B*2: " << endl;
 vecB.printVec(); //x=2, y=2
 cout << endl;
 
 vecB = 2*vecB;
 cout << "Vector 2*B: " << endl;
 vecB.printVec(); //x=4, y=4
 cout << endl;
 
 length = vecA.eNorm();
 cout << "Vector A length: " << endl;
 cout << length << endl; //sqrt(2) ~= 1.414213562373095
 cout << endl; 
 
 dotted = vecA * vecB;
 cout << "Vector A*B: " << endl;
 cout << dotted << endl; //8
 cout << endl;
 
 dotted = vecB * vecA;
 cout << "Vector B*A: " << endl;
 cout << dotted << endl; //8
 cout << endl;
 
 vecC = vecB + vecA;
 cout << "Vector C=B+A: " << endl;
 vecC.printVec(); //x=5, y=5
 cout << endl;
 
 vecC = vecA + vecB;
 cout << "Vector C=A+B: " << endl;
 vecC.printVec(); //x=5, y=5
 cout << endl;
 
 vecD = vecC - vecA;
 cout << "Vector D=C-A: " << endl;
 vecD.printVec(); //x=4, y=4
 cout << endl;
 
 vecD = vecA - vecC;
 cout << "Vector D=A-C: " << endl;
 vecD.printVec(); //x=-4, y=-4
 cout << endl;
 */


