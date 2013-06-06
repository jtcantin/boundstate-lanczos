#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Cvectors.h"

using namespace std;

int main(int argc, char** argv){
	CNvect vecA, vecB, vecC, vecD, vecE;
	double length, dotted, *coor;
	int i, dim;
	
	dim = atoi(argv[1]);
	
	coor = new double [dim];
	
	for (i=0; i<dim; i++) {
		coor[i] = atof(argv[i+2]);
	}
	
	vecA.setDim(dim);
	vecB.setDim(dim);
	vecC.setDim(dim);
	vecD.setDim(dim);
	
	vecA.setCor(coor,dim);
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
	
	vecD = vecC * vecE;
	
	return 0;
}