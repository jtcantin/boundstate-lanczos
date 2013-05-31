#include "rotmat.h"

//THINK ABOUT A CUTOFF FOR SMALL VALUES!!

//Rotate using the Euler Angles (ZY'Z" convention)
MAT rotMat3D(double *EA) {
	double cp, sp, ct, st, cx, sx;
	
	MAT mat(3,3);
	
	//Assume ZY'Z" convention, with the angles azimuth=phi(0), polar=theta(1), intrinsic=chi(2)
	cp = cos(EA[0]); //phi
	sp = sin(EA[0]);
	ct = cos(EA[1]); //theta
	st = sin(EA[1]);
	cx = cos(EA[2]); //chi
	sx = sin(EA[2]);

	mat(0,0) = (cp*ct*cx) - (sp*sx);
	mat(0,1) = (-1*cp*ct*sx) - (sp*cx);
	mat(0,2) = cp*st;
	
	mat(1,0) = (sp*ct*cx) + (cp*sx);
	mat(1,1) = (-1*sp*ct*sx) + (cp*cx);
	mat(1,2) = sp*st;
	
	mat(2,0) = -1*st*cx;
	mat(2,1) = st*sx;
	mat(2,2) = ct;

	return mat;
}

//Rotate about the z-axis (i.e. on the x-y plane)
MAT rotMat3D_Z(double ang) {
	double cA, sA;

	MAT mat(3,3);

	cA = cos(ang);
	sA = sin(ang);
	
	//cout << cA << "," << sA << ";" << ang << endl;

	mat(0,0) = cA;
	mat(0,1) = -1*sA;
	mat(0,2) = 0.;

	mat(1,0) = sA;
	mat(1,1) = cA;
	mat(1,2) = 0.;

	mat(2,0) = 0.;
	mat(2,1) = 0.;
	mat(2,2) = 1.;
	
	return mat;
}

//Rotate about the y-axis (i.e. on the z-x plane)
MAT rotMat3D_Y(double ang) {
	double cA, sA;
	
	MAT mat(3,3);
	
	cA = cos(ang);
	sA = sin(ang);
	
	mat(0,0) = cA;
	mat(0,1) = 0.;
	mat(0,2) = sA;
	
	mat(1,0) = 0.;
	mat(1,1) = 1.;
	mat(1,2) = 0.;
	
	mat(2,0) = -1*sA;
	mat(2,1) = 0.;
	mat(2,2) = cA;
	
	return mat;
}

//Rotate about the x-axis (i.e. on the z-y plane)
MAT rotMat3D_X(double ang) {
	double cA, sA;
	
	MAT mat(3,3);
	
	cA = cos(ang);
	sA = sin(ang);
	
	mat(0,0) = 1.;
	mat(0,1) = 0.;
	mat(0,2) = 0.;
	
	mat(1,0) = 0.;
	mat(1,1) = cA;
	mat(1,2) = -1*sA;
	
	mat(2,0) = 0.;
	mat(2,1) = sA;
	mat(2,2) = cA;
	
	return mat;
}
