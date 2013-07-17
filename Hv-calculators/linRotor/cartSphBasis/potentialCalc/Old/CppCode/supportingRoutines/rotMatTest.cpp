#include "vectClass.h"
#include "rotmat.h"
#include <cstdlib>
#include <iostream>
#include <cmath>

#define PI 3.14159265359 //From MMTK

int main(int argc, char **argv) {
	VECT resVec(3);
	MAT rotMat(3,3);
	VECT vec(3);
	double EA[3];
	
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	rotMat = rotMat3D_Z(45.*PI/180.);
	resVec = rotMat*vec; //Should be (0.7071, 0.7071,0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 45 << " deg about Z-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	rotMat = rotMat3D_Z(-135.*PI/180.);
	resVec = rotMat*vec; //Should be (-0.7071, -0.7071,0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << -135 << " deg about Z-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	rotMat = rotMat3D_Z(90.*PI/180.);
	resVec = rotMat*vec; //Should be (0,1,0) THINK ABOUT A CUT OFF FOR SMALL VALUES
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 90 << " deg about Z-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	rotMat = rotMat3D_Z(120.*PI/180.);
	resVec = rotMat*vec; //Should be (-0.5, 0.8660,0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 120 << " deg about Z-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	
	///////////////////////////////////////////////
	cout << endl;
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	rotMat = rotMat3D_X(45.*PI/180.);
	resVec = rotMat*vec; //Should be (1,0,0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 45 << " deg about X-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	rotMat = rotMat3D_X(-135.*PI/180.);
	resVec = rotMat*vec; //Should be (1,0,0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << -135 << " deg about X-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	rotMat = rotMat3D_X(90.*PI/180.);
	resVec = rotMat*vec; //Should be (1,0,0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 90 << " deg about X-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	rotMat = rotMat3D_X(120.*PI/180.);
	resVec = rotMat*vec; //Should be (1,0,0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 120 << " deg about X-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	///////////////////////////////////////////////
	cout << endl;
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	rotMat = rotMat3D_Y(45.*PI/180.);
	resVec = rotMat*vec; //Should be (0.7071, 0, -0.7071)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 45 << " deg about Y-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	rotMat = rotMat3D_Y(-135.*PI/180.);
	resVec = rotMat*vec; //Should be (-0.7071, 0, 0.7071)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << -135 << " deg about Y-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	rotMat = rotMat3D_Y(90.*PI/180.);
	resVec = rotMat*vec; //Should be (0,0,-1) THINK ABOUT A CUT OFF FOR SMALL VALUES
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 90 << " deg about Y-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	rotMat = rotMat3D_Y(120.*PI/180.);
	resVec = rotMat*vec; //Should be (-0.5, 0, -0.8660)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 120 << " deg about Y-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << endl;
	cout << endl;
	cout << endl;
	cout << endl;
	vec(0) = 0.;
	vec(1) = 1.;
	vec(2) = 0.;
	rotMat = rotMat3D_Z(45.*PI/180.);
	resVec = rotMat*vec; //Should be (-0.7071, 0.7071,0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 45 << " deg about Z-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 1.;
	vec(2) = 0.;
	rotMat = rotMat3D_Z(-135.*PI/180.);
	resVec = rotMat*vec; //Should be (0.7071, -0.7071,0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << -135 << " deg about Z-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 1.;
	vec(2) = 0.;
	rotMat = rotMat3D_Z(90.*PI/180.);
	resVec = rotMat*vec; //Should be (-1,0,0) THINK ABOUT A CUT OFF FOR SMALL VALUES
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 90 << " deg about Z-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 1.;
	vec(2) = 0.;
	rotMat = rotMat3D_Z(120.*PI/180.);
	resVec = rotMat*vec; //Should be (-0.8660,-0.5,0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 120 << " deg about Z-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	
	///////////////////////////////////////////////
	cout << endl;
	cout << endl;
	vec(0) = 0.;
	vec(1) = 1.;
	vec(2) = 0.;
	rotMat = rotMat3D_X(45.*PI/180.);
	resVec = rotMat*vec; //Should be (0,0.7071,0.7071)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 45 << " deg about X-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 1.;
	vec(2) = 0.;
	rotMat = rotMat3D_X(-135.*PI/180.);
	resVec = rotMat*vec; //Should be (0,-0.7071,-0.7071)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << -135 << " deg about X-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 1.;
	vec(2) = 0.;
	rotMat = rotMat3D_X(90.*PI/180.);
	resVec = rotMat*vec; //Should be (0,0,1)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 90 << " deg about X-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 1.;
	vec(2) = 0.;
	rotMat = rotMat3D_X(120.*PI/180.);
	resVec = rotMat*vec; //Should be (0,-0.5,0.8660)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 120 << " deg about X-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	///////////////////////////////////////////////
	cout << endl;
	cout << endl;
	vec(0) = 0.;
	vec(1) = 1.;
	vec(2) = 0.;
	rotMat = rotMat3D_Y(45.*PI/180.);
	resVec = rotMat*vec; //Should be (0, 1, 0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 45 << " deg about Y-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 1.;
	vec(2) = 0.;
	rotMat = rotMat3D_Y(-135.*PI/180.);
	resVec = rotMat*vec; //Should be (0, 1, 0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << -135 << " deg about Y-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 1.;
	vec(2) = 0.;
	rotMat = rotMat3D_Y(90.*PI/180.);
	resVec = rotMat*vec; //Should be (0, 1, 0) THINK ABOUT A CUT OFF FOR SMALL VALUES
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 90 << " deg about Y-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 1.;
	vec(2) = 0.;
	rotMat = rotMat3D_Y(120.*PI/180.);
	resVec = rotMat*vec; //Should be (0, 1, 0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 120 << " deg about Y-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << endl;
	cout << endl;
	cout << endl;
	cout << endl;
	vec(0) = 0.;
	vec(1) = 0.;
	vec(2) = 1.;
	rotMat = rotMat3D_Z(45.*PI/180.);
	resVec = rotMat*vec; //Should be (0,0,1)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 45 << " deg about Z-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 0.;
	vec(2) = 1.;
	rotMat = rotMat3D_Z(-135.*PI/180.);
	resVec = rotMat*vec; //Should be (0,0,1)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << -135 << " deg about Z-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 0.;
	vec(2) = 1.;
	rotMat = rotMat3D_Z(90.*PI/180.);
	resVec = rotMat*vec; //Should be (0,0,1) THINK ABOUT A CUT OFF FOR SMALL VALUES
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 90 << " deg about Z-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 0.;
	vec(2) = 1.;
	rotMat = rotMat3D_Z(120.*PI/180.);
	resVec = rotMat*vec; //Should be (0,0,1)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 120 << " deg about Z-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	
	///////////////////////////////////////////////
	cout << endl;
	cout << endl;
	vec(0) = 0.;
	vec(1) = 0.;
	vec(2) = 1.;
	rotMat = rotMat3D_X(45.*PI/180.);
	resVec = rotMat*vec; //Should be (0,-0.7071,0.7071)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 45 << " deg about X-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 0.;
	vec(2) = 1.;
	rotMat = rotMat3D_X(-135.*PI/180.);
	resVec = rotMat*vec; //Should be (0,0.7071,-0.7071)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << -135 << " deg about X-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 0.;
	vec(2) = 1.;
	rotMat = rotMat3D_X(90.*PI/180.);
	resVec = rotMat*vec; //Should be (0,-1,0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 90 << " deg about X-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 0.;
	vec(2) = 1.;
	rotMat = rotMat3D_X(120.*PI/180.);
	resVec = rotMat*vec; //Should be (0,-0.8660,-0.5)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 120 << " deg about X-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	///////////////////////////////////////////////
	cout << endl;
	cout << endl;
	vec(0) = 0.;
	vec(1) = 0.;
	vec(2) = 1.;
	rotMat = rotMat3D_Y(45.*PI/180.);
	resVec = rotMat*vec; //Should be (0.7071,0,0.7071)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 45 << " deg about Y-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 0.;
	vec(2) = 1.;
	rotMat = rotMat3D_Y(-135.*PI/180.);
	resVec = rotMat*vec; //Should be (-0.7071,0,-0.7071)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << -135 << " deg about Y-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 0.;
	vec(2) = 1.;
	rotMat = rotMat3D_Y(90.*PI/180.);
	resVec = rotMat*vec; //Should be (1, 0, 0) THINK ABOUT A CUT OFF FOR SMALL VALUES
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 90 << " deg about Y-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 0.;
	vec(1) = 0.;
	vec(2) = 1.;
	rotMat = rotMat3D_Y(120.*PI/180.);
	resVec = rotMat*vec; //Should be (0.8660,0,-0.5)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 120 << " deg about Y-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << endl;
	cout << endl;
	cout << endl;
	cout << endl;
	vec(0) = 1.;
	vec(1) = 1.;
	vec(2) = 1.;
	rotMat = rotMat3D_Z(45.*PI/180.);
	resVec = rotMat*vec; //Should be (0,1.4142,1)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angle: " << 45 << " deg about Z-axis" << endl;
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << endl;
	cout << endl;
	cout << endl;
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	EA[0] = 0.*PI/180.;
	EA[1] = 0.*PI/180.;
	EA[2] = 0.*PI/180.;
	rotMat = rotMat3D(EA);
	resVec = rotMat*vec; //Should be (1,0,0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angles: phi = " << EA[0]*180/PI << " deg, theta = " << EA[1]*180/PI << " deg, chi = " << EA[2]*180/PI << endl; 
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	EA[0] = 135.*PI/180.;
	EA[1] = 0.*PI/180.;
	EA[2] = 0.*PI/180.;
	rotMat = rotMat3D(EA);
	resVec = rotMat*vec; //Should be (-0.7071,0.7071,0)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angles: phi = " << EA[0]*180/PI << " deg, theta = " << EA[1]*180/PI << " deg, chi = " << EA[2]*180/PI << endl; 
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	EA[0] = 135.*PI/180.;
	EA[1] = 45.*PI/180.;
	EA[2] = 0.*PI/180.;
	rotMat = rotMat3D(EA);
	resVec = rotMat*vec; //Should be (-0.5,0.5,-0.7071)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angles: phi = " << EA[0]*180/PI << " deg, theta = " << EA[1]*180/PI << " deg, chi = " << EA[2]*180/PI << endl; 
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;
	
	cout << endl;
	vec(0) = 1.;
	vec(1) = 0.;
	vec(2) = 0.;
	EA[0] = 135.*PI/180.;
	EA[1] = 45.*PI/180.;
	EA[2] = 45.*PI/180.;
	rotMat = rotMat3D(EA);
	resVec = rotMat*vec; //Should be (-0.85355,-0.14645,-0.5)
	cout << "Input vector: " << "(" << vec(0) << "," << vec(1) << "," << vec(2) << ")" << endl;
	cout << "Rotation Angles: phi = " << EA[0]*180/PI << " deg, theta = " << EA[1]*180/PI << " deg, chi = " << EA[2]*180/PI << endl; 
	cout << "Output vector: " << "(" << resVec(0) << "," << resVec(1) << "," << resVec(2) << ")" << endl;

	return 0;
}

/*
 for (i=0; i<3; i++) {
 for (j=0; j<3; j++) {
 cout << "(" << i << "," << j << ")" << " = " << mat(i,j) << endl;
 }
 }*/