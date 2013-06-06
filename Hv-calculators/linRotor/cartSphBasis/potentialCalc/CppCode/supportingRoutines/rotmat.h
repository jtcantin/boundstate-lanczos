#include "vectClass.h" //This is to select which linear algebra class should be used.
#include <cstdlib>
#include <cmath>
using namespace std;

#ifndef ROTMAT_H
#define ROTMAT_H

MAT rotMat3D(double *EA);

MAT rotMat3D_Z(double ang);

MAT rotMat3D_Y(double ang);

MAT rotMat3D_X(double ang);

#endif