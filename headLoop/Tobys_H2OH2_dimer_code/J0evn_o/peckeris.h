#include "BF.h"
// JAcobi basis JCP, 101, 1343, 1994
//const double kB=3.1668288610848352283e-6; // Hartree/Kelvin
static const double kB = 3.1668153e-6; // Hartree/Kelvin
void PIE(vector &v,int size);
void PIA1(vector &v,int size);
void PIA2(vector &v,int size);
vector HpsiJac(matrix &ddr,diagmat &G11,diagmat &G22,diagmat &G33,diagmat &G12,
	       diagmat &G13,diagmat &G23,diagmat &VaplusV,int size,vector &v);
void HpsiJacII(matrix &ddr,diagmat &G11,diagmat &G22,diagmat &G33,diagmat &G12,
		 diagmat &G13,diagmat &G23,diagmat &VaplusV,int size,vector &v,
		 vector &u,matrix &Vmat1,matrix &Vmat2, matrix &Vmat3,matrix &W1,
		 matrix &W2,matrix &W3);
void Hpsipinned(matrix &ddr,diagmat &G11,diagmat &G22,diagmat &G33,diagmat &G12,
		 diagmat &G13,diagmat &G23,diagmat &VaplusV,int size,vector &v,
		 vector &u,matrix &Vmat1,matrix &Vmat2, matrix &Vmat3,matrix &W1,
		 matrix &W2,matrix &W3);
void Hpsi2d(matrix &ddr,diagmat &G22,diagmat &G33,
	    diagmat &G23,diagmat &VaplusV,int size,vector &v,
	    vector &u,matrix &Vmat1,matrix &Vmat2, matrix &Vmat3,matrix &W1,
	    matrix &W2,matrix &W3);
void HpsiJacCont(matrix &ddr,diagmat &G11,diagmat &G22,diagmat &G33,diagmat &G12,
		 diagmat &G13,diagmat &G23,diagmat &VaplusV,int size,vector &v,
		 vector &u,matrix &Vmat1,matrix &Vmat2, matrix &Vmat3,matrix &W1,
		 matrix &W2,matrix &W3,int nzeroes);
double zn(double z,double dn);
matrix BJacobi(double a,double b,int size);
double Mygamma(double x);
matrix xJacobi(double a,double b,int size);
double q(double dn,double a,double b);
double d(double dn,double a,double b);
double cdown(double dm,double dn,double a,double b);
double cboth(double dm,double b,double dn,double s);
matrix DJacobi(double a,double b,int size);
double dmn(double dm,double dn,double a,double b);
double gauss(double r,double mu,double w,double r0);
double gaussnorm(double r,double alpha,double r0);
double del(double a);
double delbar(double a);
double delta(int i,int ip);
vector Hpsi(matrix &ddr2,diagmat &Rinv,diagmat &R2inv,diagmat &R4inv2,matrix &Rinv4ij,matrix &ddr,
	    diagmat &extraV,diagmat &V,int size,vector &v,double mass);
vector HpsiBLAS(diagmat &Rinv,diagmat &R2inv,diagmat &R4inv2,
		matrix &Rinv4ij,matrix &ddr,matrix &ddr2invrddr,
		diagmat &extraV,diagmat &V,int size,vector &v,double mass);
vector HpsiHERM(matrix &ddr2,matrix &ddr,diagmat &c1,diagmat &c2,diagmat &V,int size,vector &v,double mass);
void addtobasis(vector &v,matrix &B,int n,int col);
double Aziz(double r);
double lj(double r, double sigma, double epsilon);
double HeHminus(double r);
double Hebuckypot(double r,int Ncarbons);
