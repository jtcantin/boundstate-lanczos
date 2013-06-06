#include "inter.h"
static const double  hatocm=219474.63067;
//                          219474.63067
// garshuk and light value below
//static const double  hatocm=219474.629;
static const double atob=1./0.529177249;   
//static const double autoamu=5.4857989586762187e-4;                          
static const double autoamu=5.48579903e-4;
// garshuk and light value below
//static const double autoamu=1./1822.88853006;
static const double MASS=1.;
static const int NTRAJECTORIES=1000;
static const int EQUILIBRATION=1000;
// amu    = 1.6605402e-27 Kg
// me     = 9.1093897e-31 Kg
// me/amu = 5.4857990e-4
//        = 5.48579903e-4 electron molar mass in Kg/mol
static const double R0=atob*3.820;
static const double rminimum=atob*1.278;
//static const double rminimum=2.415;
static const double Hmass=1.007825035/autoamu;
static const double Clmass=34.968852721/autoamu;
static const double hclreducedmass=Hmass*Clmass/(Hmass+Clmass);
//static const double PI=M_PI;                                          
EXTERN void pothcl2_(double *reqq,double *v);
EXTERNC void gauleg(double x1,double x2,double *x,double *w,int n);
EXTERN void FORTRAN(gaulegf)(double *x1,double *x2,double *x,double *w,int *n);
EXTERNC double plgndr(int j,int m,double x);
// prototypes
void contourout(diagmat &pot2d,vector &grid1d1,vector &grid1d2,int nsize1d1,int nsize1d2,char *filename);
void psitout(cvector &psit,vector &grid1d1,vector &grid1d2,int nsize1d1,int nsize1d2,char *filename);
vector thegrid(int nsize1d,double length,double center);
vector cosdvr(int nsize1d,double length,double center);
vector sinedvr(int nsize1d,double length,double center);
matrix tmatcosdvr(int nsize1d,vector &grid1d,double length,double center,double mass,double r0,char parity3);
matrix tmatsinedvr(int nsize1d,vector &grid1d,double length,double center,double mass,double r0,char parity3);
vector thetagrid(int nsize);
diagmat pot(vector &grid1d1,vector &grid1d2,int nsize1d1,int nsize1d2,double r0,double R,char coordinates,double De);
matrix Tmat1D(int nsize,char parity,double mass,double r0,double length,char coordinates,vector &grid);
matrix tmatleg(int nsize,double mass,double r0);
cvector initialpsi(vector &grid1d1,vector &grid1d2,int nsize1d1,int nsize1d2,double width,double center1,double center2);
double potfunc(double x1,double x2,double r0,double R,char coordinates,double De);
complex g(double Tau,double tdt);
double Pj0(int j,double x);
vector phidvr(int nsize,double length,double center,char parity);
matrix tmatphidvr(int nsize,vector &grid,double length,double center,double mass,double r0,char parity,matrix &T);
inline double Power(double x,double y) {return pow(x,y);}
inline double Cos(double x) {return cos(x);}
inline double Sin(double x) {return sin(x);}
inline double Sqrt(double x) {return sqrt(x);}
inline double krdel(int i,int j) {
  if (i ==j) return 1.;
  else return 0.;
}
void readinput(int argc,char **argv,int &nsize1d1,int &nsize1d2,int &nsize2d,
	       double &r0,double &R,double &mass,double &totaltime,int &steps,
	       double &width,double &x01,double &x02,int &nenergies,
	       double &widthfactor,char &coordinates,char &exactpropagation,
	       char &parity1,char &parity2,
	       double &length1,double &origin1,
	       double &length2,double &origin2,
	       int &mode1,int & mode2,char &modeselectedinitialstate);
diagmat potang(vector &grid1,vector &grid2,vector &grid3,int nsize1,int nsize2,
	    double r0,double R,double De);
double potfuncang(double t1,double t2,double phi,double r0,double R,double De);
void input3d(int argc,char **argv,int &nsize1,int &nsize2,int &nsize3d,
	       double &r0,double &R,double &mass,double &totaltime,int &steps,
	       double &width,double &x01,double &x02,int &nenergies,
	       double &widthfactor,char &exactpropagation,
	       char &parity1,char &parity2,int &mode1a,int &mode1b,int &mode2,
	       char &modeselectedinitialstate);
matrix phipp(int nsize,vector &grid);
matrix phippsym(int nsize1d,vector &grid1d,char parity);
void scfpropagation3d(int nsize1,int nsize2,diagmat &pot3d,matrix &Ttheta1,
matrix &Ttheta2,diagmat &I1,diagmat &I2,matrix &Tphi,vector &grid1,
vector &grid2,vector &grid3,diagmat &pot2d,diagmat &potphi,double EqI,
int mode1a,int mode1b,int mode2,char modeselectedinitialstate,int steps,
double emin,
double width,double widthfactor,int nenergies,vector &gs1d,vector &gs2d,
matrix &HU,vector &ev,diagmat &pot1da,diagmat &pot1db);
void scfpropagation4d(int nsize1,int nsize2,matrix &Ttheta1,
matrix &Ttheta2,diagmat &I1,diagmat &I2,matrix &Tphi,vector &grid1,
vector &grid2,vector &grid3,double EqI,
int mode1a,int mode1b,int mode2,int steps,
double emin,
double width,double widthfactor,int nenergies,double initialR,double r0,
double De,matrix &H,vector &eval,int Rpoints,double dt);
void input4d(int argc,char **argv,int &nsize1,int &nsize2,int &nsize3d,
	     double &r0,double &R,double &mass,double &totaltime,int &steps,
	     double &width,double &x01,double &x02,int &nenergies,
	     double &widthfactor,char &exactpropagation,
	     char &parity1,char &parity2,int &mode1a,int &mode1b,int &mode2,
	   
  char &modeselectedinitialstate,int &Rpoints);
EXTERN int time_();
diagmat Rpotavg(diagmat &pot4d,int nsize3d,int Rpoints,cvector &psiR);
diagmat getpotR(diagmat &pot4d,int index,int Rpoints,int nsize3d);
cvector gsofR(diagmat &potR,int Rpoints,double deltaR,double dimermass,
int mode2,double width);
diagmat computepotR(diagmat &pot4d,int nsize3d,cvector &psit,int Rpoints);
cmatrix constructUTR(int Rpoints,double deltaR,double dimermass,double dt);
matrix kineticCWDVR(int nsize,double dx, double mass);
diagmat getpot(diagmat &pot4d,int Rindex,int nsize3d);
void classicaltrajectory(double &R,double &v,double mass,double dt,
			 double &force,diagmat &pot4d,int nsize1,
			 int nsize2,cvector &psit,
			 vector &Rgrid,int Rpoints);
double  computeforce(double R,diagmat &potR,vector &Rgrid,int Rpoints);
void LanczosU(matrix &T,diagmat &V,cvector &v,int size,double dt);
void lancbis(int niter,vector &eval,vector &evalerr,double elmin,
double elmax,int &ngood,const vector& alpha,const vector& beta,
const vector& beta2);
EXTERN void FORTRAN(bisec)(double *lalpha,double *lbeta,double *lbeta2,
int *,double *wrk4,int *mp2,double *wrk2,int *,
double *elmin,double *elmax,int *ndis);
EXTERN void FORTRAN(isoev)(double *vs,int *mp,int *ndis);
EXTERN void FORTRAN(inverr)(double *lalpha,double *lbeta,double *wrk1,
double *wrk2,double *wrk3,int *niter,double *vs,int *mp,double *err,int
*ndis);
EXTERN void FORTRAN(trivec)(double *lalpha,double *lbeta,double *lbeta2,
			    double *wrk1,double *wrk2,int *niter,
			    double *eval,int *ngood,double *evtr,double *mamin);
void Lanczos(matrix &TR,matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
	     diagmat &I2,matrix &Tphi,diagmat &pot4d,vector &v,int sizeR, 
	     int size1,int size2,int niter,double emin,double emax);
vector Hpsi(matrix &TR,matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
	    diagmat &I2,matrix &Tphi,diagmat &pot4d,vector &v,int Rpoints, 
	    int nsize1,int nsize2);
cvector cHpsiang(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
	    diagmat &I2,matrix &Tphi,diagmat &pot3d,cvector &v, 
	    int nsize1,int nsize2);
void LanczosUang(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
		 diagmat &I2,matrix &Tphi,diagmat &V,cvector &v,
		 int size1,int size2,double dt);
void LanczosUR(matrix &T,diagmat &V,cvector &v,int size,double dt);
double averageV(double deltaR,double R,diagmat &potR,
		vector &Rgrid,int Rpoints);
void GaussianWP(double &R,double &v,double mass,double dt,
		double &force,diagmat &potR,vector &Rgrid,
		int Rpoints,complex &A,complex &gamma);
complex GCorr(double initialR,double R,double velocity0,double velocity,
	      double dimermass,vector &Rgrid,int Rpoints,complex A0,
	      complex A,complex gamma0,complex gamma);
cvector  computewavepacket(complex gamma,complex A,double R,double p,
			   vector &Rgrid,int Rpoints);
double  computeforceconst(double R,diagmat &potR,vector &Rgrid,int Rpoints);
double  getforce(double R,vector &grid1,vector &grid2,vector &grid3,
		 int nsize1,int nsize2,cvector &psit,double r0,double De);
void GWP(double &R,double &v,double mass,double dt,
	 double &force,vector &grid1,vector &grid2,
	 vector &grid3,int nsize1,int nsize2,cvector &psit,
	 complex &A,complex &gamma,double r0,double De);
double  getforceconst(double R,vector &grid1,vector &grid2,vector &grid3,
		 int nsize1,int nsize2,cvector &psit,double r0,double De);
double avgV(double R,vector &grid1,vector &grid2,vector &grid3,
		 int nsize1,int nsize2,cvector &psit,double r0,double De);
double avgV(double R,vector &grid1,vector &grid2,vector &grid3,
		 int nsize1,int nsize2,cvector &psit,double r0,double De);
void LanczosU4D(matrix &TR,matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
		 diagmat &I2,matrix &Tphi,diagmat &V,cvector &v,
		 int size1,int size2,int sizeR,double dt);
complex phasepropang(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
		     diagmat &I2,matrix &Tphi,cvector &psi,
		     int size1,int size2,double dt);
cvector corrLanczosU4D(matrix &TR,matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
		       diagmat &I2,matrix &Tphi,diagmat &V,cvector &v,
		       int size1,int size2,int sizeR,double dt,
		       double totaltime,int niter);
diagmat potangsym(vector &grid1,vector &grid2,vector &grid3,int nsize1,
		  int nsize2,double r0,double R,double De);
matrix lmatleg(int nsize);
matrix potfbrsym(int nsize,vector &grid,diagmat &pot);
matrix phipsym(int nsize1d,vector &grid1d,char parity);
double Pj1(int j,double x);
matrix thetapdvr(int nsize);
void LanczosE(matrix &TR,
	      matrix &Ttheta1,matrix &Ttheta2,matrix &Ttheta1R,
	      matrix &Ttheta2R,
	      diagmat &Inertia,
	      matrix &Tphi,matrix &cosphipp,diagmat &cosphi,matrix &sinphip,
	      matrix &cotmat1,matrix &cotmat2,matrix &delcot1,
	      matrix &delcot2,diagmat &InertiaMR,diagmat &InertiaMR2,
	      diagmat &pot4d,vector &v,int sizeR, 
	      int nsize1,int nsize2,int niter,double emin,double emax);
matrix lmatleg1(int nsize);
matrix sinphipsym(int nsize1d,vector &grid1d,char parity);
matrix cosphippsym(int nsize1d,vector &grid1d,char parity);
vector symdvr(int nsize,char parity);
matrix cotdvr(int nsize);
matrix delcotdvr(int nsize);
matrix sin2dvr(int nsize);
void LanczosUangfast(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
		     diagmat &I2,matrix &Tphi,diagmat &V,cvector &v,
		     int size1,int size2,double dt,int niter,
		     cmatrix &LanczosVectors,cvector &r,
		     cvector &work3,vector &alpha,
		     vector &beta,matrix &H,cdiagmat &ev3,
		     cvector &psiL,vector &vRe,
		     vector &vIm,vector &uRe,vector &wRe,vector &uIm,
		     vector &wIm,cvector &u
		     );
void Lanczos3D(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
	       diagmat &I2,matrix &Tphi,diagmat &V,vector &v,
	       int nsize1,int nsize2,int niter,double emin,double emax);
vector thetadvra(int nsize,matrix &T);
matrix tmatlega(int nsize,double mass,double r0,matrix &lmat);
cvector corrLanczosU3D(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
		 diagmat &I2,matrix &Tphi,diagmat &V,cvector &v,
		 int size1,int size2,double dt,double totaltime,int niter);
cvector Hpsi2dfort(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
		   diagmat &I2,double Tphiavg,diagmat &pot2d,
		   cvector &psi2d,int nsize1,int nsize2);
void addtobasis(cvector &v,cmatrix &B,int n,int col);
diagmat potangsymsak(vector &grid1,vector &grid2,vector &grid3,int nsize1,
		  int nsize2,  double r0,double R,double De);
double potfuncangsak(double t1,double t2,double phi,double r0,double R,
		     double De);
void LanczosPI(matrix &tmatR,diagmat &IroR,diagmat &IR,
	       matrix &Hrot,diagmat &inertia,
	       matrix &Tphi,matrix &delcot,matrix &cotmat,diagmat &cosphi,
	       matrix &phiA,matrix &phiB,diagmat &pot4dsym,vector &v,
	       int Rpoints,int nsize1,int nsize2,int niter,double emin,
	       double emax);
matrix CMKE(int nsize,double dx, double mass);
vector CMgrid(int size,double length,double Rmin);
cvector transf(cvector &u,int size,matrix &Tx,matrix &Ty,int direction);
cvector g1d(vector &grid1d1,int nsize1d1,double width,double center1);
cvector symm(cvector &phixEB,cvector &chixEB,cvector &phiyEB,cvector &chiyEB,
	     int index,int size);
EXTERN double FORTRAN(ran1)(int *idum);
vector pdistribution(int ntrajectories,double mass,double beta,int &idnum);
vector xdistribution(int ntraject,double wx,double cx,double gx,double beta);
double U(double x,double wx,double cx,double gx);
void classicaltrajectory(vector &x,vector &v,vector &force,double mass,
			 double dt,double wx,double cx,double gx);
matrix pharmonic(int size);
matrix p2harmonic(int size);
diagmat harmonicoscillator(int size,double wx);
matrix xharmonic(int size);
matrix symmetrizedX(int size,int size2d,int sym);
complex computetrace(const cmatrix &A,int size);
cmatrix positioncorrelationfunction(matrix &rho,matrix &x,cmatrix &u);
cmatrix positioncorrelationfunction(diagmat &rho,matrix &x,cdiagmat &u);
cmatrix permute2d(cmatrix &rho2d,int size);
matrix permute2d(matrix &rho2d,int size);
double centroidforce(double qc,int P,double wx,double cx,double gx,double &acc,int &idum,double &V,double &rho);
double Action(int P,vector &q,double wx,double cx,double gx);
double F(double x,double wx,double cx,double gx);
vector xdistribution(int ntraject,Interp &interppot,double beta);
void classicaltrajectory(vector &x,vector &v,vector &force,double mass,
			 double dt,Interp &interpforce1,Interp &interpforce2);
vector initialnecklace(double qc,int P);
vector RandomMoveAndConstraint(vector &oldx,int p,int P,int &idum,double dx,double qc);
vector centroidforce2d(double qc1,double qc2,int P,double wx,double cx,double gx,
		       double &acc,int &idum,double &pot,char particle,double &rho);
double PermutedAction2d(int P,vector &q1,vector &q2,double wx,double cx,double gx);
double Action2d(int P,vector &q1,vector &q2,double wx,double cx,double gx);
EXTERN void FORTRAN(polin2)(double *x1a,double *x2a,double *ya,int *m,int *n,
				double *x1,double *x2,double *y,double *dy);
void classicaltrajectory(vector &x,vector &v,vector &force,double mass,
			 double dt,Interp2d &forceX,
			 Interp2d &forceY);
matrix xydistribution(int ntraject,Interp2d &interppot,double beta);
double computetrace(const matrix &A,int size);
matrix FM(int P);
cmatrix traceXP(int Nkmax,double kmin,double lmin,int size,
		double wx,double cx,double gx,
		matrix &T,vector &gridx,matrix &H,double beta,
		double dk,double dl,
		matrix &x,matrix &xp,cmatrix &traceforce,cmatrix &U,cmatrix &tracex);
void classicaltrajectory1d(vector &x,vector &v,vector &force,double mass,
			 double dt,Interp &interpforce1);
vector xdistribution(int ntraject,double wx,double cx,double gx,double beta,int &idum);
void classicaltrajectory1d(double &x,double &v,double &force,double mass,
			 double dt,Interp &interpforce1);
double window(double time,double width);
cvector traceHk(int Nkmax,double kmin,int size,double wx,double cx,double gx,
		matrix &T,vector &gridx,matrix &Hxtemp,double beta,double dk,
		matrix &x,cvector &traceforce);
cvector traceHkIm(int Nkmax,double kmin,int size,double wx,double cx,double gx,
		matrix &T,vector &gridx,matrix &Hxtemp,double beta,double dk,
		matrix &x,cvector &traceforce);
cvector shifttr(double x0,int Nkmax,cvector &tracevec,double dk,double kmin);
EXTERN void FORTRAN(arpack)(double *ddr,double *G11,double *G22,double *G33,double *G12,
			    double *G13,double *G23,double *VaplusV,int *size,
int *maxnev,int *maxncv,double *emin,double *emax,int *maxn,
			    double *v, double *workl, double *workd, double *d, double *resid,double * ax,int *select,int *niter,double *tolerance,int *labelsym,double *wp, int *ngood);
