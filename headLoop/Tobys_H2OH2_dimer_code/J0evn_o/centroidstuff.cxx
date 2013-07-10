#include "BF.h"
//static const unsigned long RANDMAX = 2147483647;
matrix xharmonic(int size)
{
  matrix xexact(size,size);
  for (int i=0;i<size;i++)
    for (int ip=0;ip<size;ip++) {
      xexact(i,ip)=0.;
      if (ip == (i+1))
	xexact(i,ip)+=sqrt((double)i+1.);
      if (ip == (i-1))
	xexact(i,ip)+=sqrt((double)i); 
      xexact(i,ip)*=(1./(sqrt(2.)));
    }
  return xexact;
}
matrix pharmonic(int size)
{
  // actually -ip
  matrix xexact(size,size);
  for (int i=0;i<size;i++)
    for (int ip=0;ip<size;ip++) {
      xexact(i,ip)=0.;
      if (i == (ip+1))
	xexact(i,ip)+=sqrt((double)ip+1.);
      if (i == (ip-1))
	xexact(i,ip)-=sqrt((double)ip); 
      xexact(i,ip)*=(1./(sqrt(2.)));
    }
  return xexact;
}
matrix p2harmonic(int size)
{
  // actually -ip^2
  matrix xexact(size,size);
  for (int i=0;i<size;i++)
    for (int ip=0;ip<size;ip++) {
      xexact(i,ip)=0.;
      if (ip == (i+2))
	xexact(i,ip)+=sqrt((double)i+1.)*sqrt((double)i+2.);
      if (ip == (i-2))
	xexact(i,ip)+=sqrt((double)i)*sqrt((double)i-1); 
      if (ip == i)
	xexact(i,ip)+=(-2.*(double)i-1.);
      xexact(i,ip)*=(1./((2.)));
    }
  return xexact;
}
diagmat harmonicoscillator(int size,double wx)
{
  diagmat H0(size);
  for (int i=0; i<size;i++)
    H0(i)=wx*((double)i+.5);
  return H0;
}
void classicaltrajectory(vector &x,vector &v,vector &force,double mass,
			 double dt,double wx,double cx,double gx)

  // velocity-Verlet  algorithm of Swope et al.
  // as described in Allen and Tildesley, p. 81
{
  x(0)=x(0)+dt*v(0)+.5*dt*dt*force(0)/mass;
  x(1)=x(1)+dt*v(1)+.5*dt*dt*force(1)/mass;
  // mid velocity
  v(0)=v(0)+.5*dt*force(0)/mass;
  v(1)=v(1)+.5*dt*force(1)/mass;
  // new force
  force(0)=F(x(0),wx,cx,gx);
  force(1)=F(x(1),wx,cx,gx);
  // new velocity
  v(0)=v(0)+.5*dt*force(0)/mass;
  v(1)=v(1)+.5*dt*force(1)/mass;
  return;
}
void classicaltrajectory(vector &x,vector &v,vector &force,double mass,
			 double dt,Interp &interpforce1,Interp &interpforce2)

  // velocity-Verlet  algorithm of Swope et al.
  // as described in Allen and Tildesley, p. 81
{
  x(0)=x(0)+dt*v(0)+.5*dt*dt*force(0)/mass;
  x(1)=x(1)+dt*v(1)+.5*dt*dt*force(1)/mass;
  // mid velocity
  v(0)=v(0)+.5*dt*force(0)/mass;
  v(1)=v(1)+.5*dt*force(1)/mass;
  // new force
  force(0)=interpforce1.interp(x(0));
  force(1)=interpforce2.interp(x(1));
  // new velocity
  v(0)=v(0)+.5*dt*force(0)/mass;
  v(1)=v(1)+.5*dt*force(1)/mass;
  return;
}
void classicaltrajectory1d(vector &x,vector &v,vector &force,double mass,
			 double dt,Interp &interpforce1)
{
  x(0)=x(0)+dt*v(0)+.5*dt*dt*force(0)/mass;
  // mid velocity
  v(0)=v(0)+.5*dt*force(0)/mass;
  // new force
  force(0)=interpforce1.interp(x(0));
  // new velocity
  v(0)=v(0)+.5*dt*force(0)/mass;
  return;
}
void classicaltrajectory1d(double &x,double &v,double &force,double mass,
			 double dt,Interp &interpforce1)
{
  x=x+dt*v+.5*dt*dt*force/mass;
  // mid velocity
  v=v+.5*dt*force/mass;
  // new force
  force=interpforce1.interp(x);
  // new velocity
  v=v+.5*dt*force/mass;
  return;
}
void classicaltrajectory(vector &x,vector &v,vector &force,double mass,
			 double dt,Interp2d &forceX,
			 Interp2d &forceY)

  // velocity-Verlet  algorithm of Swope et al.
  // as described in Allen and Tildesley, p. 81
{
  double dt2=.5*dt*dt;
  double dthalf=.5*dt;
  //x(0)=x(0)+dt*v(0)+.5*dt*dt*force(0)/mass;
  //x(1)=x(1)+dt*v(1)+.5*dt*dt*force(1)/mass;
  x(0)+=dt*v(0)+dt2*force(0);
  x(1)+=dt*v(1)+dt2*force(1);

  // mid velocity
  //v(0)=v(0)+.5*dt*force(0)/mass;
  //v(1)=v(1)+.5*dt*force(1)/mass;
  v(0)+=dthalf*force(0);
  v(1)+=dthalf*force(1);
  // new force
  force(0)=forceX.interp2d(x(0),x(1));
  force(1)=forceY.interp2d(x(0),x(1));
  // new velocity
  //v(0)=v(0)+.5*dt*force(0)/mass;
  //v(1)=v(1)+.5*dt*force(1)/mass;
  v(0)+=dthalf*force(0);
  v(1)+=dthalf*force(1);

  return;
}

// maby the separable potential is replaced by an effective potential 
// with some coupling in the case of F/B
vector xdistribution(int ntraject,double wx,double cx,double gx,double beta)
{
  double oldx,newx,oldV,newV,betaV,dV,betadV,ratio;
  int idum=-1;
  vector xdist(ntraject);
  int equilibrationsteps=5000;
  int Nmax=ntraject+equilibrationsteps;
  int countpoints=-1;
  double x=0.;
  double dx=.1;
  int Naccepted=0;
  for (int n=0;n<Nmax;n++) {    
    oldx=x;
    oldV=U(x,wx,cx,gx);
    newx=oldx;
    newx+=(2.*FORTRAN(ran1)(&idum)-1.)*dx;
    newV=U(newx,wx,cx,gx); 
    // not efficient because total energy is recalculated
    dV=newV-oldV;
    betadV=beta*dV;
    if (betadV <  75.) {
      if (betadV <= 0. ) {
	x=newx;
	Naccepted+=1;
      }
      else 
	if ((exp(-betadV)) > FORTRAN(ran1)(&idum) ) {
	  x=newx;
	  Naccepted+=1;
	}
    }
    if (n>equilibrationsteps) {
      countpoints+=1;
      xdist(countpoints)=x;
    }  
    if (n%50 == 0) {
      ratio=(double)Naccepted/50.;
      //if (ratio > .5) dx*=1.05; else dx*=.95;
      // Udo's trick
      if (ratio < .25) ratio=.25;
      dx=dx*pow(2.*ratio,2.);
      Naccepted=0;
    }

  }
  return xdist;
}
vector xdistribution(int ntraject,double wx,double cx,double gx,double beta,int &idum)
{
  double oldx,newx,oldV,newV,betaV,dV,betadV,ratio;
  vector xdist(ntraject);
  int equilibrationsteps=100;
  int Nmax=ntraject+equilibrationsteps;
  int countpoints=-1;
  double x=0.;
  double dx=.1;
  int Naccepted=0;
  for (int n=0;n<Nmax;n++) {    
    oldx=x;
    oldV=U(x,wx,cx,gx);
    newx=oldx;
    newx+=(2.*FORTRAN(ran1)(&idum)-1.)*dx;
    newV=U(newx,wx,cx,gx); 
    // not effcient because total energy is recalculated
    dV=newV-oldV;
    betadV=beta*dV;
    if (betadV <  75.) {
      if (betadV <= 0. ) {
	x=newx;
	Naccepted+=1;
      }
      else 
	if ((exp(-betadV)) > FORTRAN(ran1)(&idum) ) {
	  x=newx;
	  Naccepted+=1;
	}
    }
    if (n>=equilibrationsteps) {
      countpoints+=1;
      xdist(countpoints)=x;
    }  
    if (n%50 == 0) {
      ratio=(double)Naccepted/50.;
      if (ratio > .5) dx*=1.05; else dx*=.95;
      // Udo's trick
      //if (ratio < .25) ratio=.25;
      //dx=dx*pow(2.*ratio,2.);
      Naccepted=0;
    }

  }
  return xdist;
}
vector xdistribution(int ntraject,Interp &interppot,double beta)
{
  double oldx,newx,oldV,newV,betaV,dV,betadV,ratio;
  int idum=-1;
  vector xdist(ntraject);
  int equilibrationsteps=5000;
  int Nmax=ntraject+equilibrationsteps;
  int countpoints=-1;
  double x=0.;
  double dx=.1;
  int Naccepted=0;
  for (int n=0;n<Nmax;n++) {    
    oldx=x;
    oldV=interppot.interp(x);
    newx=oldx;
    newx+=(2.*FORTRAN(ran1)(&idum)-1.)*dx;
    newV=interppot.interp(newx); 
    // not effcient because total energy is recalculated
    dV=newV-oldV;
    betadV=beta*dV;
    if (betadV <  75.) {
      if (betadV <= 0. ) {
	x=newx;
	Naccepted+=1;
      }
      else 
	if ((exp(-betadV)) > FORTRAN(ran1)(&idum) ) {
	  x=newx;
	  Naccepted+=1;
	}
    }
    if (n>equilibrationsteps) {
      countpoints+=1;
      xdist(countpoints)=x;
    }  
    if (n%50 == 0) {
      ratio=(double)Naccepted/50.;
      //if (ratio > .5) dx*=1.05; else dx*=.95;
      // Udo's trick
      if (ratio < .25) ratio=.25;
      dx=dx*pow(2.*ratio,2.);
      Naccepted=0;
    }

  }
  return xdist;
}
matrix xydistribution(int ntraject,Interp2d &interppot,double beta)
{
  double oldx,newx,oldV,newV,betaV,dV,betadV,ratio;
  double oldy,newy;
  int idum=-1;
  matrix xydist(ntraject,2);
  int equilibrationsteps=100000;
  equilibrationsteps=0;
  int Nmax=ntraject+equilibrationsteps;
  int countpoints=-1;
  double x=0.;
  double y=0.;
  double dx=.01;
  int Naccepted=0;
  double xmin=-2.;
  double xmax=2.;
  for (int n=0;n<Nmax;n++) {    
      oldx=x;
      oldy=y;
      oldV=interppot.interp2d(x,y);
      
      newx=oldx; newy=oldy;
      newx+=(2.*FORTRAN(ran1)(&idum)-1.)*dx;
      newy+=(2.*FORTRAN(ran1)(&idum)-1.)*dx;
      
      newV=interppot.interp2d(newx,newy); 
      
      dV=newV-oldV;
      betadV=beta*dV;
      if (betadV <  75.) {
	  if (betadV <= 0. ) {
	      x=newx; y=newy;
	      Naccepted+=1;
	  }
	  else 
	      if ((exp(-betadV)) > FORTRAN(ran1)(&idum) ) {
		  x=newx; y=newy;
		  Naccepted+=1;
	      }
      }
      if (n>equilibrationsteps) {
	  countpoints+=1;
	  xydist(countpoints,0)=x;
	  xydist(countpoints,1)=y;
      }  
      if (n%50 == 0) {
	  ratio=(double)Naccepted/50.;
	  if (ratio > .5) dx*=1.05; else dx*=.95;
	  Naccepted=0;
      }
  }
  
  return xydist;
}
double U(double x,double wx,double cx,double gx)
{
    return  .5*wx*wx*pow(x,2.)+cx*pow(x,3.)+gx*pow(x,4.);
}
double F(double x,double wx,double cx,double gx)
{
  return -(wx*wx*x+cx*pow(x,2.)*3.+4.*gx*pow(x,3.));
}
vector pdistribution(int ntrajectories,double mass,double beta,int &idum)
{
  vector pdist(ntrajectories);
  for (int i=0;i<ntrajectories/2;i++) {
    double r1=FORTRAN(ran1)(&idum);
    double r2=FORTRAN(ran1)(&idum);
    double x1=sqrt(-2.*log(r1))*cos(2.*M_PI*r2);
    double x2=sqrt(-2.*log(r1))*sin(2.*M_PI*r2);
    pdist(2*i)=sqrt(1./(beta*mass))*x1;
    pdist(2*i+1)=sqrt(1./(beta*mass))*x2;
  }
  return pdist;
      
}
/*
double Action(int P,vector &q,double wx,double cx,double gx)
{
  double S=0.;
  double A=(MASS*(double)P/(2.*BETA));
  double B=BETA/(double)P;
  for (int i=0;i<P-1;i++)
    S+=A*pow(q(i+1)-q(i),2.)+B*U(q(i),wx,cx,gx);
  S+=A*pow(q(0)-q(P-1),2.)+B*U(q(P-1),wx,cx,gx);
  return S;
}
*/

vector initialnecklace(double qc,int P)
{
    int i;
  vector x(P);
  double q0;
  for (i=0;i<P;i++) x(i)=qc;
  // impose centroid constraint
  q0=0.; 
  for (i=0;i<P;i++) q0+=x(i);
  q0/=(double)P;
  double deltaq=q0-qc;
  for (i=0;i<P;i++) x(i)-=(deltaq);
  return x;
}
  
vector RandomMoveAndConstraint(vector &oldx,int p,int P,int &idum,double dx,double qc)
{
    int i;
  vector newx=oldx;
  // make a move
  //for (p=0;p<P;p++)
  newx(p)+=(2.*FORTRAN(ran1)(&idum)-1.)*dx;  
  // impose centroid constraint
  double q0=0.; 
  for (i=0;i<P;i++) q0+=newx(i);
  q0/=(double)P; 
  double deltaq=q0-qc;
  for (i=0;i<P;i++) newx(i)-=(deltaq);
  return newx;
}

matrix FM(int P)
{
    int k;
    matrix fm(P,P);
    for (k=0;k<P;k++) fm(k,0)=1.;
    for (k=0;k<P/2;k++) {
	fm(2*k,P-1)=1.;
	fm(2*k+1,P-1)=-1.;
    }
    double sqrtP=(sqrt((double)P/2.));
    for (k=0;k<P;k++){
	double phase=(2.*M_PI*((double)k))/(double)P;
	for (int p=0;p<P/2-1;p++) {
	    fm(k,2*p+1)=sqrt(2.)*cos(phase*((double)p+1.));
	    fm(k,2*p+2)=-sqrt(2.)*sin(phase*((double)p+1.));
	    //fm(k,2*p)/=sqrtP;
	    //fm(k,2*p+1)/=sqrtP;
	}
    }
    return fm;
}
cvector transf(cvector &u,int size,matrix &Tx,matrix &Ty,int direction)
{
  int i,j,ip,jp;
  cvector v(size*size);
  if (direction == 1) {
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	v(i*size+j)=complex(0.,0.);
	for (ip=0;ip<size;ip++) {
	  for (jp=0;jp<size;jp++) {
	    v(i*size+j)+=complex(Tx(ip,i)*Ty(jp,j),0.)*u(ip*size+jp);
	  }
	}
      }
    }
  }
  else {
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	v(i*size+j)=complex(0.,0.);
	for (ip=0;ip<size;ip++) {
	  for (jp=0;jp<size;jp++) {
	    v(i*size+j)+=complex(Tx(i,ip)*Ty(j,jp),0.)*u(ip*size+jp);
	  }
	}
      }
    }
  }
  return v;
}
cvector symm(cvector &phixEB,cvector &chixEB,cvector &phiyEB,cvector &chiyEB,
	     int index,int size)
{
  cvector v(size*size);
  for (int i=0;i<size;i++) {
    for (int j=0;j<size;j++) {
      if (index == 1) 
	v(i*size+j)=(phixEB(i)*chiyEB(j)+chixEB(i)*phiyEB(j));
      else
	if (index == -1) 
	  v(i*size+j)=(phixEB(i)*chiyEB(j)-chixEB(i)*phiyEB(j));
	else
	  v(i*size+j)=phixEB(i)*chiyEB(j);
    }
  }
  cnormalise(v);
  return v;
}
matrix permute2d(matrix &rho2d,int size)
{
  int size2d=size*size;
  matrix a(size2d,size2d);
  for (int i=0;i<size;i++)
    for (int j=0;j<size;j++)
      for (int ip=0;ip<size;ip++)
	for (int jp=0;jp<size;jp++) 
	  a(i*size+j,ip*size+jp)=rho2d(i*size+j,jp*size+ip);
  return a; 
}
cmatrix permute2d(cmatrix &rho2d,int size)
{
  int size2d=size*size;
  cmatrix a(size2d,size2d);
  for (int i=0;i<size;i++)
    for (int j=0;j<size;j++)
      for (int ip=0;ip<size;ip++)
	for (int jp=0;jp<size;jp++) 
	  a(i*size+j,ip*size+jp)=rho2d(i*size+j,jp*size+ip);
  return a; 
}
cmatrix positioncorrelationfunction(matrix &rho,matrix &x,cmatrix &u)
{
  cmatrix A=complexm(rho)*complexm(x)*transpose(u)*complexm(x)*u;
  return A;
}
cmatrix positioncorrelationfunction(diagmat &rho,matrix &x,cdiagmat &u)
{
  cmatrix A=complexm(rho)*complexm(x)*transpose(u)*complexm(x)*u;
  return A;
}
complex computetrace(const cmatrix &A,int size)
{
  complex tr=complex(0.,0.);
  for (int i=0;i<size;i++) tr+=A(i,i);
  return tr;
}
double computetrace(const matrix &A,int size)
{
  double tr=0.;
  for (int i=0;i<size;i++) tr+=A(i,i);
  return tr;
}
matrix symmetrizedX(int size,int size2d,int sym)
{
  int i,j,ip,jp;
  if (sym == -1) {
    int sizeF=size2d;
    matrix xexactF(sizeF,sizeF);
    int indexminus=-1;
    for (i=0;i<size;i++)
      for (j=0;j<i;j++){
	indexminus+=1;
	int indexminusp=-1;
	for (ip=0;ip<size;ip++)
	  for (jp=0;jp<ip;jp++){
	    indexminusp+=1;
	    xexactF(indexminus,indexminusp)=0.;
	    if (ip == (i+1)  && jp == j) {
	      xexactF(indexminus,indexminusp)+=sqrt((double)i+1.);
	    }
	    if (jp == (i+1)  && ip == j)
	      xexactF(indexminus,indexminusp)-=sqrt((double)i+1.);
	    if (ip == (j+1)  && jp == i)
	      xexactF(indexminus,indexminusp)-=sqrt((double)j+1.);
	    if (jp == (j+1)  && ip == i) {
	      xexactF(indexminus,indexminusp)+=sqrt((double)j+1.);
	    }
	    
	    if (ip == (i-1)  && jp == j)
	      xexactF(indexminus,indexminusp)+=sqrt((double)i);
	    if (jp == (i-1)  && ip == j)
	      xexactF(indexminus,indexminusp)-=sqrt((double)i);
	    if (ip == (j-1)  && jp == i)
	      xexactF(indexminus,indexminusp)-=sqrt((double)j);
	    if (jp == (j-1)  && ip == i)
	      xexactF(indexminus,indexminusp)+=sqrt((double)j);
	    
	    xexactF(indexminus,indexminusp)*=(1./(2.*sqrt(2.)));
	  }
      }
    return xexactF;
  }
  if (sym == 1) {
    int sizeB=size2d;
    matrix xexactB(sizeB,sizeB);
    int indexplus=-1;
    for (i=0;i<size;i++)
      for (j=0;j<=i;j++){
	indexplus+=1;
	int indexplusp=-1;
	for (ip=0;ip<size;ip++)
	  for (jp=0;jp<=ip;jp++){
	    indexplusp+=1;
	    xexactB(indexplus,indexplusp)=0.;
	    if (ip == (i+1)  && jp == j)
	      xexactB(indexplus,indexplusp)+=sqrt((double)i+1.);
	    if (jp == (i+1)  && ip == j)
	      xexactB(indexplus,indexplusp)+=sqrt((double)i+1.);
	    if (ip == (j+1)  && jp == i)
	      xexactB(indexplus,indexplusp)+=sqrt((double)j+1.);
	    if (jp == (j+1)  && ip == i)
	      xexactB(indexplus,indexplusp)+=sqrt((double)j+1.);
	    
	    if (ip == (i-1)  && jp == j)
	      xexactB(indexplus,indexplusp)+=sqrt((double)i);
	    if (jp == (i-1)  && ip == j)
	      xexactB(indexplus,indexplusp)+=sqrt((double)i);
	    if (ip == (j-1)  && jp == i)
	      xexactB(indexplus,indexplusp)+=sqrt((double)j);
	    if (jp == (j-1)  && ip == i)
	      xexactB(indexplus,indexplusp)+=sqrt((double)j);
	    
	    xexactB(indexplus,indexplusp)*=(1./(2.*sqrt(2.)));
	    if (ip==jp) xexactB(indexplus,indexplusp)*=.5;
	  }
	
      } 
    return xexactB;
  }
  
  if (sym == 0) {
    matrix xexact(size2d,size2d);
    for (i=0;i<size;i++)
      for (j=0;j<size;j++) {
	for (ip=0;ip<size;ip++)
	  for (jp=0;jp<size;jp++){
	    xexact(i*size+j,ip*size+jp)=0.;
	    if (ip == (i+1)  && jp == j)
	      xexact(i*size+j,ip*size+jp)+=sqrt((double)i+1.);
	    if (ip == (i-1)  && jp == j)
	      xexact(i*size+j,ip*size+jp)+=sqrt((double)i); 
	    xexact(i*size+j,ip*size+jp)*=(1./(sqrt(2.)));
	  }
      }
    return xexact;
  }
  cerr<<"wrong symmetry"<<endl;
  matrix xexactF(size2d,size2d);
  return xexactF;
}
cmatrix traceXP(int Nkmax,double kmin,double lmin,int size,
		double wx,double cx,double gx,
		matrix &T,vector &gridx,matrix &H,double beta,
		double dk,double dl,
		matrix &x,matrix &xp,cmatrix &traceforce,cmatrix &U)
{
  // Htilde = H+ i k x/beta + i l p/beta
    int indextr;
  cmatrix  tracemat(Nkmax,Nkmax);

  for (int k=0;k<Nkmax;k++) {
    double kval=kmin+(double)k*dk;
    for (int l=0;l<Nkmax;l++) {
      double lval=lmin+(double)l*dl;
      cmatrix Htilde=complexm(H+(lval/beta)*xp,(kval/beta)*x);
      cmatrix VL(size,size);
      cmatrix VR(size,size);
      cvector evaltr=diag(Htilde,VL,VR);

      VL=inverse(VR);
      cdiagmat zev(size);
      for (indextr=0;indextr<size;indextr++) 
	  zev(indextr)=exp(-complex(beta,0.)*evaltr(indextr));

      complex tr=complex(0.,0.);
      for (indextr=0;indextr<size;indextr++) 
	  tr+=zev(indextr);
      tracemat(k,l)=tr;

      diagmat Fx(size);
      for (indextr=0;indextr<size;indextr++) 
	  Fx(indextr)=F(gridx(indextr),wx,cx,gx);


      //Htilde=VR*zev*VL;
      Htilde=VR*zev*VL*transpose(U)*complexm((T)*Fx*transpose(T))*U;
      
       tr=complex(0.,0.);
      for (indextr=0;indextr<size;indextr++) 
	  tr+=Htilde(indextr,indextr);
          
      traceforce(k,l)=tr;
    }
  }
  
  return tracemat;
}
cvector traceHk(int Nkmax,double kmin,int size,double wx,double cx,double gx,
		matrix &T,vector &gridx,matrix &Hxtemp,double beta,double dk,
		matrix &x,cvector &traceforce)
{
  cvector tracevec(Nkmax);
  for (int k=0;k<Nkmax;k++) {
    double kval=kmin+(double)k*dk;
    cmatrix a=complexm(Hxtemp,(kval/beta)*x);
    cmatrix VL(size,size);
    cmatrix VR(size,size);
    cvector evaltr=diag(a,VL,VR);
    complex tr=complex(0.,0.);
    for (int indextr=0;indextr<size;indextr++) {
      tr+=exp(-complex(beta,0.)*evaltr(indextr));
    }
    tracevec(k)=tr;
    
    /* VL=inverse(VR);
       cdiagmat zev(size);
       for (indextr=0;indextr<size;indextr++) 
       zev(indextr)=exp(-complex(beta,0.)*evaltr(indextr));
       
       diagmat Fx(size);
       for (indextr=0;indextr<size;indextr++) 
       Fx(indextr)=F(gridx(indextr),wx,cx,gx);
       a=VR*zev*VL*complexm((T)*Fx*transpose(T));
       tr=complex(0.,0.);
       for (indextr=0;indextr<size;indextr++) {
       tr+=a(indextr,indextr);
       }
       traceforce(k)=tr;     */
    
    traceforce(k)=complex(0.,dk*(double)k/beta)*tr;
  }
  return tracevec;
}
cvector traceHkIm(int Nkmax,double kmin,int size,double wx,double cx,double gx,
		matrix &T,vector &gridx,matrix &Hxtemp,double beta,double dk,
		matrix &x,cvector &traceforce)
{
  cvector tracevec(Nkmax);
  for (int k=0;k<Nkmax;k++) {
    double kval=kmin+(double)k*dk;
    cmatrix a=complexm(Hxtemp-(kval/beta)*x);
    cmatrix VL(size,size);
    cmatrix VR(size,size);
    cvector evaltr=diag(a,VL,VR);
    complex tr=complex(0.,0.);
    for (int indextr=0;indextr<size;indextr++) {
      tr+=exp(-complex(beta,0.)*evaltr(indextr));
    }
    tracevec(k)=tr;    
    
    traceforce(k)=complex(0.,dk*(double)k/beta)*tr;
  }
  return tracevec;
}
cvector shifttr(double x0,int Nkmax,cvector &tracevec,double dk,double kmin)
{
    cvector gkx0(Nkmax);
    cvector tracevectemp(Nkmax);
    for (int k=0;k<Nkmax;k++) {
	double kval=kmin+(double)k*dk;
	for (int kp=0;kp<Nkmax;kp++) {
	    double kvalp=kmin+(double)kp*dk;
	    gkx0(kp)=exp(complex(-1.*pow((double)(kval-kvalp),2.),-(double)(kval-kvalp)*x0));
	}
	tracevectemp(k)=cconj(gkx0)*tracevec;
    }
    return tracevectemp;
}
double window(double time,double width)
{
    double g=sqrt(1./(2.*M_PI*width*width))*exp(-time*time/(2.*width*width));
    return g;
}
