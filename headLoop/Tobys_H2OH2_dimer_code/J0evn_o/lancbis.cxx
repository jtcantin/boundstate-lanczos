#include "BF.h"
void lancbis(int niter,vector &eval,vector &evalerr,double elmin,
	     double elmax,int &ngood,const vector& alpha,const vector& beta,
	     const vector& beta2)
{
  // copy stuff
  int i,j,ndis;
  double* lalpha=new double[niter];
  double* lalphaql=new double[niter];
  double* lbeta=new double[niter+1];
  double* lbeta2=new double[niter+1];
  lbeta[0]=beta(0);
  lbeta2[0]=beta2(0);
  for (j=1;j<=niter;j++) {
    lalpha[j-1]=alpha(j-1);
    lbeta[j]=beta(j);
    lbeta2[j]=beta2(j);
  }
  int* mp2=new int[niter];
  double* wrk2=new double[2*niter];
  FORTRAN(bisec)(lalpha,lbeta,lbeta2,&niter,lalphaql,mp2,wrk2,&niter,
	 &elmin,&elmax,&ndis);
  delete[] wrk2;   

  double* vs=new double[ndis];
  double* vsp=new double[ndis];
  int* mp=new int[ndis];
  int* mpp=new int[ndis];
  
  for (i=0;i<ndis;i++) {
    vs[i]=lalphaql[i];
    mp[i]=mp2[i];
  }	
  FORTRAN(isoev)(vs,mp,&ndis);
  /*calculation of error */
  double* wrk1=new double[niter];
  wrk2=new double[niter+1];
  double* wrk3=new double[niter];
  double* err=new double[niter];
  double* errp=new double[niter];
  FORTRAN(inverr)(lalpha,lbeta,wrk1,wrk2,wrk3,&niter,vs,mp,err,&ndis);
  delete[] wrk2;
  delete[] wrk3;
  /* compress vs according to mp */
  int ii=0;
  for (i=0;i<ndis;i++) {
    if (mp[i] != 0 ) {
      vsp[ii]=vs[i];
      mpp[ii]=mp[i];
      errp[ii]=err[i];
      ii+=1;
    }
  }
  delete[] vs;	
  delete[] err;
  for (i=0;i<ii;i++) {
    mp[i]=mpp[i];
    eval(i)=vsp[i];
    evalerr(i)=errp[i];
  }
  delete[] mp;
  delete[] mpp;
  delete[] vsp;
  delete[] errp;
  ngood=ii;
}
