#ifndef LANCZOS_H
#define	LANCZOS_H

#include "vectClass.h" //This is to select which linear algebra class should be used.

void lancbis(int niter,VECT &eval,VECT &evalerr,double elmin,
			 double elmax,int &ngood,const VECT& alpha,const VECT& beta,
			 const VECT& beta2);

EXTERN void FORTRAN(trivec)(double *lalpha,double *lbeta,double *lbeta2,
							double *wrk1,double *wrk2,int *niter,
							double *eval,int *ngood,double *evtr,double *mamin);

EXTERN void FORTRAN(bisec)(double *lalpha,double *lbeta,double *lbeta2,
						   int *,double *wrk4,int *mp2,double *wrk2,int *,
						   double *elmin,double *elmax,int *ndis);

EXTERN void FORTRAN(isoev)(double *vs,int *mp,int *ndis);

EXTERN void FORTRAN(inverr)(double *lalpha,double *lbeta,double *wrk1,
							double *wrk2,double *wrk3,int *niter,double *vs,int *mp,double *err,int
							*ndis);

#endif