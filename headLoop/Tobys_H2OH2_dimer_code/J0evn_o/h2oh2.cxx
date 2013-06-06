#include "peckeris.h"
#include "random.h"
static double betaGlobal;
void lanczosvectors(vector &alpha,vector &beta,vector &beta2,int niter,
					vector &eval,int ngood,matrix &evtr);
void EVanalysis(vector &grid,int size,int nconv,vector &ARv,double Ri,double Rf,
				int basistype,int size3d, diagmat &Rfunc,diagmat &rfunc,
				diagmat &R2func,diagmat &r2func,diagmat &sqrtweight);
void densityanalysis(vector &grid,int size,int nconv,vector &ARv,double Ri,double Rf,
					 int basistype,int size3d, diagmat &Rfunc,diagmat &rfunc,
					 diagmat &R2func,diagmat &r2func,diagmat &sqrtweight,
					 vector &eval);
double silvera(double rval);
double buck(double rval);
double buck1(double rval);
double buckpigs(double rval);
static const double amutoau=1./5.4857989586762187e-4; // amu to au
static const double au_to_amu = 5.4857989586762187e-4; // au to amu
static const double BohrToA=0.529177249; // Bohr to angstrom
static const double HtoW=2.19474631371017e5;//hatree to reciprical centimeter
static const double WavToMHz=29979.2458;
static const double hatokJmol=2625.;
static const double MofH2 = 2.015650642; // nist value
static const double MofD2 = 2.0141017780*2.;// nist
static const double MofHe4 = 4.0026032497; // nist mass of He in amu

//void vecanalysis(vector &ARvL,int ntotbs,int ngood,int nrpont,int numbas);
EXTERN void FORTRAN(eiganl)(double *ARvL,int *ntotbs,int *neigen,int *nrpont,int *numbas,int *icode,int *Kbgind,int* Jbig,int* ib000,int* ib100,int* ib101,int* ib110,int* ib111,int* jb111,int* ip,double *eval,int* niter,double* rgrid,double* h2_2mu,char *potfil,double* Tsin,double* Tcos,int *ncgrid,int *nlgrid,double* wgtgl,double* gcgrid,double* glgrid,int *maxfac,double* fact,int *jmax,int *ik,double *vmat,int* iplt3d);

EXTERN void FORTRAN(calcnb)(int *Jbig,int *jmax,int *ik,int *ip,int *numbas);

//EXTERN void  FORTRAN(prepare)(int *maxbas,int *maxj,int *idoubas,int *Kbigmx,int *maxfac,int *mxrgrd,int *nlgrid,
//			      int *ncgrid,double *Ah2o,double *Bh2o,double *Ch2o,double *h2_2mu,double *coroff,int *icorof,
//			      double *rotoff,int *irotof,double *diacor,double *diarot,double *rkin,int *icode,int *Kbgind,
//			      double *fact,double *wgtgl,double *glgrid,double *gcgrid,double *Tsin,double *Tcos,double *vmat,int *Jbig,
//			      int *jmax,int *ik,int *ip,double *rsmall,double *rlarge,int *nrgrid,double *rgrid,
//			      int *nrpont,int *nofcor,int *nofrot,int *nchi,int *nthe,int *nrad,double *vpes,double *gchi,
//			      double *gthe,double *grad,int *numbas,int *nKbig,bool *evnk);
EXTERN void FORTRAN(prepare)(int *numbas,int *icode,int *maxj,int *jmax,int *Jbig,int *ik,int *ip,int *Kbgind,
                        double *coriol,double *rotor,double *Ah2o,double *Bh2o,double *Ch2o,double *h2_2mu,int *ijcori,
                        int *ncorio,int *ijrot,int *nrotor,int *maxfac,double *fact,double *Tcos,double *Tsin,
                        int *nlgrid,int *ncgrid,double *wgtgl,double *glgrid,double *gcgrid,int *nchi,int *nthe,
                        int *nrad,double *vpes,double *gchi,double *gthe,double *grad,int *nrpont,double *rsmall,double *rlarge,
                        double *rgrid,int *nrgrid,double *vmat,double *rkin,int *nKbig,int *ib000,int *ib100,int *ib101,int *ib110,
                        int *ib111,int *jb111,char *potfil);

//EXTERN void  FORTRAN(hv)(double *vec,double *uec,int *mxrgrd,int *maxbas,double *vmat,int *ncgrid,int *nlgrid,double *rgrid,
//			 int *nrpont,int *numbas,int *Jbig,int *ip,int *Kbgind,int *Kbigmx,int *icode,
//			 double *Tcos,double *Tsin,int *nKbig,double *coroff,int *icorof,double *rotoff,int *irotof,
//			 double *diacor,double *diarot,int *nofcor,int *nofrot,double *rkin);
EXTERN void FORTRAN(hv)(double *vec,double *uec,double *vmat,int *ncgrid,int *nlgrid,double *rgrid,
                   int *nrpont,int *numbas,int *Jbig,int *ip,int *Kbgind,int *icode,double *Tcos,
                   double *Tsin,int *nKbig,double *coriol,int *ijcori,double *rotor,int *ijrot,
                   int *ncorio,int *nrotor,double *rkin);

main(int argc,char **argv)
{
  int i,j,k,row,n;
  
  // read the input
  if (argc != 14) {
    cerr<<"usage: "<<argv[0]<<" Jbig, jmax, ik: 0 (1) for even (odd), ip: 0 (1) for even (odd) parity, rsmall, rlarge,  nrgrid, Vmax, niter, pot_file, nlgrid, ncgrid iplt3d "<<endl;
    exit(0);
  }
 /*
      write(6,*)'punch in Jbig, jmax, 0 (1) for even (odd) k',
     +          ' and 0 (1) for even (odd) parity'
      read(5,*)Jbig, jmax, ik,ip
      write(6,*)Jbig,jmax,ik,ip
      if(Jbig.ge.jmax) stop 'Jbig has to be smaller than jmax'
      if(jmax.gt.maxj) stop 'jmax exceeds maxj'
      write(6,*)'punch in rsmall, rlarge, and nrgrid'
      read(5,*)rsmall,rlarge,nrgrid
  */


  int Jbig=atoi(argv[1]);
  int jmax=atoi(argv[2]);
  int ik=atoi(argv[3]);
  int ip=atoi(argv[4]);
  double rsmall=atof(argv[5]);
  double rlarge=atof(argv[6]);
  int nrgrid=atoi(argv[7]);
  double  Vmax=atof(argv[8]);
  int niter=atoi(argv[9]);
  char* potfil=argv[10];
  int nlgrid=atoi(argv[11]);
  int ncgrid=atoi(argv[12]);
  int iplt3d=atoi(argv[13]);

  if(iplt3d != 0 && iplt3d != 1)
  {
    cerr<<"iplt3d can be only 0 or 1. it has the current value of "<<iplt3d<<endl;
    exit(0);
  }

  for (i=0;i<argc;i++)
    cout<<argv[i]<<" ";
  cout<<endl;
  cout.flush();

  // ... nrgrid is the number of intervals=
  // ... nrpont is the number of grid point, in addition to the two boundaries
  int nrpont=nrgrid-1;
  //   if(nrpont.gt.mxrgrd)stop' too many radial grid points'


  int maxbas=1000,maxj=30,idoubas=2000,mxrgrd=100;
  // ... the value of Kbig is restricted by Jbig and can't be too large
  // ... say the max of Kbig=10
  int Kbigmx=10,maxfac=50;
  // ... nlgrid=13 is the choise of AvdA
  // ... ncgrid=10 is gives unexpected null matrix elements for jmax=12
  // ... temporarily choose ncgrid=15
//int nlgrid=13,ncgrid=13;
  //int nlgrid=26,ncgrid=26;
	// ... the following are the rotational constants of H2O,
	// ... in the unit of cm-1.
  // AvdA constants:
     double Ah2o=27.8806,Bh2o=14.5216,Ch2o=9.2778;
  // Wang and Carrington constants for v2=0
//double Ah2o=27.8572,Bh2o=14.5145,Ch2o=9.2799;
// shrink all rotational constants by 10 times to have a slower rotor
//double Ah2o=2.78572,Bh2o=1.45145,Ch2o=.92799;
// shrink all rotational constants by 100 times to have a slower rotor
//double Ah2o=0.278572,Bh2o=0.145145,Ch2o=0.092799;
//double Aso2=2.027354,Bso2=0.344170,Cso2=0.29353;
  // Wang and Carrington constants for v2=1
  //double Ah2o=31.0847,Bh2o=14.6748,Ch2o=9.1361;
  // ... the following are the rotational constants of H2O,
  // ... in the unit of Hartree (atomic unit)
  // double Ah2o=1.27033e-04,Bh2o=6.61653e-05,Ch2o=4.2273e-05;
  // ... another important energy factor is hbar^2/2mu in the atomic unit
   double h2_2mu=33.20859; // h2 + h2o
 // double h2_2mu=29.86612; // h2 only
//double h2_2mu=30.8073; // h2 + so2 with 32-S and 16-O
  // ... grids for potential energy surface
  int nchi=10,nthe=19,nrad=48;

  int numbas=0;
  bool evnk=true;

  FORTRAN(calcnb)(&Jbig,&jmax,&ik,&ip,&numbas);

  double* coriol=new double[numbas*4];
  int* ijcori=new int[numbas*4*2];
  double* rotor=new double[numbas*4];
  int* ijrot=new int[numbas*4*2];
  double* rkin=new double[nrpont*(nrpont+1)/2];
  int* icode=new int[numbas];
  int* Kbgind=new int[Jbig+2];
  double* fact=new double[1+maxfac];
  double* wgtgl=new double[nlgrid];
  double* glgrid=new double[nlgrid];
  double* gcgrid=new double[ncgrid];
  double* Tsin=new double[nlgrid*ncgrid*numbas];
  double* Tcos=new double[nlgrid*ncgrid*numbas];
  double*   vmat=new double[nrpont*nlgrid*ncgrid];
  double* rgrid=new double[nrpont];
  double*   vpes=new double[(nchi+10)*(nthe+10)*nrad];
  double* gchi=new double[nchi+10];
  double*  gthe=new double[nthe+10];
  double* grad=new double[nrad];
  double* vec2=new double[nrpont*numbas];
  double*  uec2=new double[nrpont*numbas];

   // integer and double modified by prepare
  int ncorio=0;
  int nrotor=0;
//int numbas=0;
  int nKbig=0;
  int ib000=0;
  int ib100=0;
  int ib101=0;
  int ib110=0;
  int ib111=0;
// jb111 is actually the basis index for j=1, K=-1, K=1
  int jb111=0;



//FORTRAN(prepare)(&maxbas,&maxj,&idoubas,&Kbigmx,&maxfac,&mxrgrd,&nlgrid,
//  &ncgrid,&Ah2o,&Bh2o,&Ch2o,&h2_2mu,coroff,icorof,
//  rotoff,irotof,diacor,diarot,rkin,icode,Kbgind,
//  fact,wgtgl,glgrid,gcgrid,Tsin,Tcos,vmat,&Jbig,
//  &jmax,&ik,&ip,&rsmall,&rlarge,&nrgrid,rgrid,
//  &nrpont,&nofcor,&nofrot,&nchi,&nthe,&nrad,vpes,gchi,
//  gthe,grad,&numbas,&nKbig,&evnk);
  FORTRAN(prepare)(&numbas,icode,&maxj,&jmax,&Jbig,&ik,&ip,Kbgind,
                        coriol,rotor,&Ah2o,&Bh2o,&Ch2o,&h2_2mu,ijcori,
                        &ncorio,ijrot,&nrotor,&maxfac,fact,Tcos,Tsin,
                        &nlgrid,&ncgrid,wgtgl,glgrid,gcgrid,&nchi,&nthe,
                        &nrad,vpes,gchi,gthe,grad,&nrpont,&rsmall,&rlarge,
                        rgrid,&nrgrid,vmat,rkin,&nKbig,&ib000,&ib100,&ib101,
                        &ib110,&ib111,&jb111,potfil);


  // for (i=0;i<mxrgrd*nlgrid*ncgrid;i++)
  //  cout<<i<<" "<<vmat[i]<<endl;





  double* vec=new double[nrpont*numbas];
  double* uec=new double[nrpont*numbas];

  //test Hv


  int ntotbs=nrpont*numbas;

//for (int ib=0;ib<ntotbs;ib++) {
//  vec[ib]=1./sqrt((double)(ntotbs));
//  uec[ib]=0.;
//}


//  FORTRAN(hv)(vec,uec,vmat,&ncgrid,&nlgrid,rgrid,
//                 &nrpont,&numbas,&Jbig,&ip,Kbgind,icode,Tcos,
//                 Tsin,&nKbig,coriol,ijcori,rotor,ijrot,
//                 &ncorio,&nrotor,rkin);

  
  //ouput the result of Hv

//for (int ib=0;ib<ntotbs;ib++) 
//  cout<<"hv"<<ib+1<<" "<<uec[ib]<<" "<<vec[ib]<<endl;
//cout.flush();

//  exit(0);
  

  double emin=-100.;
  cout<<"emin= "<<emin<<endl;
  double emax=200.;
  //double emax=-emin;
  
  Rand *randomseed = new Rand(1);
  
  int ngood;
  
  vector evalerr(niter);
  vector eval(niter);
  vector alpha(niter);
  vector beta(niter+1);
  vector beta2(niter+1);
  
  double* rvec=new double[nrpont*numbas];  
  for (i=0;i<(nrpont*numbas);i++) rvec[i]=0.;

  for (int ib=0;ib<ntotbs;ib++) {
    uec[ib]=0.;
    vec[ib]=1./sqrt((double)ntotbs);
    rvec[ib]=0.;
  }
  //vec[0]=1.;

  cout<<"start iterations"<<endl;
  for (j=1;j<=niter;j++) {    

    for (int ib=0;ib<ntotbs;ib++) 
      uec[ib]=0.;

//    FORTRAN(hv)(vec,uec,&mxrgrd,&maxbas,vmat,&ncgrid,&nlgrid,rgrid,
//	      &nrpont,&numbas,&Jbig,&ip,Kbgind,&Kbigmx,icode,
//	      Tcos,Tsin,&nKbig,coroff,icorof,rotoff,irotof,
//	      diacor,diarot,&nofcor,&nofrot,rkin);
    FORTRAN(hv)(vec,uec,vmat,&ncgrid,&nlgrid,rgrid,
                   &nrpont,&numbas,&Jbig,&ip,Kbgind,icode,Tcos,
                   Tsin,&nKbig,coriol,ijcori,rotor,ijrot,
                   &ncorio,&nrotor,rkin);

    for (int ib=0;ib<ntotbs;ib++) 
      rvec[ib]+=uec[ib];
    
    alpha(j-1)=0.;
    for (int ib=0;ib<ntotbs;ib++) 
      alpha(j-1)+=vec[ib]*rvec[ib];

    for (int ib=0;ib<ntotbs;ib++) 
      rvec[ib]-=(alpha(j-1)*vec[ib]);
    
    beta2(j)=0.;
    for (int ib=0;ib<ntotbs;ib++) 
      beta2(j)+=rvec[ib]*rvec[ib];

    beta(j)=sqrt(beta2(j));

    for (int ib=0;ib<ntotbs;ib++) 
      rvec[ib]=(1./beta(j))*rvec[ib]; // to get v
    for (int ib=0;ib<ntotbs;ib++) 
      vec[ib]=(-beta(j))*vec[ib]; // prepare r check minus sign!!!
    
    for (int ib=0;ib<ntotbs;ib++) {
      uec[ib]=vec[ib];     // swapping
      vec[ib]=rvec[ib];
      rvec[ib]=uec[ib];
    }
    cout<<j<<" "<<alpha(j-1)<<" "<<beta(j)<<endl;
    if (j%100 == 0)
      cout<<"iteration "<<j<<endl;
  }                  
  lancbis(niter,eval,evalerr,emin,emax,ngood,alpha,beta,beta2);
  cout<<" ngood = "<<ngood<<endl;
  cout<<"E0= "<<eval(0)<<endl;
  // lanczos report:
  ofstream lancout("boundstates.out");
  ofstream lanczpeout("states_zpe.out");
  for (i=0;i<ngood;i++) {
    lancout<<eval(i)<<" "<<evalerr(i)<<endl;
    lanczpeout<<(eval(i)-eval(0))<<endl;
    }
  lancout.flush();
  lancout.close();
 
  // without eigenvectors
  //exit(0);

  // eigevectors

  if (ngood>5) ngood=5;
    
  matrix evtr(niter,ngood);
  lanczosvectors(alpha,beta,beta2,niter,eval,ngood,evtr);
    
  
  vector ARvL(ntotbs*ngood);
		
  // reset vectors


  for (int ib=0;ib<ntotbs;ib++) {
    uec[ib]=0.;
    vec[ib]=1./sqrt((double)ntotbs);
    rvec[ib]=0.;
  }

  for (j=0;j<niter;j++) {
    alpha(j)=0.;
    beta(j)=0.;
    beta2(j)=0.;

  }
  beta(niter)=0.;
  beta2(niter)=0.;



  ofstream lvecout("lv");
  vector cumulnorm(ngood);
  for (j=1;j<=niter;j++) {	
    lvecout<<j<<" ";
    for (n=0;n<ngood;n++) {
      double coeff2=pow(evtr(j-1,n),2.);
      lvecout<<coeff2<<" ";
      cumulnorm(n)+=coeff2;
      lvecout<<cumulnorm(n)<<" ";
    }
    lvecout<<endl;
  }
  

	
  // lanczos vector coefficent matrix
  for (n=0;n<ngood;n++) cumulnorm(n)=0.;

  for (j=1;j<=niter;j++) {	
		
    // tranform vector	
    for (n=0;n<ngood;n++) {
      double treshold=pow(evtr(j-1,n),2.);
      cumulnorm(n)+=treshold;
      for (row=0;row<ntotbs;row++){
	double coeff=evtr(j-1,n)*vec[row];	  
	if (cumulnorm(n) <(1.-1.e-16))
	  ARvL(row+ntotbs*n)+=coeff;	
      }
    }
		
    for (int ib=0;ib<ntotbs;ib++) 
      uec[ib]=0.;

//    FORTRAN(hv)(vec,uec,&mxrgrd,&maxbas,vmat,&ncgrid,&nlgrid,rgrid,
//	      &nrpont,&numbas,&Jbig,&ip,Kbgind,&Kbigmx,icode,
//	      Tcos,Tsin,&nKbig,coroff,icorof,rotoff,irotof,
//	      diacor,diarot,&nofcor,&nofrot,rkin);
    FORTRAN(hv)(vec,uec,vmat,&ncgrid,&nlgrid,rgrid,
                   &nrpont,&numbas,&Jbig,&ip,Kbgind,icode,Tcos,
                   Tsin,&nKbig,coriol,ijcori,rotor,ijrot,
                   &ncorio,&nrotor,rkin);

    for (int ib=0;ib<ntotbs;ib++) 
      rvec[ib]+=uec[ib];
    
    alpha(j-1)=0.;
    for (int ib=0;ib<ntotbs;ib++) 
      alpha(j-1)+=vec[ib]*rvec[ib];

    for (int ib=0;ib<ntotbs;ib++) 
      rvec[ib]-=(alpha(j-1)*vec[ib]);
    
    beta2(j)=0.;
    for (int ib=0;ib<ntotbs;ib++) 
      beta2(j)+=rvec[ib]*rvec[ib];

    beta(j)=sqrt(beta2(j));

    for (int ib=0;ib<ntotbs;ib++) 
      rvec[ib]=(1./beta(j))*rvec[ib]; // to get v
    for (int ib=0;ib<ntotbs;ib++) 
      vec[ib]=(-beta(j))*vec[ib]; // prepare r check minus sign!!!
    
    for (int ib=0;ib<ntotbs;ib++) {
      uec[ib]=vec[ib];     // swapping
      vec[ib]=rvec[ib];
      rvec[ib]=uec[ib];
    }
    cout<<j<<" "<<alpha(j-1)<<" "<<beta(j)<<endl;
    if (j%100 == 0)
      cout<<"iteration "<<j<<endl;
  }    

  
  // eigenvectors for good eigenvalues are in ARvL(row+ntotbs*n)

  // for the groudn state, eg, ARvL(row+ntotbs*0) is the vector in the original basis
  // calculate radial distribution function
  ofstream psi2rout("psi2r");
  ofstream psirout("psir");

  double rstep=(rlarge - rsmall)/nrgrid;
  rstep=rstep*0.529177249;

  for (i=0;i<nrpont;i++) { //R loop
    double rvalue=rgrid[i];
//  psi2rout<<rvalue<<" ";
    psi2rout<<rvalue*0.529177249<<" ";
    psirout<<rvalue<<" ";
    for (n=0;n<ngood;n++) {
      double psi2r=0.;
      double psir=0.;
      for (int ib=0;ib<numbas;ib++) {
	int row=i*numbas+ib;
	psi2r+=pow(ARvL(row+ntotbs*n),2.);
	psir+=ARvL(row+ntotbs*n);

      }
//    psi2rout<<psi2r<<" ";
      psi2rout<<psi2r/rstep<<" ";
      psirout<<psir<<" ";
    }
    psi2rout<<endl;
    psirout<<endl;
  }

  vecanalysis(ARvL,ntotbs,ngood,nrpont,numbas,icode,Kbgind,Jbig,ib000,ib100,ib101,ib110,ib111,jb111,ip,eval,niter,rgrid,h2_2mu,potfil,Tsin,Tcos,ncgrid,nlgrid,wgtgl,gcgrid,glgrid,maxfac,fact,jmax,ik,vmat,iplt3d);
  //cout<<ARvL(0)<<" "<<ARvL(ntotbs*ngood-1)<<" "<<icode[0]<<" "<<icode[numbas-1]<<endl;
	  
  cout<<"end of program"<<endl;
  

}
void lanczosvectors(vector &alpha,vector &beta,vector &beta2,int niter,
					vector &eval,int ngood,matrix &evtr)
{
    // copy stuff
    int i,j,ndis;
    double* lalpha=new double[niter];
    double* lbeta=new double[niter+1];
    double* lbeta2=new double[niter+1];
    double* leval=new double[ngood];
    lbeta[0]=0.;
    lbeta2[0]=0.;
    for (j=1;j<=niter;j++) {
		lalpha[j-1]=alpha(j-1);
		lbeta[j]=beta(j);
		lbeta2[j]=beta2(j);
    }
    for (j=0;j<ngood;j++) leval[j]=eval(j);
	
    double* wrk1=new double[niter];
    double* wrk2=new double[niter];
    double* levtr=new double[niter*ngood];
    double* mamin=new double[ngood];
	
    FORTRAN(trivec)(lalpha,lbeta,lbeta2,wrk1,wrk2,&niter,leval,&ngood,levtr,mamin);
	
    for (i=0;i<niter;i++)
		for (j=0;j<ngood;j++) 
			evtr(i,j)=levtr[i+j*niter];
    return;
}
void vecanalysis(vector &ARvL,int ntotbs,int ngood,int nrpont,int numbas,int* icode,int* Kbgind,int Jbig,int ib000,int ib100,int ib101,int ib110,int ib111,int jb111,int ip,vector &eval,int niter,double* rgrid,double h2_2mu,char* potfil,double* Tsin,double* Tcos,int ncgrid,int nlgrid,double* wgtgl,double* gcgrid,double* glgrid,int maxfac,double* fact,int jmax,int ik, double* vmat, int iplt3d){
	FORTRAN(eiganl)(ARvL.TheVector,&ntotbs,&ngood,&nrpont,&numbas,icode,Kbgind,&Jbig,&ib000,&ib100,&ib101,&ib110,&ib111,&jb111,&ip,eval.TheVector,&niter,rgrid,&h2_2mu,potfil,Tsin,Tcos,&ncgrid,&nlgrid,wgtgl,gcgrid,glgrid,&maxfac,fact,&jmax,&ik,vmat,&iplt3d);
 return;
}
