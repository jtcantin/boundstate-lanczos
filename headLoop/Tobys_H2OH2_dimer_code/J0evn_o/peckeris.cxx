#include "peckeris.h"
double delta(int i,int ip)
{
  if (i == ip) return 1.;
  else return 0.;
}
vector Hpsi(matrix &ddr2,diagmat &Rinv,diagmat &R2inv,diagmat &R4inv2,matrix &Rinv4ij,matrix &ddr,
	    diagmat &extraV,diagmat &V,int size,vector &v,double mass)
{
  int i,j,k,l,kp;
  vector u(size*size*size);
  vector work1(size*size*size);
  vector work2(size*size*size);

  // operate with ti term
  
  for (l=0;l<3;l++) {
    // operate with d/dr_k (matrix term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {
	  work1((i*size+j)*size+k)=0.;
	  for (kp=0;kp<size;kp++) {
	    if (l==0)
	      work1((i*size+j)*size+k)+=ddr(k,kp)*v((i*size+j)*size+kp);
	    if (l==1)
	      work1((i*size+j)*size+k)+=ddr(i,kp)*v((kp*size+j)*size+k);
	    if (l==2)
	      work1((i*size+j)*size+k)+=ddr(j,kp)*v((i*size+kp)*size+k);	    
	  }	 
	}
      }
    }

    // operate with d/dr_j  (matrix term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) { 
	  work2((i*size+j)*size+k)=0.;
	  for (kp=0;kp<size;kp++) {
	    if (l==0) 
	      work2((i*size+j)*size+k)+=ddr(j,kp)*work1((i*size+kp)*size+k);
	    if (l==1) 
	      work2((i*size+j)*size+k)+=ddr(k,kp)*work1((i*size+k)*size+kp);
	    if (l==2) 
	      work2((i*size+j)*size+k)+=ddr(i,kp)*work1((kp*size+k)*size+k);	    
	  }
	}
      }
    }

    // operate with 1/2Rj (diagonal term) and reuse d/drk
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) { 
	  if (l==0) 
	    work2((i*size+j)*size+k)-=R2inv(j)*work1((i*size+j)*size+k);
	  if (l==1) 
	    work2((i*size+j)*size+k)-=R2inv(k)*work1((i*size+j)*size+k);
	  if (l==2) 
	    work2((i*size+j)*size+k)-=R2inv(i)*work1((i*size+j)*size+k);	  
	}
      }
    }

    // operate with d/dr_j (matrix term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) { 
	  work1((i*size+j)*size+k)=0.;
	  for (kp=0;kp<size;kp++) {
	    if (l==0) 
	      work1((i*size+j)*size+k)+=ddr(j,kp)*v((i*size+kp)*size+k);
	    if (l==1) 
	      work1((i*size+j)*size+k)+=ddr(k,kp)*v((i*size+k)*size+kp);
	    if (l==2) 
	      work1((i*size+j)*size+k)+=ddr(i,kp)*v((kp*size+j)*size+k);	    
	  }	  
	}
      }
     }
    // operate with 1/2Rk (diagonal term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	 for (k=0;k<size;k++) { 
	   if (l==0)
	     work2((i*size+j)*size+k)-=R2inv(k)*work1((i*size+j)*size+k);
	   if (l==1)
	     work2((i*size+j)*size+k)-=R2inv(i)*work1((i*size+j)*size+k);
	   if (l==2)
	     work2((i*size+j)*size+k)-=R2inv(j)*work1((i*size+j)*size+k);
	 }
      }
    }
    // operate with 1/4RjRk (diagonal term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) { 
	  if (l==0)
	    work2((i*size+j)*size+k)+=Rinv4ij(j,k)*v((i*size+j)*size+k);
	  if (l==1)
	    work2((i*size+j)*size+k)+=Rinv4ij(k,i)*v((i*size+j)*size+k);
	  if (l==2)
	    work2((i*size+j)*size+k)+=Rinv4ij(i,j)*v((i*size+j)*size+k);
	}
      }
    }
    // operate with extraV (diagonal term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {
	  if (l==0)
	    work2((i*size+j)*size+k)*=extraV((i*size+j)*size+k);
	  if (l==1)
	    work2((i*size+j)*size+k)*=extraV((j*size+k)*size+i);
	  if (l==2)
	    work2((i*size+j)*size+k)*=extraV((k*size+i)*size+j);	  
	}
      }
    }
    // operate with 1/4R2 (diagonal term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {
	  if (l==0)
	    work2((i*size+j)*size+k)-=R4inv2(i)*v((i*size+j)*size+k);
	  if (l==1)
	    work2((i*size+j)*size+k)-=R4inv2(j)*v((i*size+j)*size+k);
	  if (l==2)
	    work2((i*size+j)*size+k)-=R4inv2(k)*v((i*size+j)*size+k);
	}
      }
    }
    // operate with operate with d/dr_i (matrix term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {
	  work1((i*size+j)*size+k)=0.;
	  for (kp=0;kp<size;kp++) {
	    if (l==0)
	      work1((i*size+j)*size+k)+=ddr(i,kp)*v((kp*size+j)*size+k);
	    if (l==1)
	      work1((i*size+j)*size+k)+=ddr(j,kp)*v((i*size+kp)*size+k);
	    if (l==2)
	      work1((i*size+j)*size+k)+=ddr(k,kp)*v((i*size+j)*size+kp);
	  }
	}
      }
    }
    // operate with 1/1Ri
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {
	  if (l==0)
	    work2((i*size+j)*size+k)+=Rinv(i)*work1((i*size+j)*size+k);
	  if (l==1)
	    work2((i*size+j)*size+k)+=Rinv(j)*work1((i*size+j)*size+k);
	  if (l==2)
	    work2((i*size+j)*size+k)+=Rinv(k)*work1((i*size+j)*size+k);	  
	}
      }
    }
    // operate with operate with d2/dr_idr_i (matrix term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {
	  work1((i*size+j)*size+k)=0.;
	  for (kp=0;kp<size;kp++) {
	    if (l==0)
	      work1((i*size+j)*size+k)+=ddr2(i,kp)*v((kp*size+j)*size+k);
	    if (l==1)
	      work1((i*size+j)*size+k)+=ddr2(j,kp)*v((i*size+kp)*size+k);
	    if (l==2)
	      work1((i*size+j)*size+k)+=ddr2(k,kp)*v((i*size+j)*size+kp);
	  }
	}
      }
    }
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {
	    work2((i*size+j)*size+k)+=work1((i*size+j)*size+k);	  
	}
      }
    }
    
    //multiply by -hbar^2/m and storage into u
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {
	  work2((i*size+j)*size+k)/=(-mass);
	  u((i*size+j)*size+k)+=work2((i*size+j)*size+k);
	}
      }
    }  
  }
  u=u+V*v;
  return u;
}

vector HpsiBLAS(diagmat &Rinv,diagmat &R2inv,diagmat &R4inv2,
		matrix &Rinv4ij,matrix &ddr,matrix &ddr2invrddr,
		diagmat &extraV,diagmat &V,int size,vector &v,double mass)
{
  int i,j,k,l;
  vector u(size*size*size);

  matrix Umat(size*size,size);
  matrix Vmat(size*size,size);
  matrix W1(size*size,size);
  matrix W1p(size*size,size);
  matrix W2(size*size,size);

  // operate with ti term
  
  for (l=0;l<3;l++) {
    // operate with d/dr_k (matrix term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {	  
	  // prepare for blas for matrix vector product
	  // fill the V matrix
	  if (l==0) // k sum
	    Vmat(i*size+j,k)=v((i*size+j)*size+k);
	  if (l==1) // i sum
	    Vmat(j*size+k,i)=v((i*size+j)*size+k);
	  if (l==2) // j sum
	    Vmat(i*size+k,j)=v((i*size+j)*size+k);
	}
      }
    }
    // blas for matrix vector product
    W1=Vmat*transpose(ddr); // i,j,k order for l=0

    // operate with d/dr_j  (matrix term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) { 	  
	  // do some reordering
	  if (l==0) // j sum
	    W1p(i*size+k,j)=W1(i*size+j,k);
	  if (l==1) // k sum
	    W1p(i*size+j,k)=W1(j*size+k,i);
	  if (l==2) // i sum
	    W1p(j*size+k,i)=W1(i*size+k,j);
	}
      }
    }
    // blas operation
    W2=W1p*transpose(ddr); // i,k,j order for l=0
    

    // operate with 1/2Rj (diagonal term) and reuse d/drk
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) { 	  
	  if (l==0) 
	    W2(i*size+k,j)-=R2inv(j)*W1p(i*size+k,j);
	  if (l==1)
	    W2(i*size+j,k)-=R2inv(k)*W1p(i*size+j,k);
	  if (l==2)
	    W2(j*size+k,i)-=R2inv(i)*W1p(j*size+k,i);
	}
      }
    }
   
    // operate with d/dr_j (matrix term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) { 	  
	    if (l==0) // j sum
	      Vmat(i*size+k,j)=v((i*size+j)*size+k);
	    if (l==1) // k sum
	      Vmat(i*size+j,k)=v((i*size+j)*size+k);
	    if (l==2) // i sum
	      Vmat(j*size+k,i)=v((i*size+j)*size+k);	    
	}
      }
    }
    // blas operation 
    W1=Vmat*transpose(ddr); // i,k,j order for l=0

    // operate with 1/2Rk (diagonal term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	 for (k=0;k<size;k++) { 	   
	   if (l==0) 
	     W2(i*size+k,j)-=R2inv(k)*W1(i*size+k,j);
	   if (l==1)
	     W2(i*size+j,k)-=R2inv(i)*W1(i*size+j,k);
	   if (l==2)
	     W2(j*size+k,i)-=R2inv(j)*W1(j*size+k,i);
	 }
      }
    }

    // operate with 1/4RjRk (diagonal term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) { 
	  if (l==0)
	    W2(i*size+k,j)+=Rinv4ij(j,k)*v((i*size+j)*size+k);
	  if (l==1)
	    W2(i*size+j,k)+=Rinv4ij(k,i)*v((i*size+j)*size+k);
	  if (l==2)
	    W2(j*size+k,i)+=Rinv4ij(i,j)*v((i*size+j)*size+k);
	}
      }
    }

    // operate with extraV (diagonal term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {
	  if (l==0)
	    W2(i*size+k,j)*=extraV((i*size+j)*size+k);
	  if (l==1)
	    W2(i*size+j,k)*=extraV((j*size+k)*size+i);
	  if (l==2)
	    W2(j*size+k,i)*=extraV((k*size+i)*size+j);	  
	}
      }
    }
    
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {	  
	  if (l==0) // i sum
	    Vmat(j*size+k,i)=v((i*size+j)*size+k);
	  if (l==1) // j sum
	    Vmat(i*size+k,j)=v((i*size+j)*size+k);
	  if (l==2) // k sum
	    Vmat(i*size+j,k)=v((i*size+j)*size+k);	  
	}
      }
    }
    // operate with operate with d2/dr_idr_i + 1/r_i d/dr_i -  1/4R2 (matrix term)
    // blas operation
    W1=Vmat*transpose(ddr2invrddr); // j,k,i order for l=0
    
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {
	    if (l==0)
	      W2(i*size+k,j)+=W1(j*size+k,i);
	    if (l==1)
	      W2(i*size+j,k)+=W1(i*size+k,j);
	    if (l==2)
	      W2(j*size+k,i)+=W1(i*size+j,k);
	}
      }
    }
    
    // storage into u
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {
	  if (l==0)
	    u((i*size+j)*size+k)+=W2(i*size+k,j);
	  if (l==1)
	    u((i*size+j)*size+k)+=W2(i*size+j,k);
	  if (l==2)
	    u((i*size+j)*size+k)+=W2(j*size+k,i);
	}
      }
    }  
  }
  u=(-1./mass)*u+V*v;
  return u;
}
vector HpsiHERM(matrix &ddr2,matrix &ddr,diagmat &c1,
		diagmat &c2,diagmat &V,int size,vector &v,double mass)
{
  // t term is 
  // kin2:= (c1(ri,rj,rk)+c2(ri,rj,rk)*drj)*dri
  //      + (c1(rj,rk,ri)+c2(rj,rk,ri)*drk)*drj
  //      + (c1(rk,ri,rj)+c2(ri,rk,rj)*dri)*drk;
  //   where
  //   c1 :=(ri,rj,rk) -> (1/ri-1/4*(rk^2+ri^2-rj^2)/(rk^2*ri)-1/4*(ri^2+rj^2-rk^2)/(ri*rj^2));
  //   and 
  //   c2:=(ri,rj,rk) ->1/2*(ri^2+rj^2-rk^2)/(ri*rj);
  
  int i,j,k,l;
  vector u(size*size*size);

  matrix Umat(size*size,size);
  matrix Vmat(size*size,size);
  matrix W1(size*size,size);
  matrix W1p(size*size,size);
  matrix W2(size*size,size);

  // operate with ti term
  
  for (l=0;l<3;l++) {
    // operate with d/dr_i (matrix term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {	  
	  // prepare for blas for matrix vector product
	  // fill the V matrix
	  if (l==0) // k sum
	    Vmat(j*size+k,i)=v((i*size+j)*size+k);
	  if (l==1) // i sum
	    Vmat(i*size+k,j)=v((i*size+j)*size+k);
	  if (l==2) // j sum
	    Vmat(i*size+j,k)=v((i*size+j)*size+k);
	}
      }
    }
    // blas for matrix vector product
    W1=Vmat*transpose(ddr); // j,k,i order for l=0
    //W1=Vmat*(ddr);

    // operate with d/dr_j  (matrix term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) { 	  
	  // do some reordering
	  if (l==0) // j sum
	    W1p(i*size+k,j)=W1(j*size+k,i);
	  if (l==1) // k sum
	    W1p(i*size+j,k)=W1(i*size+k,j);
	  if (l==2) // i sum
	    W1p(j*size+k,i)=W1(i*size+j,k);
	}
      }
    }
    // blas operation
    W1=W1p*transpose(ddr); // i,k,j order for l=0
    //W1=W1p*(ddr);// i,k,j order for l=0
    
    // kin2:= (c1(ri,rj,rk)+c2(ri,rj,rk)*drj)*dri
    //      + (c1(rj,rk,ri)+c2(rj,rk,ri)*drk)*drj
    //      + (c1(rk,ri,rj)+c2(ri,rk,rj)*dri)*drk;
    // complete the operation

    // operate with 1/2Rj (diagonal term) and reuse d/drk
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) { 	  
	  if (l==0) 
	    W2(i*size+k,j)=(c1((i*size+j)*size+k)*W1p(i*size+k,j)
			    +c2((i*size+j)*size+k)*W1(i*size+k,j));
	  if (l==1)
	    W2(i*size+j,k)=(c1((j*size+k)*size+i)*W1p(i*size+j,k)
			    +c2((j*size+k)*size+i)*W1(i*size+j,k));
	  if (l==2)
	    W2(j*size+k,i)=(c1((k*size+i)*size+j)*W1p(j*size+k,i)
			    +c2((k*size+i)*size+j)*W1(j*size+k,i));
	}
      }
    }  
    // operate with d^2/dr_i^2 (matrix term)
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) { 	  
	    if (l==0) // j sum
	      Vmat(j*size+k,i)=v((i*size+j)*size+k);
	    if (l==1) // k sum
	      Vmat(i*size+k,j)=v((i*size+j)*size+k);
	    if (l==2) // i sum
	      Vmat(i*size+j,k)=v((i*size+j)*size+k);	    
	}
      }
    }
    // blas operation 
    W1=Vmat*transpose(ddr2); // j,k,i order for l=0    
    // storage into u
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {
	  if (l==0)
	    u((i*size+j)*size+k)+=(W2(i*size+k,j)+W1(j*size+k,i));
	  if (l==1)
	    u((i*size+j)*size+k)+=(W2(i*size+j,k)+W1(i*size+k,j));
	  if (l==2)
	    u((i*size+j)*size+k)+=(W2(j*size+k,i)+W1(i*size+j,k));
	}
      }
    }  
  }
  u=(-1./mass)*u+V*v;
  return u;
}
vector HpsiJac(matrix &ddr,diagmat &G11,diagmat &G22,diagmat &G33,diagmat &G12,
	       diagmat &G13,diagmat &G23,diagmat &VaplusV,int size,vector &v)
{
  // H= 1/2 d/drj^t Gjk d/drk +Va +V
  int i,j,k,l,lp;
  vector u(size*size*size);
  matrix Umat(size*size,size);
  matrix Vmat(size*size,size);
  matrix W1(size*size,size);
  matrix W1p(size*size,size);
  matrix W2(size*size,size);

  // operate with ti term  
  for (l=0;l<3;l++) {
    for (lp=0;lp<3;lp++) {
      
      for (i=0;i<size;i++) {
	for (j=0;j<size;j++) {
	  for (k=0;k<size;k++) {	  
	    // prepare for blas for matrix vector product
	    // fill the V matrix
	    if (lp==0) // k sum
	      Vmat(j*size+k,i)=v((i*size+j)*size+k);
	    if (lp==1) // i sum
	      Vmat(i*size+k,j)=v((i*size+j)*size+k);
	    if (lp==2) // j sum
	      Vmat(i*size+j,k)=v((i*size+j)*size+k);
	  }
	}
      }
      // blas for matrix vector product
      W1=Vmat*transpose(ddr); // j,k,i order for lp=0
      for (i=0;i<size;i++) {
	for (j=0;j<size;j++) {
	  for (k=0;k<size;k++) {	
	    if (lp==0) {
	      if (l==0)
		W1(j*size+k,i)*=G11((i*size+j)*size+k);	      
	      if (l==1)
		W1(j*size+k,i)*=G12((i*size+j)*size+k);
	      if (l==2)
		W1(j*size+k,i)*=G13((i*size+j)*size+k);
	    }
	    if (lp==1) {
	      if (l==0)
		W1(i*size+k,j)*=G12((i*size+j)*size+k);	      
	      if (l==1)
		W1(i*size+k,j)*=G22((i*size+j)*size+k);
	      if (l==2)
		W1(i*size+k,j)*=G23((i*size+j)*size+k);
	    }
	    if (lp==2) {
	      if (l==0)
		W1(i*size+j,k)*=G13((i*size+j)*size+k);	      
	      if (l==1)
		W1(i*size+j,k)*=G23((i*size+j)*size+k);
	      if (l==2)
		W1(i*size+j,k)*=G33((i*size+j)*size+k);
	    }
	  }
	}
      }
      for (i=0;i<size;i++) {
	for (j=0;j<size;j++) {
	  for (k=0;k<size;k++) {	
	    if (lp==0) {
	      if (l==0)
		Vmat(j*size+k,i)=W1(j*size+k,i);	      
	      if (l==1)
		Vmat(i*size+k,j)=W1(j*size+k,i);
	      if (l==2)
		Vmat(i*size+j,k)=W1(j*size+k,i);
	    }
	    if (lp==1) {
	      if (l==0)
		Vmat(j*size+k,i)=W1(i*size+k,j);	      
	      if (l==1)
		Vmat(i*size+k,j)=W1(i*size+k,j);
	      if (l==2)
		Vmat(i*size+j,k)=W1(i*size+k,j);
	    }
	    if (lp==2) {
	      if (l==0)
		Vmat(j*size+k,i)=W1(i*size+j,k);	      
	      if (l==1)
		Vmat(i*size+k,j)=W1(i*size+j,k);
	      if (l==2)
		Vmat(i*size+j,k)=W1(i*size+j,k);
	    }
	  }
	}
      }
      W1=Vmat*(ddr);
      for (i=0;i<size;i++) {
	for (j=0;j<size;j++) {
	  for (k=0;k<size;k++) {	
	      if (l==0)
		u((i*size+j)*size+k)+=W1(j*size+k,i);	      
	      if (l==1)
		u((i*size+j)*size+k)+=W1(i*size+k,j);
	      if (l==2)
		u((i*size+j)*size+k)+=W1(i*size+j,k);	   
	  }
	}
      }
    }
  }
  u=(.5)*u+VaplusV*v;
  return u;
}
void HpsiJacII(matrix &ddr,diagmat &G11,diagmat &G22,diagmat &G33,diagmat &G12,
		 diagmat &G13,diagmat &G23,diagmat &VaplusV,int size,vector &v,
		 vector &u,matrix &Vmat1,matrix &Vmat2, matrix &Vmat3,matrix &W1,
		 matrix &W2,matrix &W3)
{
  // H= 1/2 d/drj^t Gjk d/drk +Va +V
  int i,j,k,l,lp;
  //  return;

  // operate with ti term  
  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {	  
	// prepare for blas for matrix vector product
	// fill the V matrix
	u((i*size+j)*size+k)=0.;
	Vmat1(j*size+k,i)=v((i*size+j)*size+k);
	Vmat2(i*size+k,j)=v((i*size+j)*size+k);
	Vmat3(i*size+j,k)=v((i*size+j)*size+k);
      }
    }
  }
  W1=Vmat1*transpose(ddr);
  W2=Vmat2*transpose(ddr); 
  W3=Vmat3*transpose(ddr); 
  
  for (l=0;l<3;l++){  
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {	
	  if (l==0) {
	    Vmat1(j*size+k,i)=W1(j*size+k,i)*G11((i*size+j)*size+k);
	    Vmat2(j*size+k,i)=W2(i*size+k,j)*G12((j*size+i)*size+k);
	    Vmat3(j*size+k,i)=W3(i*size+j,k)*G13((i*size+j)*size+k);
	  }
	  if (l==1) {
	    Vmat1(i*size+k,j)=W1(j*size+k,i)*G12((i*size+j)*size+k);
	    Vmat2(i*size+k,j)=W2(i*size+k,j)*G22((i*size+j)*size+k);
	    Vmat3(i*size+k,j)=W3(i*size+j,k)*G23((i*size+j)*size+k);
	  }
	  if (l==2) {
	    Vmat1(i*size+j,k)=W1(j*size+k,i)*G13((k*size+j)*size+i);
	    Vmat2(i*size+j,k)=W2(i*size+k,j)*G23((i*size+k)*size+j);
	    Vmat3(i*size+j,k)=W3(i*size+j,k)*G33((i*size+j)*size+k);
	  }
	}
      }
    }
    Vmat1=Vmat1+Vmat2+Vmat3;
    Vmat1=Vmat1*(ddr);

    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {	
	  if (l==0)
	    u((i*size+j)*size+k)+=Vmat1(j*size+k,i);	      
	  if (l==1)
	    u((i*size+j)*size+k)+=Vmat1(i*size+k,j);
	  if (l==2)
	    u((i*size+j)*size+k)+=Vmat1(i*size+j,k);	   
	}
      }
    }
  }
  u=(.5)*u+VaplusV*v;
  return;
}
void Hpsipinned(matrix &ddr,diagmat &G11,diagmat &G22,diagmat &G33,diagmat &G12,
		 diagmat &G13,diagmat &G23,diagmat &VaplusV,int size,vector &v,
		 vector &u,matrix &Vmat1,matrix &Vmat2, matrix &Vmat3,matrix &W1,
		 matrix &W2,matrix &W3)
{
  // H= 1/2 d/drj^t Gjk d/drk +Va +V
  int i,j,k,l,lp;
  // return;
  // operate with ti term  
  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {	  
	// prepare for blas for matrix vector product
	// fill the V matrix
	u((i*size+j)*size+k)=0.;
	Vmat1(j*size+k,i)=v((i*size+j)*size+k);
	Vmat2(i*size+k,j)=v((i*size+j)*size+k);
	Vmat3(i*size+j,k)=v((i*size+j)*size+k);
      }
    }
  }
  W1=Vmat1*transpose(ddr);
  W2=Vmat2*transpose(ddr); 
  W3=Vmat3*transpose(ddr); 
  
  for (l=0;l<3;l++){  
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {	
	  if (l==0) {
	    Vmat1(j*size+k,i)=W1(j*size+k,i)*G11((i*size+j)*size+k);
	    Vmat2(j*size+k,i)=W2(i*size+k,j)*G12((j*size+i)*size+k);
	    Vmat3(j*size+k,i)=W3(i*size+j,k)*G13((i*size+j)*size+k);
	  }
	  if (l==1) {
	    Vmat1(i*size+k,j)=W1(j*size+k,i)*G12((i*size+j)*size+k);
	    Vmat2(i*size+k,j)=W2(i*size+k,j)*G22((i*size+j)*size+k);
	    Vmat3(i*size+k,j)=0.;
	  }
	  if (l==2) {
	    Vmat1(i*size+j,k)=W1(j*size+k,i)*G13((k*size+j)*size+i);
	    Vmat2(i*size+j,k)=0.;
	    Vmat3(i*size+j,k)=W3(i*size+j,k)*G33((i*size+j)*size+k);
	  }
	}
      }
    }
    Vmat1=Vmat1+Vmat2+Vmat3;
    Vmat1=Vmat1*(ddr);

    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {	
	  if (l==0)
	    u((i*size+j)*size+k)+=Vmat1(j*size+k,i);	      
	  if (l==1)
	    u((i*size+j)*size+k)+=Vmat1(i*size+k,j);
	  if (l==2)
	    u((i*size+j)*size+k)+=Vmat1(i*size+j,k);	   
	}
      }
    }
  }
  u=(.5)*u+VaplusV*v;
  return;
}
void Hpsi2d(matrix &ddr,diagmat &G22,diagmat &G33
	    ,diagmat &G23,diagmat &VaplusV,int size,vector &v,
	    vector &u,matrix &Vmat1,matrix &Vmat2, matrix &Vmat3,matrix &W1,
	    matrix &W2,matrix &W3)
{
  // H= 1/2 d/drj^t Gjk d/drk +Va +V
  int i,j,k,l,lp;
  //  return;
  
  // operate with ti term  
  for (j=0;j<size;j++) {
    for (k=0;k<size;k++) {	  
      // prepare for blas for matrix vector product
      // fill the V matrix
      u((j)*size+k)=0.;
      Vmat2(k,j)=v((j)*size+k);
      Vmat3(j,k)=v((j)*size+k);
    }
  }
  W2=Vmat2*transpose(ddr); 
  W3=Vmat3*transpose(ddr); 
  
  for (l=1;l<3;l++){  
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {	
	if (l==1) {
	  Vmat2(k,j)=W2(k,j)*G22((j)*size+k);
	  Vmat3(j,k)=0.;
	}
	if (l==2) {
	  Vmat3(j,k)=W3(j,k)*G33((j)*size+k);
	  Vmat2(j,k)=0.;
	}
      }
    }
    Vmat1=Vmat2+Vmat3;
    Vmat1=Vmat1*(ddr);

    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {		      
	if (l==1)
	  u((j)*size+k)+=Vmat1(k,j);
	if (l==2)
	  u((j)*size+k)+=Vmat1(j,k);	   
      }
    }
  }
  u=(.5)*u+VaplusV*v;
  return;
}
void HpsiJacCont(matrix &ddr,diagmat &G11,diagmat &G22,diagmat &G33,diagmat &G12,
		 diagmat &G13,diagmat &G23,diagmat &VaplusV,int size,vector &v,
		 vector &u,matrix &Vmat1,matrix &Vmat2, matrix &Vmat3,matrix &W1,
		 matrix &W2,matrix &W3,int nzeroes)
{
  // H= 1/2 d/drj^t Gjk d/drk +Va +V
  int i,j,k,l,lp;
  int sizeC=size-nzeroes;
  //  return;
  vector vtemp(size*size*size);
  matrix EV(size,size);
  for (i=0;i<size;i++) 
    for (j=0;j<size;j++) EV(i,j)=Vmat1(i,j);

  // transform v to the DVR
  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {
	Vmat1(i*size+j,k)=0.;
	for (l=0;l<size;l++) {	  
	  Vmat1(i*size+j,k)+=EV(i,l)*v((l*size+j)*size+k);
	}
      }
    }
  }
  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {
	W1((i*size+k),j)=0.;
	for (l=0;l<size;l++) {	  
	  W1((i*size+k),j)+=EV(j,l)*Vmat1(i*size+l,k);
	}
      }
    }  
  }
  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {
	vtemp((i*size+j)*size+k)=0.;
	for (l=0;l<size;l++) {	  
	  vtemp((i*size+j)*size+k)+=EV(k,l)*W1(i*size+l,j);
	}
      }
    }  
  }


  // operate with ti term  
  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {	  
	// prepare for blas for matrix vector product
	// fill the V matrix
	u((i*size+j)*size+k)=0.;
	Vmat1(j*size+k,i)=vtemp((i*size+j)*size+k);
	Vmat2(i*size+k,j)=vtemp((i*size+j)*size+k);
	Vmat3(i*size+j,k)=vtemp((i*size+j)*size+k);
      }
    }
  }
  W1=Vmat1*transpose(ddr);
  W2=Vmat2*transpose(ddr); 
  W3=Vmat3*transpose(ddr); 
  
  for (l=0;l<3;l++){  
    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {	
	  if (l==0) {
	    Vmat1(j*size+k,i)=W1(j*size+k,i)*G11((i*size+j)*size+k);
	    Vmat2(j*size+k,i)=W2(i*size+k,j)*G12((j*size+i)*size+k);
	    Vmat3(j*size+k,i)=W3(i*size+j,k)*G13((i*size+j)*size+k);
	  }
	  if (l==1) {
	    Vmat1(i*size+k,j)=W1(j*size+k,i)*G12((i*size+j)*size+k);
	    Vmat2(i*size+k,j)=W2(i*size+k,j)*G22((i*size+j)*size+k);
	    Vmat3(i*size+k,j)=W3(i*size+j,k)*G23((i*size+j)*size+k);
	  }
	  if (l==2) {
	    Vmat1(i*size+j,k)=W1(j*size+k,i)*G13((k*size+j)*size+i);
	    Vmat2(i*size+j,k)=W2(i*size+k,j)*G23((i*size+k)*size+j);
	    Vmat3(i*size+j,k)=W3(i*size+j,k)*G33((i*size+j)*size+k);
	  }
	}
      }
    }
    Vmat1=Vmat1+Vmat2+Vmat3;
    Vmat1=Vmat1*(ddr);

    for (i=0;i<size;i++) {
      for (j=0;j<size;j++) {
	for (k=0;k<size;k++) {	
	  if (l==0)
	    u((i*size+j)*size+k)+=Vmat1(j*size+k,i);	      
	  if (l==1)
	    u((i*size+j)*size+k)+=Vmat1(i*size+k,j);
	  if (l==2)
	    u((i*size+j)*size+k)+=Vmat1(i*size+j,k);	   
	}
      }
    }
  }
  // transform u from the DVR
  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {
	Vmat1(i*size+j,k)=0.;
	for (l=0;l<size;l++) {	  
	  Vmat1(i*size+j,k)+=EV(l,i)*u((l*size+j)*size+k);
	}
      }
    }
  }
  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {
	W1((i*size+k),j)=0.;
	for (l=0;l<size;l++) {	  
	  W1((i*size+k),j)+=EV(l,j)*Vmat1(i*size+l,k);
	}
      }
    }  
  }
  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {
	u((i*size+j)*size+k)=0.;
	for (l=0;l<size;l++) {	  
	  u((i*size+j)*size+k)+=EV(l,k)*W1(i*size+l,j);
	}
      }
    }  
  }

  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {
	Vmat1(i*size+j,k)=0.;
	for (l=0;l<size;l++) {	  
	  Vmat1(i*size+j,k)+=EV(i,l)*v((l*size+j)*size+k);
	}
      }
    }
  }
  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {
	W1((i*size+k),j)=0.;
	for (l=0;l<size;l++) {	  
	  W1((i*size+k),j)+=EV(j,l)*Vmat1(i*size+l,k);
	}
      }
    }  
  }
  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {
	vtemp((i*size+j)*size+k)=0.;
	for (l=0;l<size;l++) {	  
	  vtemp((i*size+j)*size+k)+=EV(k,l)*W1(i*size+l,j);
	}
      }
    }  
  }
  vtemp=VaplusV*vtemp;

  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {
	Vmat1(i*size+j,k)=0.;
	for (l=0;l<size;l++) {	  
	  Vmat1(i*size+j,k)+=EV(l,i)*vtemp((l*size+j)*size+k);
	}
      }
    }
  }
  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {
	W1((i*size+k),j)=0.;
	for (l=0;l<size;l++) {	  
	  W1((i*size+k),j)+=EV(l,j)*Vmat1(i*size+l,k);
	}
      }
    }  
  }
  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      for (k=0;k<size;k++) {
	vtemp((i*size+j)*size+k)=0.;
	for (l=0;l<size;l++) {	  
	  vtemp((i*size+j)*size+k)+=EV(l,k)*W1(i*size+l,j);
	}
      }
    }  
  }

  u=(.5)*u+vtemp;
  return;
}
double gaussnorm(double r,double alpha,double r0)
{
  return sqrt(1./(2.*M_PI*alpha))*exp(-.5*(r-r0)*(r-r0)/alpha);
}
double gauss(double r,double mu,double w,double r0)
{
  return exp(-.5*mu*w*w*(r-r0)*(r-r0));
}
matrix DJacobi(double a,double b,int size)
{
  int m,n;
  matrix D(size,size);
  for (m=0;m<size;m++) {
    double dm=(double)m;
    for (n=0;n<size;n++) {
      double dn=(double)n;
      if (m >= n)
	D(m,n)=dmn(dm,dn,a,b);
      if (m <= n)
	D(m,n)=.5*cdown(dm,dn,a,b)*(del(a)-pow(-1.,dm-dn)*del(b))-dmn(dn,dm,a,b);
    }
  }
  return D;
}
matrix BJacobi(double a,double b,int size)
{
  int m,n;
  matrix D(size,size);
  for (m=0;m<size;m++) {
    double dm=(double)m;
    for (n=0;n<size;n++) {
      double dn=(double)n;
	D(m,n)=.5*(a+b+2.*dn+2.)*q(dn,a,b)*delta(m+1,n)
	  +d(dn,a,b)*delta(m,n)
	  -.5*(a+b+2.*dn)*q(dn+1.,a,b)*delta(m,n+1);
    }
  }
  return D;
}
double dmn(double dm,double dn,double a,double b)
{
  double d;
  d=.25*(delbar(b)*pow(-1.,dm-dn)*cboth(dm,a,dn,b)-
	      delbar(a)*cboth(dm,b,dn,a));
  return d;
}
double del(double a)
{
  if (a == 0.) 
    return 1.;
  else
    return 0.;
}
double delbar(double a)
{
  return (1.-del(a));
}
double cdown(double dm,double dn,double a,double b)
{
  return sqrt((a+b+2.*dm+1.)*(a+b+2.*dn+1.));
}
double cboth(double dm,double b,double dn,double a)
{
  double c;
  
  c=cdown(dm,dn,a,b)*sqrt( zn(dn+1.,dm-dn)/zn(a+b+dn+1.,dm-dn)*
			   zn(b+dn+1.,dm-dn)/zn(a+dn+1.,dm-dn));
  return c;
}
double zn(double z,double dn)
{
  double zout,lz;
  if (dn > 0.) {
    lz=lgamma(z+dn)-lgamma(z);
    zout=exp(lz);
    int n=(int)dn;
//     zout=z;
//     for (int i=1;i<n;i++)
//       zout*=(z+(double)i);
  }
  if (fabs(dn) <= 1.e-14)
    zout=1.;
  if (dn < 0.) {
    lz=lgamma(z-dn)-lgamma(z);
    zout=exp(lz);
//     zout=1./zn(z+dn,-dn);
  }
  return zout;
}
matrix xJacobi(double a,double b,int size)
{
  int m,n;
  double xel;
  matrix x(size,size);
  for (m=0;m<size;m++) {
    double dm=(double)m;
    for (n=0;n<size;n++) {
      double dn=(double)n;
      if (m >= n)
	xel=d(dn,a,b)*delta(m,n)+q(dn,a,b)*delta(m+1,n)+q(dm,a,b)*delta(m,n+1);
      else
	xel=d(dm,a,b)*delta(m,n)+q(dm,a,b)*delta(m,n+1)+q(dn,a,b)*delta(m+1,n);
      x(m,n)=xel;
    }
  }
  return x;
}
double d(double dn,double a,double b)
{
  if (a == b) return 0.;
  else
    return (b*b-a*a)/(pow(a+b+2.*dn+1.,2.)-1.);
}
double q(double dn,double a,double b)
{
  if (a == b && a==0.)
    return dn/sqrt(4.*dn*dn-1.);
  else
    return (2./(a+b+2.*dn))*
      sqrt(dn*(a+dn)*(b+dn)*(a+b+dn)/(pow(a+b+2.*dn,2.)-1.));
}
void PIA1(vector &v,int size) {
  int i,j,k;
  vector u(size*size*size);    
  for (i=0;i<size;i++) 
    for (j=0;j<size;j++) 
      for (k=0;k<size;k++) {   
	u((i*size+j)*size+k)=(.16666666666666666666)*(v((i*size+j)*size+k)+v((j*size+i)*size+k)
	  +v((k*size+j)*size+i)+v((k*size+i)*size+j)
	  +v((i*size+k)*size+j)+v((j*size+k)*size+i));
      }
  v=u;
  return;
}
void PIA2(vector &v,int size) {
  int i,j,k;
  vector u(size*size*size);    
  for (i=0;i<size;i++) 
    for (j=0;j<size;j++) 
      for (k=0;k<size;k++) {   
// 	//do S3 instead
	u((i*size+j)*size+k)=(.16666666666666666666)*(v((i*size+j)*size+k)  // E
					+v((j*size+k)*size+i) // (123)
					+v((k*size+i)*size+j) // (132)
					-v((j*size+i)*size+k) // (12)
					-v((k*size+j)*size+i) // (13)
					-v((i*size+k)*size+j)); // (23)

      }
  v=u;
  return;
}
void PIE(vector &v,int size) {
  int i,j,k;
  vector u(size*size*size);    
  for (i=0;i<size;i++) 
    for (j=0;j<size;j++) 
      for (k=0;k<size;k++) {   
// 	u((i*size+j)*size+k)=(1./(6.))*(2.*v((i*size+j)*size+k)  // E
// 				      -v((j*size+k)*size+i) // (123)
// 				      +0.*v((i*size+k)*size+j) // (23)
// 				      +2.*v((i*size+j)*size+k) // (E*)
// 				      -v((j*size+k)*size+i) // (123)*
// 				      +0.*v((i*size+k)*size+j)); // (23)*
	// using S3 below 
	u((i*size+j)*size+k)=(.33333333333333333333)*(2.*v((i*size+j)*size+k)  // E
					-v((j*size+k)*size+i) // (123)
					-v((k*size+i)*size+j)); // (132)
      }
  v=u;
  return;
}
double Aziz(double r)
{
  // parameters for the LM2M2 potential of Aziz JCP, 94, 8047 (1991)
  double AStar=1.89635353e5;
  double alphaStar=10.70203539;
  double c6=1.34687065;
  double c8=0.41308398;
  double c10=0.17060159;
  double C6=1.461; // au
  double C8=14.11; // au
  double C10=183.5; // au
  double betaStar=-1.90740649;
  double beta=-0.21631; // A^-2
  double D=1.4088;
  double epsilonoverk=10.97; // K
  double epsilon=epsilonoverk*kB;
  double rm=2.9695; // A
  rm*=atob;
  double sigma= 2.6417278; // A
  // add-on portions
  double Aa=0.0026;
  double x1=1.003535949;
  double x2=1.454790369;

  // HFD potential
  // aziz and chen, JCP, 67, 5719 (1977)
  // V(x) = epsilon V*(x)
  // V*(x) = A exp(-alpha x) - (C6/x^6+C8/x^8+C10/x^10) F(x)
  // F(x) = exp(-(D/x-1)^2) for x<D and F(x)=1 for >=D


  double x=r/rm;
  double B=2.*M_PI/(x2-x1);
  double VaStar=Aa*(sin(B*(x-x1)-M_PI/2.)+1.);
  if (x <x1 || x > x2) VaStar=0.;

  double F= exp(-pow((D/x-1.),2.));
  if (x >= D) F=1.;

  
  double VbStar=AStar*exp(-alphaStar*x+betaStar*x*x)-(c6/pow(x,6.)+c8/pow(x,8.)+c10/pow(x,10.))*F;
  
  double V=epsilon*(VaStar+VbStar);
  return V;
}
double HeHminus(double r)
{
	double a1=1.303648;
	double a2=-1.297418;
	double a3=0.7503146;
	double a4=-0.00001989;
	double a5=0.000012769;
	double a6=0.313155;
	double a7=13.;
	double a8=0.00001467;
	double a9=0.736841;
	double V;
	if ( r< 10.)
		V=a1*pow(r,a2)*exp(-a3*r)+a4;
	if (r >=10. && r<= 20.)
		V=a5*pow(1.-exp(-a6*(r-a7)),2.)-a8;
	if (r >20.)
		V=-a9*pow(1./r,4.);
	return V;
}
double Hebuckypot(double r,int Ncarbons)
{
  double epsilonoverk=16.92;//kelvin
  double epsilon=epsilonoverk*kB;
  double sigma=2.98; sigma*=atob;
  double R=0.4583*sqrt((double)Ncarbons); R*=atob;
  double V=((double)Ncarbons)*epsilon*(1./(r*R))*
    ((pow(sigma,12.)/5.)*(pow(r-R,-10.)-pow(r+R,-10.))
     -(pow(sigma,6.)/2.)*(pow(r-R,-4.)-pow(r+R,-4.)));
  if (V >= 1.E+13) V=1.E+13;
  if (r<=R) V=1.E+13;
  return V;
}

void addtobasis(vector &v,matrix &B,int n,int col)
{
  for (int i=0;i<n;i++) B(i,col)=v(i);
  return;
}
double lj(double r, double sigma, double epsilon)
{
  double v=4.*epsilon*(pow(sigma/r,12.)-pow(sigma/r,6.));
  return v;
}
matrix CMKE(int nsize,double dx, double mass)
{
  matrix tmat(nsize,nsize);
  int ii,iip;
  for (int i=0;i<nsize;i++) {
    for (int j=0;j<nsize;j++) {
      ii=i-(nsize-1)/2;
      iip=j-(nsize-1)/2;
      if (i == j)
        tmat(i,j)=M_PI*M_PI/(6.*mass*dx*dx);
      else if((ii-iip)%2 == 0)
        tmat(i,j)=1./(((double) (ii-iip)*(ii-iip))*mass*
                      dx*dx);
      else
        tmat(i,j)=(-1./(((double) (ii-iip)*(ii-iip))*mass
                        *dx*dx));
    }
  }
  return tmat;
}
vector CMgrid(int size,double length,double Rmin)
{
  vector grid(size);
  double deltaR=(length)/((double)size-1.);
  vector Rgrid(size);
  for (int i=0;i<size;i++) 
    grid(i)=Rmin+(double)i*deltaR;
  return grid;
}
// int arpackinterface(matrix &ddr,diagmat &G11,diagmat &G22,diagmat &G33,diagmat &G12,
// 		     diagmat &G13,diagmat &G23,diagmat &VaplusV,int size,
// 		     int niter,double emin,
// 		     double emax,int maxnev,double tolerance,int labelsym,
// 		     vector &v,vector &d)
// {
//     int ngood;
//     int maxn=size*size*size;
//     int maxncv=2*maxnev;
//     vector workl(maxncv*(maxncv+8));
//     vector workd(3*maxn); 
//     vector resid(maxn);
//     vector ax(maxn);
//     int *select=new int[maxncv];
//     vector workperm(maxn);
    
//     FORTRAN(arpack)(ddr.TheMatrix,G11.TheMatrix,G22.TheMatrix,G33.TheMatrix,G12.TheMatrix,
// 		    G13.TheMatrix,G23.TheMatrix,VaplusV.TheMatrix,&size,
// 		    &maxnev,&maxncv,&emin,&emax,&maxn,
// 		    v.TheVector, workl.TheVector, workd.TheVector, d.TheVector, resid.TheVector,
// 		    ax.TheVector,select,&niter,&tolerance,&labelsym,workperm.TheVector,&ngood);
    
    
//     return ngood;
// }
