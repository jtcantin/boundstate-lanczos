#include <iostream>
#include <cmath>
#include "vectClass.h" //This is to select which linear algebra class should be used.

#include "peckeris.h"

static double betaGlobal;
void lanczosvectors(VECT &alpha,VECT &beta,VECT &beta2,int niter,
					VECT &eval,int ngood,MAT &evtr);

static const double amutoau=1./5.4857989586762187e-4; // amu to au
static const double au_to_amu = 5.4857989586762187e-4; // au to amu
static const double BohrToA=0.529177249; // Bohr to angstrom
static const double HtoW=2.19474631371017e5;//hatree to reciprical centimeter
static const double WavToMHz=29979.2458;
static const double hatokJmol=2625.;
static const double MofH2 = 2.015650642; // nist value
static const double MofD2 = 2.0141017780*2.;// nist
static const double MofHe4 = 4.0026032497; // nist mass of He in amu


int main(int argc,char **argv) {
	int i,j,k,row,n;
	
	int niter=atoi(argv[2]); //Number of iterations
	
	
	
	
	//Set up what is needed for any Hv calculations, such as pre-calculating the potential and any other matricies
	FORTRAN(prepare)(&numbas,icode,&maxj,&jmax,&Jbig,&ik,&ip,Kbgind,
					 coriol,rotor,&Ah2o,&Bh2o,&Ch2o,&h2_2mu,ijcori,
					 &ncorio,ijrot,&nrotor,&maxfac,fact,Tcos,Tsin,
					 &nlgrid,&ncgrid,wgtgl,glgrid,gcgrid,&nchi,&nthe,
					 &nrad,vpes,gchi,gthe,grad,&nrpont,&rsmall,&rlarge,
					 rgrid,&nrgrid,vmat,rkin,&nKbig,&ib000,&ib100,&ib101,
					 &ib110,&ib111,&jb111,potfil);
	
	
	double* vec=new double[nrpont*numbas];
	double* uec=new double[nrpont*numbas];
	
	//test Hv
	
	double emin=-100.;
	cout<<"emin= "<<emin<<endl;
	double emax=200.;
	//double emax=-emin;
	
	int ntotbs=nrpont*numbas;
	
	int ngood;
	
	VECT evalerr(niter);
	VECT eval(niter);
	VECT alpha(niter);
	VECT beta(niter+1);
	VECT  beta2(niter+1);
	
	double* rvec=new double[nrpont*numbas];  
	for (i=0;i<(nrpont*numbas);i++) rvec[i]=0.;
	
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Lanczos Loop 1 - Get Tm
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout<<"Lanczos Loop 1 - Get Tm beginning."<<endl;
	
	
	// Initialize Vectors
	for (int ib=0;ib<ntotbs;ib++) {
		uec[ib]=0.;
		vec[ib]=1./sqrt((double)ntotbs); //v_1 (normalized equal superposition)
		rvec[ib]=0.;
	}
	cout << "Lanczos Coefficients"<<endl;
	cout<<"Iteration j"<<" "<<"alpha_j"<<" "<<"beta_j"<<endl;
	
	//Lanczos Loop, run over given number of iterations
	for (j=1;j<=niter;j++) {    
		
		//u_j = 0. -> NOTE: Eliminate this by changing Hv output to be double*; the current method is a waste of computational time
		for (int ib=0;ib<ntotbs;ib++) 
			uec[ib]=0.;
		
		//u_j = u_j + H*v_j
		FORTRAN(hv)(vec,uec,vmat,&ncgrid,&nlgrid,rgrid,
					&nrpont,&numbas,&Jbig,&ip,Kbgind,icode,Tcos,
					Tsin,&nKbig,coriol,ijcori,rotor,ijrot,
					&ncorio,&nrotor,rkin); 
		
		//r_j = r_j + u_j; ; results in r_j = (-beta_(j-1) * v_(j-1)) + (H*v_j) = H*v_j - beta_(j-1) * v_(j-1)
		for (int ib=0;ib<ntotbs;ib++) 
			rvec[ib]+=uec[ib];
		
		//alpha_j = v_j * r_j
		alpha(j-1)=0.;
		for (int ib=0;ib<ntotbs;ib++) 
			alpha(j-1)+=vec[ib]*rvec[ib];
		
		//r_j = r_j - alpha_j * v_j
		for (int ib=0;ib<ntotbs;ib++) 
			rvec[ib]-=(alpha(j-1)*vec[ib]);
		
		//beta_j = ||r_j||
		beta2(j)=0.;
		for (int ib=0;ib<ntotbs;ib++) 
			beta2(j)+=rvec[ib]*rvec[ib];
		
		beta(j)=sqrt(beta2(j));
		
		
		//r_j = r_j/beta_j (normalize r_j) = v_(j+1) (after swap)
		for (int ib=0;ib<ntotbs;ib++) 
			rvec[ib]=(1./beta(j))*rvec[ib];
		
		//v_j = -beta_j * v_j = r_(j+1) (after swap)
		for (int ib=0;ib<ntotbs;ib++) 
			vec[ib]=(-beta(j))*vec[ib]; // prepare r check minus sign!!!
		
		//Swap vectors
		for (int ib=0;ib<ntotbs;ib++) {
			uec[ib]=vec[ib];	//u_(j+1) = v_j -> NOTE: u_(j+1) is set to zero every loop; this is a waste of time!
			vec[ib]=rvec[ib];	//v_(j+1) = r_j; results in v_(j+1) = r_j/beta_j
			rvec[ib]=uec[ib];   //r_(j+1) = u_j = v_j; results in r_(j+1) = -beta_j * v_j -> NOTE: this indirect notation is a waste of computational time!
		}
		cout<<j<<" "<<alpha(j-1)<<" "<<beta(j)<<endl;
		/*if (j%100 == 0)
			cout<<"iteration "<<j<<endl;*/
	}                  
	
	cout << "Lanczos Loop 1 - Get Tm FINISHED." << endl;
	cout << "--------------------------------------------------------------------------------" << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Determine Good Eigenvalues
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	cout << endl;
	cout << "Determination of Good Eigenvalues beginning." << endl;
	
	//Determine Good Eigenvalues
	lancbis(niter,eval,evalerr,emin,emax,ngood,alpha,beta,beta2);
	
	cout << endl;
	cout << "Number of Good Eigenvalues is: " << ngood << endl;
	cout << "Groundstate Eigenvalue is: " << eval(0) << " kJ/mol" << endl;
	
	// Output Eigenvalues and Errors
	ofstream lancout("boundstates.out");
	ofstream lanczpeout("states_zpe.out");
	
	//Write header for Eigenvalues and Errors
	lancout << "The following are the eigenvalues and errors for the simulation " << sim_descr << endl;
	lancout << "Number of Eigenvalues: " << ngood << endl;
	lancout << "State | Eigenvalue (kJ/mol) | Error (kJ/mol?)" << endl;
	
	//Write header for Relative Energies
	lanczpeout << "The following are the relative energies for the simulation " << sim_descr << endl;
	lanczpeout << "Number of Eigenvalues: " << ngood << endl;
	lanczpeout << "Zero Point Energy = " << eval(0) << " kJ/mol" << endl;
	lanczpeout << "State | Relative Energy (kJ/mol)" << endl;
	
	//Output Values
	for (i=0;i<ngood;i++) {
		lancout << i << " " << eval(i) << " " << evalerr(i) << endl;
		lanczpeout << i << " " << (eval(i)-eval(0)) << endl;
    }
	
	lancout.flush();
	lancout.close();
	
	lanczpeout.flush();
	lanczpeout.close();
	
	cout << "Eigenvalues have been calculated." << endl;
	cout << "The Eigenvalues and their errors have been stored in " << "boundstates.out" << endl;
	cout << "The Relative Energies have been stored in " << "states_zpe.out" << endl;
	cout << "--------------------------------------------------------------------------------" << endl;
	
	// If only the eigenvalues are wished to be calculated, end the program here.
	if !(calcEigVec) {
		cout << "LANCZOS ALGORITHM COMPLETED" << endl;
		
		exit(0);
	}
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Determine the Tm Eigenvectors
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//If not all of the eigenvectors are desired, generate only the desired number.
	if (numEigVecWord != "all") {
		ngood = numEigVec;
	}
    
	cout << endl;
	cout << "Eigenvector calculation beginning." << endl;
	
	//Calculate the eigenvectors of Tm
	MAT evtr(niter,ngood);
	
	cout << "Tm eigenvector calculation beginning." << endl;
	
	//Calculate Tm eigenvectors and store in evtr
	lanczosvectors(alpha,beta,beta2,niter,eval,ngood,evtr); 
	
	cout << "Tm eigenvector calculation finished." << endl;
	
	//Output Tm eigenvectors to a file
	ofstream lvecout("lv");
	VECT cumulnorm(ngood);
	double coeff2;
	
	lvecout << "The following are the Coefficients squared of the Eigenvectors of Tm for the simulation " << sim_descr << endl;
	lvecout << "Number of Eigenvectors: " << ngood << endl;
	lvecout << "Iteration | Coefficient Squared | Cumulative Norm " << endl;
	
	for (j=1; j<=niter; j++) {	
		lvecout << j << " ";
		
		for (n=0; n<ngood; n++) {
			//Calculate the square of the jth coefficient of each nth vector
			coeff2 = pow(evtr(j-1,n),2.);
			lvecout << coeff2 << " ";
			
			//Cumulative norm
			cumulnorm(n) += coeff2;
			lvecout << cumulnorm(n) << " ";
		}
		
		lvecout << endl;
	}
    cout << "The Coefficients squared of the Eigenvectors of Tm have been stored in " << "lv" << endl;
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Lanczos Loop 2 - Get Eigenvectors of H
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	cout<<"Lanczos Loop 2 - Get Eigenvectors of H beginning."<<endl;
	
	//Reinitialize u_j, v_j, and r_j
	for (int ib=0;ib<ntotbs;ib++) {
		uec[ib]=0.;
		vec[ib]=1./sqrt((double)ntotbs);
		rvec[ib]=0.;
	}
	
	//Reinitialize alpha_j, beta_j, and beta2_j
	for (j=0;j<niter;j++) {
		alpha(j)=0.;
		beta(j)=0.;
		beta2(j)=0.;
		
	}
	beta(niter)=0.; //beta_j and beta2_j run up to niter+1
	beta2(niter)=0.;
	
	VECT ARvL(ntotbs*ngood); //This stores the final Eigenvectors (Ritz vectors)
		
	//Reset the cumulative norm
	for (n=0;n<ngood;n++) cumulnorm(n)=0.;
	
	double treshold;
	double coeff;
	
	cout << "Lanczos Coefficients" << endl;
	cout<<"Iteration j"<<" "<<"alpha_j"<<" "<<"beta_j"<<endl;
	
	//Begin the 2nd Lanczos Loop
	for (j=1; j<=niter; j++) {	
		
		// Effectively calculate ARvL_rn = Ym_rn = Vm_rj * Em_jn as Ym_rn += Em_jn * v_j(r) summed over j
		// Em is the matrix holding the Tm eigenvectors
		for (n=0;n<ngood;n++) {
			
			treshold = pow(evtr(j-1,n),2.);
			cumulnorm(n) += treshold;
			
			for (row=0;row<ntotbs;row++){
				
				coeff=evtr(j-1,n)*vec[row];	  //coeff = Em_jn * v_j(r)
				
				if (cumulnorm(n) < (1. - 1.e-16)) {//If the cumulative norm is effectively one, don't add the coefficient anymore
					ARvL(row+ntotbs*n) += coeff; //NOTE: This is summed over j, not row
				}
			}
		}
		
		//u_j = 0. -> NOTE: Eliminate this by changing Hv output to be double*; the current method is a waste of computational time
		for (int ib=0;ib<ntotbs;ib++) 
			uec[ib]=0.;
		
		//u_j = u_j + H*v_j
		FORTRAN(hv)(vec,uec,vmat,&ncgrid,&nlgrid,rgrid,
					&nrpont,&numbas,&Jbig,&ip,Kbgind,icode,Tcos,
					Tsin,&nKbig,coriol,ijcori,rotor,ijrot,
					&ncorio,&nrotor,rkin);
		
		//r_j = r_j + u_j; ; results in r_j = (-beta_(j-1) * v_(j-1)) + (H*v_j) = H*v_j - beta_(j-1) * v_(j-1)
		for (int ib=0;ib<ntotbs;ib++) 
			rvec[ib]+=uec[ib];
		
		//alpha_j = v_j * r_j
		alpha(j-1)=0.;
		for (int ib=0;ib<ntotbs;ib++) 
			alpha(j-1)+=vec[ib]*rvec[ib];
		
		//r_j = r_j - alpha_j * v_j
		for (int ib=0;ib<ntotbs;ib++) 
			rvec[ib]-=(alpha(j-1)*vec[ib]);
		
		//beta_j = ||r_j||
		beta2(j)=0.;
		for (int ib=0;ib<ntotbs;ib++) 
			beta2(j)+=rvec[ib]*rvec[ib];
		
		beta(j)=sqrt(beta2(j));
		
		
		//r_j = r_j/beta_j (normalize r_j) = v_(j+1) (after swap)
		for (int ib=0;ib<ntotbs;ib++) 
			rvec[ib]=(1./beta(j))*rvec[ib];
		
		//v_j = -beta_j * v_j = r_(j+1) (after swap)
		for (int ib=0;ib<ntotbs;ib++) 
			vec[ib]=(-beta(j))*vec[ib]; // prepare r check minus sign!!!
		
		//Swap vectors
		for (int ib=0;ib<ntotbs;ib++) {
			uec[ib]=vec[ib];	//u_(j+1) = v_j -> NOTE: u_(j+1) is set to zero every loop; this is a waste of time!
			vec[ib]=rvec[ib];	//v_(j+1) = r_j; results in v_(j+1) = r_j/beta_j
			rvec[ib]=uec[ib];	//r_(j+1) = u_j = v_j; results in r_(j+1) = -beta_j * v_j -> NOTE: this indirect notation is a waste of computational time!
		}
		
		cout << j << " " << alpha(j-1) << " " << beta(j) << endl;
		/*if (j%100 == 0)
		 cout<<"iteration "<<j<<endl;*/
	}    
	
	
	// eigenvectors for good eigenvalues are in ARvL(row+ntotbs*n)
	
	// for the ground state, eg, ARvL(row+ntotbs*0) is the vector in the original basis

	//cout<<ARvL(0)<<" "<<ARvL(ntotbs*ngood-1)<<" "<<icode[0]<<" "<<icode[numbas-1]<<endl;
	
	//Output eigenvectors to file.	
	ofstream eigVecFile("eigVecFile.out")
	
	eigVecFile << "The following are the Ritz vectors of H for the simulation " << sim_descr << endl;
	eigVecFile << "Number of Ritz vectors: " << ngood << endl;
	eigVecFile << "Basis Index | Vector 1 | Vector 2 ... " << endl;
	
	for (row=0; row<ntotbs; row++) {
		eigVecFile << row << " ";
		
		for (n=0; n<ngood; n++) {
			eigVecFile << ARvL(row+ntotbs*n) << " ";
		}
		
		eigVecFile << endl;
	}
	
	cout << "The Eigenvectors of H have been stored in " << "eigVecFile.out" << endl;
	
	cout << "Lanczos Loop 2 - Get Eigenvectors of H FINISHED." << endl;
	cout << "--------------------------------------------------------------------------------" << endl;
	
	cout << "Boundstate Lanczos Calculator has finished executing." << endl;
	cout << "Have a nice day." << endl;
	
	return 0;
}


void lanczosvectors(VECT &alpha, VECT &beta, VECT &beta2, int niter, VECT &eval, int ngood, MAT &evtr) {
	int i,j,ndis;
	
    // copy stuff
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
	
	//Calculate the Tm eigenvectors
    FORTRAN(trivec)(lalpha,lbeta,lbeta2,wrk1,wrk2,&niter,leval,&ngood,levtr,mamin);
	
	//Copy the eigenvectors to the matrix evtr
    for (i=0;i<niter;i++)
		for (j=0;j<ngood;j++) 
			evtr(i,j)=levtr[i+j*niter];
    return;
}

