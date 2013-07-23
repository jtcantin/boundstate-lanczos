#include "lanczos.h"

void lanczosvectors(VECT &alpha,VECT &beta,VECT &beta2,int niter,
					VECT &eval,int ngood,MAT &evtr);
/*
static const double amutoau=1./5.4857989586762187e-4; // amu to au
static const double au_to_amu = 5.4857989586762187e-4; // au to amu
static const double BohrToA=0.529177249; // Bohr to angstrom
static const double HtoW=2.19474631371017e5;//hatree to reciprical centimeter
static const double WavToMHz=29979.2458;
static const double hatokJmol=2625.;
static const double MofH2 = 2.015650642; // nist value
static const double MofD2 = 2.0141017780*2.;// nist
static const double MofHe4 = 4.0026032497; // nist mass of He in amu
*/

int main(int argc,char **argv) {
	int i,j,row,n;
	
	int niter, numEigVec;
	double emin, emax;
	string	numEigVecWord, HvCalculatorSwitch;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Get information from Lanczos Algorithm Input File
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	string lanczosInputFilename, junk, line, outputDirPath, HvInputFilename;
	ifstream lanczosInputFile;
	
	lanczosInputFilename = argv[1];
	HvInputFilename = argv[2];
	
	lanczosInputFile.open(lanczosInputFilename.c_str(), ios::in);
	if (lanczosInputFile.is_open()) {
		
		//Get rid of comment lines;
		getline(lanczosInputFile, line);
		getline(lanczosInputFile, line);
		
		//Gather data for input; each line has the name of the value separated from the value with a space
		//Always ignore the name of the value.
		
		//Number of Iterations
		lanczosInputFile >> junk;
		lanczosInputFile >> niter;
		
		//EigenvalueOpenLowerLimit(kJ/mol)
		lanczosInputFile >> junk;
		lanczosInputFile >> emin;
		
		//EigenvalueClosedUpperLimit(kJ/mol)
		lanczosInputFile >> junk;
		lanczosInputFile >> emax;
		
		//NumberOfEigenvectors(all/partial/none)
		lanczosInputFile >> junk;
		lanczosInputFile >> numEigVecWord;
		
		//NumberOfEigenvectors
		lanczosInputFile >> junk;
		lanczosInputFile >> numEigVec;
		
		//HvCalculator
		lanczosInputFile >> junk;
		lanczosInputFile >> HvCalculatorSwitch;
		
		//Output directory path
		lanczosInputFile >> junk;
		lanczosInputFile >> outputDirPath;

	}
	else {
		cerr << "ERROR: Lanczos input file '" << lanczosInputFilename << "' could not be opened." << endl;
		exit(1);
	}
	
	if ((numEigVecWord != "all") && (numEigVecWord != "partial") && (numEigVecWord != "none")) {
		cerr << "ERROR: Invalid value of " << numEigVecWord << " for NumberOfEigenvectors(all/partial/none) " << endl;
		exit(1);
	}
	
	//Set output precision
	cout << scientific << setprecision(15) << endl;
	
	
	cout << "Eigenvalue and Eigenvector Lanczos calculator begun." << endl;
	cout << "The number of Lanczos iterations will be: " << niter << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Run Hv_prep
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	lanczosStor *lanczos_data = new lanczosStor();
	generalStor *general_data;
	//generalStor *general_data = new generalStor();
	
	void (*HvPrepPtr)(int, char**, generalStor*, lanczosStor*) = NULL;
	void (*HvPtr)(int, char**, generalStor*, lanczosStor*, double*, double*) = NULL;
	
	HvInterfaceSetup(HvCalculatorSwitch, &general_data, &HvPrepPtr, &HvPtr);
	
	//Set up what is needed for any Hv calculations, such as pre-calculating the potential and any other matricies
	(*HvPrepPtr)(argc, argv, general_data, lanczos_data);
	
	int ntotbs = lanczos_data->total_basis_size;
	string sim_descr = lanczos_data->sim_descr;
	string sim_descr_short = lanczos_data->sim_descr_short;
	
	cout << "The total basis size is: " << ntotbs << endl;
	cout << "The current simulation has been described as: " << sim_descr << endl;
	
	//test Hv
	
	cout << "Eigenvalues will be searched for in the interval (" << emin << ", " << emax << "]; all in kJ/mol." << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Create directory structure
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Make a directory to store the output files and copy over the input parameters
	
	//Get the current date and time (method from http://stackoverflow.com/questions/997946/c-get-current-time-and-date )
	time_t t = time(0); //Get current time
	struct tm *now = localtime(&t);
	char time_charArray[80];
	string time_string, commandString, dir_string;
	int retVal;
	
	strftime(time_charArray, sizeof(time_charArray), "%F_%H_%M_%S", now);
	
	time_string = time_charArray;
	dir_string = outputDirPath + sim_descr_short + "_" + time_string;
	commandString = "mkdir -p " + dir_string;
	
	retVal = system(commandString.c_str());
	
	if (retVal!=0) {
		cerr << "ERROR: Storage directory could not be created with the command: " << commandString << endl;
		exit(1);
	}
	
	cout << "The output directory is: " << dir_string << endl;
	
	//copy over input files for reference
	commandString = "cp " + lanczosInputFilename + " " + dir_string;
	system(commandString.c_str());
	
	commandString = "cp " + HvInputFilename + " " + dir_string;
	system(commandString.c_str());
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Lanczos Loop 1 - Get Tm
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	
	int ngood; //Number of good eigenvalues
	
	//Vectors for the eigenvalues and their errors
	VECT eval(niter);
	VECT evalerr(niter);
	
	//Vectors for the Lanczos algorithm Tm elements
	VECT alpha(niter);
	VECT beta(niter+1);
	VECT  beta2(niter+1);
	
	//Allocate r_j, v_j and u_j for the Lanczos Algorithm
	double* rvec=new double[ntotbs];  	
	double* vec=new double[ntotbs];
	double* uec=new double[ntotbs];	
	
	cout<<"'Lanczos Loop 1 - Get Tm' beginning."<<endl;
	
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
//#pragma omp parallel for default(shared) private (ib) schedule(guided)
		for (int ib=0;ib<ntotbs;ib++) 
			uec[ib]=0.;
		
		//u_j = u_j + H*v_j
		(*HvPtr)(argc, argv, general_data, lanczos_data, vec, uec);
		
		//r_j = r_j + u_j; ; results in r_j = (-beta_(j-1) * v_(j-1)) + (H*v_j) = H*v_j - beta_(j-1) * v_(j-1)
//#pragma omp parallel for default(shared) private (ib) schedule(guided)
		for (int ib=0;ib<ntotbs;ib++) 
			rvec[ib]+=uec[ib];
		
		//alpha_j = v_j * r_j
		alpha(j-1)=0.;
		for (int ib=0;ib<ntotbs;ib++) 
			alpha(j-1)+=vec[ib]*rvec[ib];
		
		//r_j = r_j - alpha_j * v_j
//#pragma omp parallel for default(shared) private (ib) schedule(guided)
		for (int ib=0;ib<ntotbs;ib++) 
			rvec[ib]-=(alpha(j-1)*vec[ib]);
		
		//beta_j = ||r_j||
		beta2(j)=0.;
		for (int ib=0;ib<ntotbs;ib++) 
			beta2(j)+=rvec[ib]*rvec[ib];
		
		beta(j)=sqrt(beta2(j));
		
		
		//r_j = r_j/beta_j (normalize r_j) = v_(j+1) (after swap)
//#pragma omp parallel for default(shared) private (ib) schedule(guided)
		for (int ib=0;ib<ntotbs;ib++) 
			rvec[ib]=(1./beta(j))*rvec[ib];
		
		//v_j = -beta_j * v_j = r_(j+1) (after swap)
//#pragma omp parallel for default(shared) private (ib) schedule(guided)
		for (int ib=0;ib<ntotbs;ib++) 
			vec[ib]=(-beta(j))*vec[ib]; // prepare r check minus sign!!!
		
		//Swap vectors
//#pragma omp parallel for default(shared) private (ib) schedule(guided)
		for (int ib=0;ib<ntotbs;ib++) {
			
			//Swap r_j and v_j, using u_j as storage space
			uec[ib]=vec[ib];	//u_(j+1) = v_j
			vec[ib]=rvec[ib];	//v_(j+1) = r_j; results in v_(j+1) = r_j/beta_j
			rvec[ib]=uec[ib];   //r_(j+1) = u_j = v_j; results in r_(j+1) = -beta_j * v_j
		}
		cout << j << " " << alpha(j-1) << " " << beta(j) << endl;
		/*if (j%100 == 0)
			cout<<"iteration "<<j<<endl;*/
	}                  
	
	cout << "'Lanczos Loop 1 - Get Tm' FINISHED." << endl;
	cout << "--------------------------------------------------------------------------------" << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Determine Good Eigenvalues
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	cout << endl;
	cout << "'Determination of Good Eigenvalues' beginning." << endl;
	
	//Determine Good Eigenvalues
	lancbis(niter,eval,evalerr,emin,emax,ngood,alpha,beta,beta2);
	
	cout << endl;
	cout << "Number of Good Eigenvalues is: " << ngood << endl;
	cout << "Groundstate Eigenvalue is: " << eval(0) << " kJ/mol" << endl;
	cout << "Groundstate Eigenvalue Error is: " << evalerr(0) << " kJ/mol?" << endl;
	
	// Output Eigenvalues and Errors
	ofstream lancout((dir_string + "/boundstates.txt").c_str());
	ofstream lanczpeout((dir_string + "/states_zpe.txt").c_str());
	
	lancout << scientific << setprecision(15) << endl;
	lanczpeout << scientific << setprecision(15) << endl;
	
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
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Determine the Tm Eigenvectors
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Determine how many, if any, eigenvectors are to be calculated
	
	// If only the eigenvalues are wished to be calculated, end the program here.
	if (numEigVecWord == "none") { 
			cout << "LANCZOS ALGORITHM COMPLETED - no eigenvectors calculated as based on user input." << endl;
			exit(0);
	}
	
	//If not all of the eigenvectors are desired, generate only the desired number.
	else if (numEigVecWord == "partial") {	
			ngood = numEigVec;
	} 
	
	//Last option for numEigVecWord is "all", which means I leave ngood alone
	else if (numEigVecWord == "all") {	
	}
	else {
		cerr << "Unknown value for the Number of Eigenvectors to be calculated: " << numEigVecWord << endl;
	}

	
	    
	cout << endl;
	cout << "'Eigenvector calculation' beginning." << endl;
	
	//Calculate the eigenvectors of Tm
	MAT evtr(niter,ngood);
	
	cout << "'Tm eigenvector calculation' beginning." << endl;
	
	//Calculate Tm eigenvectors and store in evtr
	lanczosvectors(alpha,beta,beta2,niter,eval,ngood,evtr); 
	
	cout << "'Tm eigenvector calculation' finished." << endl;
	
	//Output Tm eigenvectors to a file
	ofstream lvecout((dir_string + "/lv").c_str());
	VECT cumulnorm(ngood);
	double coeff2;
	
	lvecout << scientific << setprecision(15) << endl;
	
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
	
	cout << "'Lanczos Loop 2 - Get Eigenvectors of H' beginning." << endl;
	
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
	beta(niter)=0.; //beta_j and beta2_j are of length niter+1
	beta2(niter)=0.;
	
	VECT ARvL(ntotbs*ngood); //This stores the final Eigenvectors (Ritz vectors)
		
	//Reset the cumulative norm
	for (n=0;n<ngood;n++) cumulnorm(n)=0.;
	
	double treshold;
	double coeff;
	
	cout << "Lanczos Coefficients" << endl;
	cout << "Iteration j" << " " << "alpha_j" << " " << "beta_j" << endl;
	
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
		(*HvPtr)(argc, argv, general_data, lanczos_data, vec, uec);
		
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
			
			//Swap r_j and v_j, using u_j as storage space
			uec[ib]=vec[ib];	//u_(j+1) = v_j
			vec[ib]=rvec[ib];	//v_(j+1) = r_j; results in v_(j+1) = r_j/beta_j
			rvec[ib]=uec[ib];	//r_(j+1) = u_j = v_j; results in r_(j+1) = -beta_j * v_j 
		}
		
		cout << j << " " << alpha(j-1) << " " << beta(j) << endl;
		/*if (j%100 == 0)
		 cout<<"iteration "<<j<<endl;*/
	}    
	
	
	// eigenvectors for good eigenvalues are in ARvL(row+ntotbs*n)
	
	// for the ground state, eg, ARvL(row+ntotbs*0) is the vector in the original basis

	//cout<<ARvL(0)<<" "<<ARvL(ntotbs*ngood-1)<<" "<<icode[0]<<" "<<icode[numbas-1]<<endl;
	
	//Output eigenvectors to file.	
	ofstream eigVecFile((dir_string + "/eigVecFile.txt").c_str());
	
	eigVecFile << scientific << setprecision(15) << endl;
	
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
	
	cout << "'Lanczos Loop 2 - Get Eigenvectors of H' FINISHED." << endl;
	cout << "--------------------------------------------------------------------------------" << endl;
	
	cout << "Boundstate Lanczos Calculator has finished executing." << endl;
	cout << "Have an exciting day." << endl;
	
	return 0;
}


void lanczosvectors(VECT &alpha, VECT &beta, VECT &beta2, int niter, VECT &eval, int ngood, MAT &evtr) {
	int i,j;
	
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
}

