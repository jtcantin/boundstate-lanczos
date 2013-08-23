#include "linRotCartSphHvFullVPrepRoutines.h"

using namespace std;

#define m_shift(X) (X+l_max)

double* calc_Vlmlpmp_NoQuad(interfaceStor *interface) {
	int m, mp, n, a, b, i, j, k, np;
	int nmp, nm, nnp, nn;
	double potentialCeiling;
	
	//Extract out desired variables from the interface storage structure
	quadStor *gaussQuad;
	pointPotentialStorH2 *atomPotentials;
	tesseralStor *tesseralHarmonics;
	tesseralStor *tesseralHarmonics2PI;
	lmFBR *lmBasis;
	fiveDGrid *cartGrid;
	
	int l_max, **qNum, length;
	double (*linearMoleculePotential)(interfaceStor*, H2_orient*) = NULL;
	
	gaussQuad = interface->quadrature;
	atomPotentials = interface->potential;
	lmBasis = interface->lmBasis;
	cartGrid = interface->grids;
	tesseralHarmonics = interface->tesseral;
	tesseralHarmonics2PI = interface->tesseral2PI;
	
	
	l_max = lmBasis->lmax;
	qNum = lmBasis->qNum;
	length = lmBasis->length;
	
	potentialCeiling = atomPotentials->potentialCeiling;
	linearMoleculePotential = interface->fcnPointers->linearMoleculePotential;
	
	//Get number of m values
	nmp = 2*l_max + 1;
	nm = 2*l_max + 1;
	
	nn = length;
	nnp = nn;
				   
	
	//Quadrature Variables
	double *cosThetaAbscissae, *cosPhiAbscissae, *wa, *wb, *phiAbscissae, *PIphiAbscissae, *thetaAbscissae;
	int na, nb;
	
	cosThetaAbscissae = gaussQuad->GLabscissae;
	thetaAbscissae = gaussQuad->GLacosAbscissae;
	wa = gaussQuad->GLweights;
	na = gaussQuad->GLnum;
	
	cosPhiAbscissae = gaussQuad->GCabscissae;
	phiAbscissae = gaussQuad->GCacosAbscissae;
	PIphiAbscissae = gaussQuad->GCPIacosAbscissae;
	wb = gaussQuad->GCweights;
	nb = gaussQuad->GCnum;
	//
	
	//Tesseral Harmonics Variables
	double **L_lpmp1, **S_mp1, **S_m1, **L_lm1;
	
	L_lpmp1 = tesseralHarmonics->L_lpmp;
	S_mp1 = tesseralHarmonics->S_mp;
	S_m1 = tesseralHarmonics->S_m;
	L_lm1 = tesseralHarmonics->L_lm;
	
	double **L_lpmp2, **S_mp2, **S_m2, **L_lm2;
	
	L_lpmp2 = tesseralHarmonics2PI->L_lpmp;
	S_mp2 = tesseralHarmonics2PI->S_mp;
	S_m2 = tesseralHarmonics2PI->S_m;
	L_lm2 = tesseralHarmonics2PI->L_lm;
	//
	
	//Cartesian Grid Variables
	
	int ni, nj, nk;
	double *xGrid, *yGrid, *zGrid;
	
	ni = cartGrid->nx;
	xGrid = cartGrid->x_Grid;
	
	nj = cartGrid->ny;
	yGrid = cartGrid->y_Grid;
	
	nk = cartGrid->nz;
	zGrid = cartGrid->z_Grid;
	
	//
	
	int matrixSize = ni*nj*nk*nn*nnp;
	int count = 0;
	int step = matrixSize/100;
	int ind, ind2, ind3;
	
	double *potentialMatrix = new double [matrixSize];
	double V_ijkab_1, V_ijkab_2;
	VECT CMpos(3);
	H2_orient linearMolecule;
	
	cout << "Total potential matrix size: " << matrixSize << endl;
	cout << "Note: the following statements are only an estimate of the progress." << endl;
	
	//Precompute the spherical harmonics terms and the weights
	
	double *harmFactor = new double [nn*nnp*na*nb];
	double *harmFactorPI = new double [nn*nnp*na*nb];
	
#pragma omp parallel for default(shared) private(i,j,k,n,np,m,mp,ind,ind2,ind3,linearMolecule, CMpos, a, b, V_ijkab_1, V_ijkab_2) schedule(guided) collapse(2)
	for (n=0; n<nn; n++) {
		for (np=0; np<nnp; np++) {
			
			ind2 = (n*nnp + np)*na;
			
			m = qNum[n][1];
			mp = qNum[np][1];
			
			for (a=0; a<na; a++) {				
				for (b=0; b<nb; b++) {
					harmFactor[(ind2 + a)*nb + b]   = wa[a] * wb[b] * L_lm1[n][a] * S_m1[m_shift(m)][b] * L_lpmp1[a][np] * S_mp1[b][m_shift(mp)];
					harmFactorPI[(ind2 + a)*nb + b] = wa[a] * wb[b] * L_lm2[n][a] * S_m2[m_shift(m)][b] * L_lpmp2[a][np] * S_mp2[b][m_shift(mp)];
				}
			}
		}
	}
	
	//Precompute the potential matrix	
	
	double *V_ijkab_1_mat = new double [ni*nj*nk*na*nb];
	double *V_ijkab_2_mat = new double [ni*nj*nk*na*nb];
	
#pragma omp parallel for default(shared) private(i,j,k,n,np,m,mp,ind,ind2,ind3,linearMolecule, CMpos, a, b, V_ijkab_1, V_ijkab_2) schedule(guided) collapse(3)	
	for (i=0; i<ni; i++) {
		for (j=0; j<nj; j++) {
			for (k=0; k<nk; k++) {
				ind = ((i*nj + j)*nk + k)*na;
				
				CMpos.DIM(3);
				CMpos.COOR(0) = xGrid[i];
				CMpos.COOR(1) = yGrid[j];
				CMpos.COOR(2) = zGrid[k];
				
				linearMolecule.CM = &CMpos;
				
				for (a=0; a<na; a++) {
					linearMolecule.theta = thetaAbscissae[a];
					
					for (b=0; b<nb; b++) {
						
						//Get the potential for phi in [0, pi)
						linearMolecule.phi = phiAbscissae[b]; 								
						V_ijkab_1 = (*linearMoleculePotential)(interface, &linearMolecule);
						
						if (V_ijkab_1 >= potentialCeiling) {
							V_ijkab_1 = potentialCeiling;
						}
						
						V_ijkab_1_mat[(ind + a)*nb + b] = V_ijkab_1;
						
						
						//Get the potential for phi in [pi, 2pi)								
						linearMolecule.phi = PIphiAbscissae[b];
						V_ijkab_2 = (*linearMoleculePotential)(interface, &linearMolecule);
						
						if (V_ijkab_2 >= potentialCeiling) {
							V_ijkab_2 = potentialCeiling;
						}
						
						V_ijkab_2_mat[(ind + a)*nb + b] = V_ijkab_2;
						
						
					}
				}
				
				//Give feedback as to progress
//				if (count % step == 0) {
//					cout << double(count)/double(matrixSize) << endl;
//				}
//				count++;
			}
		}
	}
	
	
	//Calculate <lm|V|l'm'>(x,y,z) = sum(a,b)[ wa * wb * ( L_lm(theta) * S_m(phi) * V(theta,phi;x,y,z) * L_lpmp(theta) * S_mp(phi) + L_lm(theta) * S_m(phi2) * V(theta,phi2;x,y,z) * L_lpmp(theta) * S_mp(phi2) )
	//      where phi = acos(cosPhiAbscissae[b])] and phi2 = 2*PI - acos(cosPhiAbscissae[b]); you need both to get the full range of phi [0,2pi)
#pragma omp parallel for default(shared) private(i,j,k,n,np,m,mp,ind,ind2,ind3,linearMolecule, CMpos, a, b, V_ijkab_1, V_ijkab_2) schedule(guided) collapse(5)
	for (i=0; i<ni; i++) {
		for (j=0; j<nj; j++) {
			for (k=0; k<nk; k++) {
				for (n=0; n<nn; n++) {
					for (np=0; np<nnp; np++) {
						
						ind = (((i*nj + j)*nk + k)*nn + n)*nnp + np;
						ind2 = (n*nnp + np)*na;
						ind3 = ((i*nj + j)*nk + k)*na;
						
						potentialMatrix[ind] = 0.0;						
					
						for (a=0; a<na; a++) {
							for (b=0; b<nb; b++) {
								potentialMatrix[ind] += harmFactor[(ind2 + a)*nb + b] * V_ijkab_1_mat[(ind3 + a)*nb + b] + harmFactorPI[(ind2 + a)*nb + b] * V_ijkab_2_mat[(ind3 + a)*nb + b];
							}
						}
						
						//Give feedback as to progress
						if (count % step == 0) {
							cout << double(count)/double(matrixSize) << endl;
						}
						count++;
					}
				}
			}
		}
	}
	
	return potentialMatrix;
}

void HvPrep_Internal_NoQuad(int argc, char **argv, interfaceStor *interface, lanczosStor *lanczos) {
	int na, nb;
	
	cout << endl;
	cout << "Hv Preparation STARTED." << endl;
	cout << "////////////////////////////////////////////////////////////////////////////" << endl;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Get required information from input file 
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	string inputFilename;
	ifstream inputFile;
	
	int nx, ny, nz, l_max, thetaPoints, phiPoints, pnx, pny, pnz;
	double x_max, y_max, z_max, px_max, py_max, pz_max;
	double rotationalConstant, totalMass, ceilingPotential;
	string geometryFilename, line, junk, simulationFilename, quadConvergeStudy, symFlag;
	
	inputFilename = argv[2];
	
	inputFile.open(inputFilename.c_str(), ios::in);
	if (inputFile.is_open()) {
		//Get rid of comment lines;
		getline(inputFile, line);
		getline(inputFile, line);
		
		//Gather data for input; each line has the name of the value separate from the value with a space
		//Always ignore the name of the value.
		inputFile >> junk;
		inputFile >> nx;		
		
		inputFile >> junk;
		inputFile >> x_max;
		
		inputFile >> junk;
		inputFile >> ny;
		
		inputFile >> junk;
		inputFile >> y_max;
		
		inputFile >> junk;
		inputFile >> nz;
		
		inputFile >> junk;
		inputFile >> z_max;
		
		inputFile >> junk;
		inputFile >> l_max;
		
		inputFile >> junk;
		inputFile >> thetaPoints;
		
		inputFile >> junk;
		inputFile >> phiPoints;
		
		inputFile >> junk;
		inputFile >> pnx;
		
		inputFile >> junk;
		inputFile >> px_max;
		
		inputFile >> junk;
		inputFile >> pny;
		
		inputFile >> junk;
		inputFile >> py_max;
		
		inputFile >> junk;
		inputFile >> pnz;
		
		inputFile >> junk;
		inputFile >> pz_max;
		
		inputFile >> junk;
		inputFile >> totalMass;
		
		inputFile >> junk;
		inputFile >> rotationalConstant;
		
		inputFile >> junk;
		inputFile >> geometryFilename;
		
		inputFile >> junk;
		inputFile >> simulationFilename;
		
		inputFile >> junk;
		inputFile >> ceilingPotential;
		
		inputFile >> junk;
		inputFile >> quadConvergeStudy;
		
		inputFile >> junk;
		inputFile >> symFlag;
		
	}
	else {
		cerr << "Input file '" << inputFilename << "' could not be opened." << endl;
		exit(1);
	}
	
	inputFile.close();
	
	//Check to ensure the geometry file can be opened before doing the rest of the calculations
	//NULL is to allow for cases where the geometry file is not needed
	inputFile.open(geometryFilename.c_str(), ios::in);
	if (!(inputFile.is_open()) && (geometryFilename != "NULL")) {
		cerr << "Atom geometry file '" << geometryFilename << "' could not be opened." << endl;
		exit(1);
	}
	inputFile.close();
	
	//Check that the system dimensions are reasonable
	if (x_max > px_max || y_max > py_max || z_max > pz_max) {
		cerr << "Error, the system's physical dimensions must be the same size or smaller than the potential universe dimensions." << endl;
		exit(1);
	}
	
	if (x_max*px_max*y_max*py_max*z_max*pz_max < DBL_EPSILON || nx*ny*nz*pnx*pny*pnz == 0) {
		cerr << "Error, none of the system dimensions nor the number of points can be zero and none should be smaller than the machine epsilon: " << DBL_EPSILON << endl;
		exit(1);
	}
	
	cout << "Hv Input file parameters read." << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Generate the x, y, and z bases and Kinetic Energy Operators
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double *xKinmat, *xGrid;
	
	cartKinGrid(x_max, nx, totalMass, &xKinmat, &xGrid); //grid[i], kinMat[(i*nx) + j]
	
	double *yKinmat, *yGrid;
	
	if ((fabs(y_max-x_max)<DBL_EPSILON) && (ny==nx)) { //If the y-grid is the same as the x-grid, reuse the matrix and save memory and time
		yKinmat = xKinmat;
		yGrid = xGrid;
	}
	else {
		cartKinGrid(y_max, ny, totalMass, &yKinmat, &yGrid); //grid[i], kinMat[(i*ny) + j]
	}
	
	double *zKinmat, *zGrid;
	
	if ((fabs(z_max-x_max)<DBL_EPSILON) && (nz==nx)) { //If the z-grid is the same as the x-grid or the y-grid, reuse the matrix and save memory and time
		zKinmat = xKinmat;
		zGrid = xGrid;
	}
	else if ((fabs(z_max-y_max)<DBL_EPSILON) && (nz==ny)) {
		zKinmat = yKinmat;
		zGrid = yGrid;
	}
	else {
		cartKinGrid(z_max, nz, totalMass, &zKinmat, &zGrid); //grid[i], kinMat[(i*nz) + j]
	}
	
	//Store everything for interface
	fiveDGrid *gridStor = new fiveDGrid();
	
	gridStor->x_Grid = xGrid;
	gridStor->nx = nx;
	gridStor->x_max = x_max;
	gridStor->xKinMat = xKinmat;
	
	gridStor->y_Grid = yGrid;
	gridStor->ny = ny;
	gridStor->y_max = y_max;
	gridStor->yKinMat = yKinmat;
	
	gridStor->z_Grid = zGrid;
	gridStor->nz = nz;
	gridStor->z_max = z_max;
	gridStor->zKinMat = zKinmat;	
	
	cout << "Cartesian grid and Kinetic Energy Operators generated." << endl;
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Generate the |lm> basis and Rotational Kinetic Energy Operator
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Generate the |lm> basis
	int **qNum, length, **index, *dims;
	dims = new int [2];
	genIndices_lm(l_max, &qNum, &length, &index, dims);
	
	//Calculate the rotational Kinetic Energy operator
	double *rotKinEngOp;
	
	rotKinEngOp = rotKinEng(qNum, length, rotationalConstant); //This vector is based on the [n] composite basis
	
	//Store everything for interface
	lmFBR *lmBasis = new lmFBR();
	
	lmBasis->lmax = l_max;
	
	lmBasis->qNum = qNum;
	lmBasis->length = length;
	
	lmBasis->index = index;
	lmBasis->dims = dims;
	
	lmBasis->rotKinMat = rotKinEngOp;
	
	cout << "|lm> basis and Rotational Kinetic Energy Operator generated." << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Pre-compute everything required for the potential
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Calculate abscissae and weights for the quadratures
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	quadStor *quadrature;
	
	quadrature = QuadraturePrep(thetaPoints, phiPoints);
	
	cout << "Gauss-Legendre and Gauss-Chebyshev quadrature weights and abscissae calculated." << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Calculate Tesseral Harmonics
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	tesseralStor *tessHarmonics, *tessHarmonics2PI;
	
	na = thetaPoints;
	nb = phiPoints;
	
	TesseralPrep(na, nb, quadrature, lmBasis, &tessHarmonics, &tessHarmonics2PI);
		
	cout << "Tesseral Harmonics terms calculated." << endl;
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Pre-Calculate the Potential
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	cout << "Pre-calculating the potential:" << endl;
	cout << "----------------------------------------------------------------------------" << endl;
	
	int numDim = 3; //There are 3 spatial dimensions
	double *gridMax = new double [3];
	int *gridPoints = new int [3];
	
	gridMax[0] = px_max;
	gridPoints[0] = pnx;
	
	gridMax[1] = py_max;
	gridPoints[1] = pny;
	
	gridMax[2] = pz_max;
	gridPoints[2] = pnz;
	
	
	pointPotentialStorH2 *partialPotential;
	
	partialPotential = (*(interface->fcnPointers->preCalcPotential))(numDim, gridMax, gridPoints, geometryFilename, interface, argc, argv);
	
	partialPotential->potentialCeiling = ceilingPotential;
	
	//Put everything in the interface sructure
	interface->grids = gridStor;
	interface->lmBasis = lmBasis;
	interface->quadrature = quadrature;
	interface->tesseral = tessHarmonics;
	interface->tesseral2PI = tessHarmonics2PI;
	interface->potential = partialPotential;
	
	//Put what is needed for the Lanczos algorithm in lanczos	
	lanczos->total_basis_size = nx*ny*nz*length;
	lanczos->sim_descr = simulationFilename; //Keep it at this for now.
	lanczos->sim_descr_short = simulationFilename;
	
	//Check if a convergence study is desired
	if (quadConvergeStudy == "TRUE") {
		quadratureConvergenceStudy_NoQuad(interface, lanczos);
	}
	else {
		cout << "No n_theta, n_phi convergence study performed, using values in : " << inputFilename << endl;
	}
	
	//Pre-Calculate <lm|V|lpmp>
	interface->potential->fullPotential = calc_Vlmlpmp_NoQuad(interface);
	
	cout << "Potentials Pre-calculated." << endl;
	
	interface->lmBasis->symmeterizer = calcSym(interface, symFlag);
	

	cout << "////////////////////////////////////////////////////////////////////////////" << endl;
	
	//Free up memory for later use
	delete [] interface->potential->CMpotential;
	delete [] interface->potential->H_potential;
	delete interface->potential->potentialUniverse;
	
	cout << "Hv Preparation FINISHED." << endl;
}

double* calcSym(interfaceStor *interface, string symFlag) {
	int nn = interface->lmBasis->length;
	int **qNum = interface->lmBasis->qNum;
	int n, l;
	
	double *symmeterizer = new double [nn];
	
	if (symFlag == "All") {
		for (n=0; n<nn; n++) {
			symmeterizer[n] = 1.0;
		}
	}
	else if (symFlag == "Even") {
		for (n=0; n<nn; n++) {
			l = qNum[n][0];
			
			if (l % 2 == 0) {
				symmeterizer[n] = 1.0;
			}
			else {
				symmeterizer[n] = 0.0;
			}
		}
	}
	else if (symFlag == "Odd") {
		for (n=0; n<nn; n++) {
			l = qNum[n][0];
			
			if (l % 2 == 1) {
				symmeterizer[n] = 1.0;
			}
			else {
				symmeterizer[n] = 0.0;
			}
		}
	}
	else {
		cerr << "ERROR: Symmeterizer Flag not recognized." << endl;
		exit(1);
	}

	
	return symmeterizer;
}
	

void quadratureConvergenceStudy_NoQuad(interfaceStor *interface, lanczosStor *lanczos) {

	int i,a,b,j;
	int basisSize;
	int na, nb;
	
	int nx, ny, nz, l_max, npx, npy, npz;
	
	basisSize = lanczos->total_basis_size;
	
	na = interface->quadrature->GLnum;
	nb = interface->quadrature->GCnum;
	
	nx = interface->grids->nx;
	ny = interface->grids->ny;
	nz = interface->grids->nz;
	
	l_max = interface->lmBasis->lmax;
	
	npx = interface->potential->potentialUniverse->grid_num[0];
	npy = interface->potential->potentialUniverse->grid_num[1];
	npz = interface->potential->potentialUniverse->grid_num[2];
	
	
	double* vec=new double[basisSize];
	double* uec;
	double expectVal = 0.0;
	
	quadStor *quadrature;
	tesseralStor *tessHarmonics, *tessHarmonics2PI;
	
	string outputFilename = lanczos->sim_descr_short;
	
	ofstream outputFile(outputFilename.c_str());

	cout << "Quadrature Convergence Study beginning." << endl;
	
	for (i=0; i<basisSize; i++) {
		vec[i] = 1.0/sqrt((double)basisSize);
	}
	
	outputFile << "#nx= " << nx << " ny= " << ny << " nz= " << nz;
	outputFile << " l_max= " << l_max << " nThetaMax= " << na << " nPhiMax= " << nb;
	outputFile << "npx= " << npx << " npy= " << npy << " npz= " << npz;
	outputFile << endl;
	
	outputFile << "#ThetaPoints" << " " << "PhiPoints" << " " << "<v|H|v>" << endl;
	
	j = 0;
	
	for (a=1; a<=na; a++) {
		for (b=1; b<=nb; b++) {
			j++;
			
			quadrature = QuadraturePrep(a, b);
			
			TesseralPrep(a, b, quadrature, interface->lmBasis, &tessHarmonics, &tessHarmonics2PI);
			
			delete interface->quadrature;
			delete interface->tesseral;
			delete interface->tesseral2PI;
			
			interface->quadrature = quadrature;
			interface->tesseral = tessHarmonics;
			interface->tesseral2PI = tessHarmonics2PI;
			
			delete interface->potential->fullPotential;
			
			//Pre-Calculate <lm|V|lpmp>
			interface->potential->fullPotential = calc_Vlmlpmp_NoQuad(interface);
			
			uec = Hv_5D_oneCompositeIndex_NoQuad(interface, vec);
			
			for (i=0; i<basisSize; i++) {
				expectVal += vec[i]*uec[i];
			}
			
			outputFile << a << " " << b << " " << expectVal << endl;
			//outputFile << a << " " << expectVal << endl;
			
			if (j % 10 == 0) {
				cout << "Working... at " << j << " out of " << na*nb << endl;
			}
			
			expectVal = 0.0;
			
		}
	}
	
	
	
}

// Calculate |u> = V|v>; <xyzlm|u> = sum(l'm')[ <xyzlm|V|xyzl'm'><xyzl'm'|v> ]; V is a matrix diagonal in x,y,z with <lm|V|l'm'>(x,y,z) on the diagonal
double* Vv_5D_oneCompositeIndex_NoQuad(interfaceStor *interface, double *v_ijknp) {
	int i, j, k, n, np;
	int ni, nj, nk, nn, nnp;
	
	ni = interface->grids->nx;
	nj = interface->grids->ny;
	nk = interface->grids->nz;
	
	nn = interface->lmBasis->length;
	nnp = nn;
	
	double *u_ijkn = new double [ni*nj*nk*nn];
	double *V_npnkji;
	
	V_npnkji = interface->potential->fullPotential;
	
#pragma omp parallel for default(shared) private(i,j,k,n,np) schedule(guided) collapse(4)
	for (i=0; i<ni; i++) {
		for (j=0; j<nj; j++) {
			for (k=0; k<nk; k++) {
				for (n=0; n<nn; n++) {
					u_ijkn[((n*nk + k)*nj + j)*ni + i] = 0.0;
					
					for (np=0; np<nnp; np++) {
						u_ijkn[((n*nk + k)*nj + j)*ni + i] += V_npnkji[(((i*nj + j)*nk + k)*nn + n)*nnp + np] * v_ijknp[((np*nk + k)*nj + j)*ni + i];
					}
				}
			}
		}
	}
	
	return u_ijkn;
}

double* Hv_5D_oneCompositeIndex_NoQuad(interfaceStor *interface, double *v_ipjkn) {
	int i, j, k, n;
	
	int ni, nj, nk, nn, basis_size;	
	ni = interface->grids->nx;
	nj = interface->grids->ny;
	nk = interface->grids->nz;
	
	nn = interface->lmBasis->length;
	
	basis_size = ni*nj*nk*nn;	
	
	double *Tv_ijkn;
	double *Vv_ijkn;
	double *S_n;
	double *Hv_ijkn = new double [basis_size];
	
	//Calculate the kinetic energy terms
	Tv_ijkn = Tv_5D_oneCompositeIndex(interface, v_ipjkn);
	
	//cout << "Tv finished" << endl;
	
	//Calculate the potential energy terms
	Vv_ijkn = Vv_5D_oneCompositeIndex_NoQuad(interface, v_ipjkn);
	
	//cout << "Vv finished" << endl;
	
	//Get symmeterizer
	S_n = interface->lmBasis->symmeterizer;
	
	//Sum Tv_ijkn and Vv_ijkn to get Hv_ijkn
	/*//#pragma omp parallel for default(shared) private(p) schedule(guided) //Parallelization justified as no processor will access the same memory location at the same time (i.e. p is different for each processor)
	for (p=0; p<basis_size; p++) {
		Hv_ijkn[p] = Tv_ijkn[p] + Vv_ijkn[p]; 
	//	Hv_ijkn[p] = Tv_ijkn[p]; // Free system - only for debugging purposes
	}*/
	
	//Sum Tv_ijkn and Vv_ijkn to get Hv_ijkn and then multiply by S_n to get all, even or odd l
	for (i=0; i<ni; i++) {
		for (j=0; j<nj; j++) {
			for (k=0; k<nk; k++) {
				for (n=0; n<nn; n++) {
					Hv_ijkn[((n*nk + k)*nj + j)*ni + i] = Tv_ijkn[((n*nk + k)*nj + j)*ni + i] + Vv_ijkn[((n*nk + k)*nj + j)*ni + i];
					Hv_ijkn[((n*nk + k)*nj + j)*ni + i] *= S_n[n];
				}
			}
		}
	}
	
	delete [] Tv_ijkn;
	delete [] Vv_ijkn;
	
	return Hv_ijkn;
}

void Hv_Prep_linRotCartSph_NoQuad(int argc, char **argv, generalStor *general_data, lanczosStor *lanczos_data) {
	
	interfaceStor *Hv_data;
	Hv_data = reinterpret_cast<interfaceStor*> (general_data);
	
	HvPrep_Internal_NoQuad(argc, argv, Hv_data, lanczos_data);
	
};

void Hv_linRotCartSph_NoQuad(int argc, char **argv, generalStor *general_data, lanczosStor *lanczos_data, double *vec, double *uec) {
	int i;
	double *uec1;
	
	interfaceStor *Hv_data;
	Hv_data = reinterpret_cast<interfaceStor*> (general_data);
	
	uec1 = Hv_5D_oneCompositeIndex_NoQuad(Hv_data, vec);
	
	for (i=0; i<lanczos_data->total_basis_size; i++) {
		uec[i] += uec1[i];
	}
	
	delete [] uec1;
	
};
