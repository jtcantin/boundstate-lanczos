
#include "AlaviHvRoutines.h"

using namespace std;

#define m_shift(X) (X+l_max)




//This calculates the vector multiplied by the potential, for a given xyz and either acos(cosPhiAbscissae) or 2PI - acos(cosPhiAbscissae)
//Test this with both a potential = 1 and potential = cos(theta) or cos(phi)
double* calc_ulm(double x, double y, double z, double *v_lpmp, interfaceStor *interface, int rangeFlag) {
	int m, mp, n, a, b;
	int anmp, nmp, nm, bna, anb, mna;
	double potentialCeiling;
	
	//Extract out desired variables from the interface storage structure
	quadStor *gaussQuad;
	pointPotentialStorH2 *atomPotentials;
	tesseralStor *tesseralHarmonics;
	lmFBR *lmBasis;
	
	int l_max, **qNum, length;
	
	gaussQuad = interface->quadrature;
	atomPotentials = interface->potential;
	lmBasis = interface->lmBasis;
	
	l_max = lmBasis->lmax;
	qNum = lmBasis->qNum;
	length = lmBasis->length;
	
	potentialCeiling = atomPotentials->potentialCeiling;
	
	//Determine whether I am calculating phi = acos(cosPhiAbscissae) or phi = 2PI - acos(cosPhiAbscissae)
	if (rangeFlag == 0) {
		tesseralHarmonics = interface->tesseral;
	}
	else {
		tesseralHarmonics = interface->tesseral2PI;
	}
	
	//Get number of m values
	nmp = 2*l_max + 1;
	nm = 2*l_max + 1;
	
	//Quadrature Variables
	double *cosThetaAbscissae, *cosPhiAbscissae, *wa, *wb;
	int na, nb;
	
	cosThetaAbscissae = gaussQuad->GLabscissae;
	wa = gaussQuad->GLweights;
	na = gaussQuad->GLnum;
	
	cosPhiAbscissae = gaussQuad->GCabscissae;
	wb = gaussQuad->GCweights;
	nb = gaussQuad->GCnum;
	//
	
	//Tesseral Harmonics Variables
	double **L_lpmp, **S_mp, **S_m, **L_lm;
	
	L_lpmp = tesseralHarmonics->L_lpmp;
	S_mp = tesseralHarmonics->S_mp;
	S_m = tesseralHarmonics->S_m;
	L_lm = tesseralHarmonics->L_lm;
	//
	
	//Memory allocation and assignment
	//Loop 1
	double *vt_mpa;
	vt_mpa = new double [na*nmp];
	
	//Loop 2
	double	*u_ab;
	u_ab = new double [na*nb];
	
	//Loop 3
	double *ut_ab;
	ut_ab = new double [na*nb];
	VECT CMpos(3);
	H2_orient linearMolecule;	
	double V_ab;
	
	//Loop 4
	double *ut_ma;
	ut_ma = new double [nm*na];
	
	//Loop 5
	double *u_lm;
	u_lm = new double [length];
	
#pragma omp parallel default(shared) private (a,mp,anmp,n,b,bna,m,mna,anb,linearMolecule,CMpos,V_ab)
	{
	
	//Loop 1 - vt_mpa = L_lpmp(q_a) * v_lpmp; FLOPS = na * length = na * (l_max+1)^2
#pragma omp for schedule(guided) collapse(2)
	for (a = 0; a<na; a++) {
		for (mp=-l_max; mp<=l_max; mp++) {
			anmp = a * nmp;			
			vt_mpa[anmp + m_shift(mp)] = 0.0; //Initialize vector
			
			for (n=0; n<length; n++) {
				vt_mpa[anmp + m_shift(mp)] += L_lpmp[a][n] * v_lpmp[n];
			}
		}
	}
	
	//Loop 2 - u_ab = S_mp(Pb) * vt_mpa; FLOPS = na * nb * nm = na * nb * (2l_max+1)
#pragma omp for schedule(guided) collapse(2)
	for (b=0; b<nb; b++) {
		for  (a=0; a<na; a++){
			bna = b * na;
			anmp = a * nmp;			
			u_ab[bna + a] = 0.0;
			
			for (mp=-l_max; mp<=l_max; mp++) {
				u_ab[bna + a] += S_mp[b][m_shift(mp)] * vt_mpa[anmp + m_shift(mp)]; //Pb------------------------------
			}
		}
	}
	
	
	//Loop 3 - ut_ab = V_ab(xyz, Qa, Pb) * u_ab; FLOPS = na * nb
	CMpos.DIM(3);
	CMpos.COOR(0) = x;
	CMpos.COOR(1) = y;
	CMpos.COOR(2) = z;
	
	linearMolecule.CM = &CMpos;
	
#pragma omp for schedule(guided) collapse(2)
	for (b=0; b<nb; b++) {
		for (a=0; a<na; a++) {
			bna = b * na;
			anb = a * nb;
			
			linearMolecule.theta = acos(cosThetaAbscissae[a]);
			//Determine whether I am calculating phi = acos(cosPhiAbscissae) or phi = 2PI - acos(cosPhiAbscissae) //Pb------------------------------
			if (rangeFlag == 0) { //This if statement may be SLOW!!!!!!!!!!!!!
				linearMolecule.phi = acos(cosPhiAbscissae[b]); 
			}
			else {
				linearMolecule.phi = 2*PI - acos(cosPhiAbscissae[b]);
			}
			
			//Calculate potential at x, y, z, theta, phi EDIT
			V_ab = Alavi_H2_Eng_Point(interface, &linearMolecule);
			
			//V_ab = 0.0; //Set to zero for debugging purposes
			
			if (V_ab >= potentialCeiling) {
				V_ab = potentialCeiling;
			}
			
			ut_ab[anb + b] = V_ab * u_ab[bna + a];
		}
	}
	
	//Loop 4 - ut_ma = wb * S_m(Pb) * ut_ab; FLOPS = na * nb * nm = na * nb * (2l_max+1)
#pragma omp for schedule(guided) collapse(2)
	for (m=-l_max; m<=l_max; m++) {
		for (a=0; a<na; a++) {
			mna = m_shift(m) * na;
			ut_ma[mna + a] = 0.0;
			anb = a * nb;
			
			for (b=0; b<nb; b++) {
				ut_ma[mna + a] += wb[b] * S_m[m_shift(m)][b] * ut_ab[anb + b]; //Pb------------------------------
			}
		}
	}
	
	//Loop 5 - u_lm = wa * L_lm(Qa) * ut_ma; FLOPS = na * length = na * (l_max+1)^2
#pragma omp for schedule(guided)
	for (n=0; n<length; n++) {
		u_lm[n] = 0.0;
		
		m = qNum[n][1];
		mna = m_shift(m) * na;
		
		for (a=0; a<na; a++) {
			u_lm[n] = wa[a] * L_lm[n][a] * ut_ma[mna + a];
		}
	}
	}
	
	delete [] ut_ma;
	delete [] u_ab;
	delete [] ut_ab;
	delete [] vt_mpa;
	
	return u_lm;
	
}

void HvPrep_Internal(int argc, char **argv, interfaceStor *interface, lanczosStor *lanczos) {
	int a, b, na, nb, nm, n, m, l;
	
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
	string geometryFilename, line, junk, simulationFilename;
	
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
	}
	else {
		cerr << "Input file '" << inputFilename << "' could not be opened." << endl;
		exit(1);
	}
	
	inputFile.close();
	
	//Check to ensure the geometry file can be opened before doing the rest of the calculations
	inputFile.open(geometryFilename.c_str(), ios::in);
	if (!(inputFile.is_open())) {
		cerr << "Atom geometry file '" << geometryFilename << "' could not be opened." << endl;
		exit(1);
	}
	inputFile.close();
	
	//Check that the system dimensions are reasonable
	if (x_max > px_max || y_max > py_max || z_max > pz_max) {
		cerr << "Error, the system dimensions must be the same size or smaller than the potential universe dimensions." << endl;
		exit(1);
	}
	
	if (x_max*px_max*y_max*py_max*z_max*pz_max < DBL_EPSILON || nx*ny*nz*pnx*pny*pnz == 0) {
		cerr << "Error, none of the system dimensions nor the number of points can be zero and none should be smaller than the machine epsilon: " << DBL_EPSILON << endl;
		exit(1);
	}
	
	cout << "Input file parameters read." << endl;
	
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
	
	nm = 2*l_max + 1;
	
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
	//Gauss-Legendre
	double minVal, maxVal, *cosThetaAbscissae, *cosThetaWeights;
	minVal = -1.0;
	maxVal = 1.0;
	
	cosThetaAbscissae = new double [thetaPoints];
	cosThetaWeights = new double [thetaPoints];
	na = thetaPoints;
	nb = phiPoints;
	
	gauleg(minVal, maxVal, cosThetaAbscissae, cosThetaWeights, thetaPoints);
	
	//Gauss-Chebyshev
	double *cosPhiAbscissae, *cosPhiWeights;
	
	gaussChebyshev(phiPoints, &cosPhiAbscissae, &cosPhiWeights);
	
	//Store everything for interface
	quadStor *quadrature = new quadStor();
	
	quadrature->GCabscissae = cosPhiAbscissae;
	quadrature->GCweights = cosPhiWeights;
	quadrature->GCnum = phiPoints;
	
	quadrature->GLabscissae = cosThetaAbscissae;
	quadrature->GLweights = cosThetaWeights;
	quadrature->GLnum = thetaPoints;	
	
	cout << "Gauss-Legendre and Gauss-Chebyshev quadrature weights and abscissae calculated." << endl;
	
	//Calculate the Tesseral Harmonics terms and rearrange appropriately for phi = acos(cosPhiAbscissae)
	tesseralStor *tessHarmonics = new tesseralStor();
	double *stor;
	
	tessHarmonics->na = na;
	tessHarmonics->nb = nb;
	tessHarmonics->lmax = l_max;
	
	//Memory allocation
	
	//Sec 1	
	tessHarmonics->L_lpmp = new double* [na];
	
	//Sec 2
	tessHarmonics->S_mp = new double* [nb];
	
	for (b=0; b<nb; b++) {
		tessHarmonics->S_mp[b] = new double [nm];
	}
	
	//Sec3
	tessHarmonics->L_lm = new double* [length];
	
	for (n=0; n<length; n++) {
		tessHarmonics->L_lm[n] = new double [na];
	}
	
	//Sec 4
	tessHarmonics->S_m = new double* [nm];
	
	for (m=-l_max; m<=l_max; m++) {
		tessHarmonics->S_m[m_shift(m)] = new double [nb];
	}
	
	//Sec 5
	tesseralStor *tessHarmonics2PI = new tesseralStor();
	
	tessHarmonics2PI->na = na;
	tessHarmonics2PI->nb = nb;
	tessHarmonics2PI->lmax = l_max;
	
	tessHarmonics2PI->L_lpmp = new double* [na];
	
	//Sec 6
	tessHarmonics2PI->S_mp = new double* [nb];
	
	for (b=0; b<nb; b++) {
		tessHarmonics2PI->S_mp[b] = new double [nm];
	}
	
	//Sec 7
	tessHarmonics2PI->L_lm = new double* [length];
	
	for (n=0; n<length; n++) {
		tessHarmonics2PI->L_lm[n] = new double [na];
	}
	
	//Sec 8
	tessHarmonics2PI->S_m = new double* [nm];
	
	for (m=-l_max; m<=l_max; m++) {
		tessHarmonics2PI->S_m[m_shift(m)] = new double [nb];
	}
	
	
#pragma omp parallel default(shared) private(a,b,m,l,n,stor)
	{
	
	//Sec 1	
#pragma omp for schedule(guided)
	for (a=0; a<na; a++) {
		tessHarmonics->L_lpmp[a] = normAssocLegendrePoly(qNum, length, cosThetaAbscissae[a]);
	}
	
	//Sec 2
#pragma omp for schedule(guided)
	for (b=0; b<nb; b++) {		
		stor = tesseralTrigTerm(qNum, length, acos(cosPhiAbscissae[b])); //Phi --------------------------------------------------------
		
		for (m=-l_max; m<=l_max; m++) { //Reshift indices around
			l=l_max;
			n = index[l][m_shift(m)];
			
			if (n<0) {
				cerr << "Error in finding the index within HvPrep_Internal for S_mp" << endl;
				exit(1);
			}
			
			tessHarmonics->S_mp[b][m_shift(m)] = stor[n];			
		}
		
		delete [] stor;
	}
	
	//Sec3
#pragma omp for schedule(guided)
	for (a=0; a<na; a++) {
		stor = normAssocLegendrePoly(qNum, length, cosThetaAbscissae[a]);
		
		for (n=0; n<length; n++) {
			tessHarmonics->L_lm[n][a] = stor[n];
		}
		delete [] stor;
	}
	
	//Sec 4	
#pragma omp for schedule(guided)
	for (b=0; b<nb; b++) {
		stor = tesseralTrigTerm(qNum, length, acos(cosPhiAbscissae[b])); //Phi --------------------------------------------------------
		
		for (m=-l_max; m<=l_max; m++) {
			l=l_max;
			n = index[l][m_shift(m)];
			
			if (n<0) {
				cerr << "Error in finding the index within HvPrep_Internal for S_m" << endl;
				exit(1);
			}
			
			tessHarmonics->S_m[m_shift(m)][b] = stor[n];
		}
		delete [] stor;
	}
	
	//Calculate the Tesseral Harmonics terms and rearrange appropriately for phi = 2PI - acos(cosPhiAbscissae)
	
	//Sec 5	
#pragma omp for schedule(guided)
	for (a=0; a<na; a++) {
		tessHarmonics2PI->L_lpmp[a] = normAssocLegendrePoly(qNum, length, cosThetaAbscissae[a]);
	}
	
	//Sec 6	
#pragma omp for schedule(guided)
	for (b=0; b<nb; b++) {		
		stor = tesseralTrigTerm(qNum, length, 2.0*PI - acos(cosPhiAbscissae[b])); //Phi --------------------------------------------------------
		
		for (m=-l_max; m<=l_max; m++) { //Reshift indices around
			l=l_max;
			n = index[l][m_shift(m)];
			
			if (n<0) {
				cerr << "Error in finding the index within HvPrep_Internal for S_mp" << endl;
				exit(1);
			}
			
			tessHarmonics2PI->S_mp[b][m_shift(m)] = stor[n];			
		}
		delete [] stor;
	}
	
	//Sec 7
#pragma omp for schedule(guided)
	for (a=0; a<na; a++) {
		stor = normAssocLegendrePoly(qNum, length, cosThetaAbscissae[a]);
		
		for (n=0; n<length; n++) {
			tessHarmonics2PI->L_lm[n][a] = stor[n];
		}
		delete [] stor;
	}
	
	//Sec 8
#pragma omp for schedule(guided)
	for (b=0; b<nb; b++) {
		stor = tesseralTrigTerm(qNum, length, 2.0*PI - acos(cosPhiAbscissae[b])); //Phi --------------------------------------------------------
		
		for (m=-l_max; m<=l_max; m++) {
			l=l_max;
			n = index[l][m_shift(m)];
			
			if (n<0) {
				cerr << "Error in finding the index within HvPrep_Internal for S_m" << endl;
				exit(1);
			}
			
			tessHarmonics2PI->S_m[m_shift(m)][b] = stor[n];
		}
		delete [] stor;
		
	}
	}
	
	//Everything already stored for interface in tessHarmonics2PI and tessHarmonics
	
	
	cout << "Tesseral Harmonics terms calculated." << endl;
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Pre-Calculate the Potential EDIT
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	cout << "Pre-calculating the partial potentials:" << endl;
	cout << "----------------------------------------------------------------------------" << endl;
	
	//Get the grid for the potential
	universeProp *potentialUniverse;
	int numDim = 3; //There are 3 spatial dimensions
	double *gridMax = new double [3];
	int *gridPoints = new int [3];
	
	gridMax[0] = px_max;
	gridPoints[0] = pnx;
	
	gridMax[1] = py_max;
	gridPoints[1] = pny;
	
	gridMax[2] = pz_max;
	gridPoints[2] = pnz;
	
	potentialUniverse = generateGrid(numDim, gridMax, gridPoints);
	
	//Get the system geometry for TIP4P molecules
	sysAtoms *atomGeo = new sysAtoms();
	
	getTIP4Patoms(&(atomGeo->atomType), &(atomGeo->atomPos), &(atomGeo->nAtoms), geometryFilename);
	
	//Get the partial potentials
	double *CMpotential, *Hpotential;
	
	CMpotential = new double [potentialUniverse->sysSize];
	Hpotential = new double [potentialUniverse->sysSize];
	
	Alavi_TIP4P_point_Eng(CMpotential, Hpotential, potentialUniverse, atomGeo);
	
	cout << "----------------------------------------------------------------------------" << endl;
	
	//Store everything for interface
	pointPotentialStorH2 *partialPotential = new pointPotentialStorH2();
	
	partialPotential->CMpotential = CMpotential;
	partialPotential->H_potential = Hpotential;
	partialPotential->potentialUniverse = potentialUniverse;
	partialPotential->potentialCeiling = ceilingPotential;
	
	//delete [] atomGeo->atomType;
	//	delete [] atomGeo->atomPos;
	delete atomGeo;
	
	cout << "Partial potentials calculated." << endl;
	
	cout << "////////////////////////////////////////////////////////////////////////////" << endl;
	
	
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
	
	cout << "Hv Preparation FINISHED." << endl;
}
