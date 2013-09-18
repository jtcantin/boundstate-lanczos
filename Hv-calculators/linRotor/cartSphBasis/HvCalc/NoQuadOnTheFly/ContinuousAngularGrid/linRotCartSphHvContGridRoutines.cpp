#include "linRotCartSphHvContGridRoutines.h"

using namespace std;

#define m_shift(X) (X+l_max)

double* calc_Vlmlpmp_NoQuad_ContGrid(interfaceStor *interface) {
	int m, mp, n, a, b, i, j, k, np;
	int nnp, nn;
	
	//Extract out desired variables from the interface storage structure
	quadStor *gaussQuad;
	pointPotentialStorH2 *atomPotentials;
	tesseralStor *tesseralHarmonics;
	tesseralStor *tesseralHarmonics2PI;
	lmFBR *lmBasis;
	fiveDGrid *cartGrid;
	
	int l_max, **qNum, length;
	
	gaussQuad = interface->quadrature;
	atomPotentials = interface->potential;
	lmBasis = interface->lmBasis;
	cartGrid = interface->grids;
	tesseralHarmonics = interface->tesseral;
	tesseralHarmonics2PI = interface->tesseral2PI;
	
	
	l_max = lmBasis->lmax;
	qNum = lmBasis->qNum;
	length = lmBasis->length;
		
	nn = length;
	nnp = nn;
				   
	
	//Quadrature Variables
	double *wa, *wb;
	int na, nb;
	
	wa = gaussQuad->GLweights;
	na = gaussQuad->GLnum;
	
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
	
	ni = cartGrid->nx;
	
	nj = cartGrid->ny;
	
	nk = cartGrid->nz;
	
	//
	
	int matrixSize = ni*nj*nk*nn*nnp;
	int count = 0;
	int step = matrixSize/100;
	int ind, ind2, ind3;
	
	double *potentialMatrix = new double [matrixSize];
	
	cout << "Total potential matrix size: " << matrixSize << endl;
	//cout << "Note: the following statements are only an estimate of the progress." << endl;
	
	double *harmFactor = new double [nn*nnp*na*nb];
	double *harmFactorPI = new double [nn*nnp*na*nb];
	
	double *V_ijkab_1_mat = atomPotentials->fullPotential2;
	double *V_ijkab_2_mat = atomPotentials->fullPotential3;
	
#pragma omp parallel default(shared) private(i,j,k,n,np,m,mp,ind,ind2,ind3,a, b) 
	{
	//Precompute the spherical harmonics terms and the weights
#pragma omp for schedule(guided) collapse(4)
	for (n=0; n<nn; n++) {
		for (np=0; np<nnp; np++) {
			for (a=0; a<na; a++) {				
				for (b=0; b<nb; b++) {
					ind2 = ((n*nnp + np)*na + a)*nb + b;
					
					m = qNum[n][1];
					mp = qNum[np][1];
					
					harmFactor[ind2]   = wa[a] * wb[b] * L_lm1[n][a] * S_m1[m_shift(m)][b] * L_lpmp1[a][np] * S_mp1[b][m_shift(mp)];
					harmFactorPI[ind2] = wa[a] * wb[b] * L_lm2[n][a] * S_m2[m_shift(m)][b] * L_lpmp2[a][np] * S_mp2[b][m_shift(mp)];
				}
			}
		}
	}
	
	//Calculate <lm|V|l'm'>(x,y,z) = sum(a,b)[ wa * wb * ( L_lm(theta) * S_m(phi) * V(theta,phi;x,y,z) * L_lpmp(theta) * S_mp(phi) + L_lm(theta) * S_m(phi2) * V(theta,phi2;x,y,z) * L_lpmp(theta) * S_mp(phi2) )
	//      where phi = acos(cosPhiAbscissae[b])] and phi2 = 2*PI - acos(cosPhiAbscissae[b]); you need both to get the full range of phi [0,2pi)
	// NOTE: You can't collapse(7) this loop as the a and b loops are sums.
#pragma omp for schedule(guided) collapse(5)
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
//						if (count % step == 0) {
//							cout << double(count)/double(matrixSize) << endl;
//						}
//						count++;
					}
				}
			}
		}
	}
	}
	
	delete [] harmFactor;
	delete [] harmFactorPI;
	delete [] V_ijkab_1_mat;
	delete [] V_ijkab_2_mat;
	
	return potentialMatrix;
}