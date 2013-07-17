#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#define EPS 2.0e-16

using namespace std;

void gauleg(double x1,double x2,double *x,double *w,int n)
{
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;
	
	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++)  {
		z=cos(M_PI*((double)i-0.25)/((double)n+0.5));
		do {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*(double)j-1.0)*z*p2-((double)j-1.0)*p3)/(double)j;
			}
			pp=(double)n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		x[i-1]=xm-xl*z;
		x[n-i]=xm+xl*z;
		w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n-i]=w[i-1];
	}
}

#undef EPS

int main(int argc, char **argv) {
	double x1, x2, *x, *w;
	int i, l;
	
	x1 = atof(argv[1]);
	x2 = atof(argv[2]);
	l = atoi(argv[3]);
	
	x = new double [l];
	w = new double [l];
	
	gauleg(x1, x2, x, w, l);
	
	cout << fixed << setprecision(15);
	
	for (i=0; i<l; i++) {
		cout << "Root: " << i+1 << " " << x[i] << endl;
		cout << "Weight: " << i+1 << " " << w[i] << endl;
	}
	
	for (i=0; i<l; i++) {
		cout << x[i] << endl;
	}
}