#include<stdio.h>

#include <math.h>
#define EPS 1.0e-4
#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5 

/* JP's declaration
double qromo(func,a,b)
double a,b;
double (*func)();
*/	
double qromo(double (*func)(double), double a, double b)
{
  /******JP's declarations
	void polint();
	void nrerror();
  	double midpnt();
  ******/
        void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
        void nrerror(char error_text[]);
	double midpnt(double (*func)(double), double a, double b, int n);
	int j;
	double zero=0.0,ss=1.0,dss=2.0,h[JMAXP+1],s[JMAXP];


	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=midpnt(func,a,b,j);
		if (j >= K) {
/*	fprintf(stderr,"QROMO %d in %f %f %d %d %f %f\n",j,
			h[j-K],s[j-K],j,K,ss,dss);
*/
			polint(&h[j-K],&s[j-K],K,zero,&ss,&dss);

			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=h[j]/9.0;
	}
	nrerror("Too many steps in routing qromo");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K
