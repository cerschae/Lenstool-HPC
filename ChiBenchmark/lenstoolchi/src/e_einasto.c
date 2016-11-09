/*fichier source pour le profil de Einasto
ROBERT Frédéric 17 Juin 2013*/

/****************************************************************/
/*              name:           e_einasto                       */
/*              author:         Robert Frédéric                 */
/*              date:              06/2013                      */
/*              place:         Lyon, CRAL                       */
/****************************************************************/


#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<stdlib.h>
#include<string.h>


#define  ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

static double d(double n);

//einasto
#define CMAX 20
#define LMAX 80
extern float Tab1[LMAX][CMAX];
extern float Tab2[LMAX][CMAX];
extern float Tab3[LMAX][CMAX];

// From Numerical recipes
static double gammln(double xx);
static double gammp(double a, double x);
static void gser(double *gamser, double a, double x, double *gln);
static void gcf(double *gammcf, double a, double x, double *gln);


static double d(double n)
{
    return( 3.*n-1./3 + 8./1215./n + 184./229635./n/n/ + 1048./31000725./n/n/n - 17557576./1242974068875./n/n/n/n);
}

float einasto_masse( double r, double rs, int n, double rhos)
{
    /* D'après l'article analytical properties of einasto dark matter haloes */

    	float r_2=(rs*pow(2*n,n))/pow(d(n),n);
    	float rho_2=(rhos*exp(d(n)))/exp(2*n);
	float h=r_2/pow(2*n,n);
	float rho_0=rho_2*exp(2*n);
	float x=r/h;
	float rapport=pow(r/r_2,2);
	float log_rapport=log10(rapport);
	int val_a_chercher=floor((log_rapport+5)/0.1);

	float masse=(sqrt(n)*rho_0*pow(h,3))/(2*pow((2*M_PI),n-2))*pow(x,3)*Tab1[val_a_chercher][2*n-1];

	return masse;
}

float einasto_sigma( double r, double rs, int n, double rhos)
{
    float r_2=(rs*pow(2*n,n))/pow(d(n),n);
    float rho_2=(rhos*exp(d(n)))/exp(2*n);
    float h=r_2/pow(2*n,n);
    float rho_0=rho_2*exp(2*n);
    float x=r/h;
    float rapport=pow(r/r_2,2);
    float log_rapport=log10(rapport);
    int val_a_chercher=floor((log_rapport+5)/0.1);

    float sigma=(sqrt(n)*rho_0*h)/(pow((2*M_PI),n-1))*x*Tab2[val_a_chercher][2*n-1];
    return sigma;
}




float einasto_kappa( double r, double rs, int n, double rhos, double kappa_crit)
{
   	float r_2=(rs*pow(2*n,n))/pow(d(n),n);
    	float rho_2=(rhos*exp(d(n)))/exp(2*n);
    	float h=r_2/pow(2*n,n);
    	float rho_0=rho_2*exp(2*n);
    	float x=r/h;
    	float rapport=pow(r/r_2,2);
    	float log_rapport=log10(rapport);
    	int val_a_chercher=floor((log_rapport+5)/0.1);

	float kappa_x=((kappa_crit)/(2*pow(2*M_PI,n-1)*sqrt(n)*exp(gammln(n))))*x*Tab2[val_a_chercher][2*n-1];
	return kappa_x;
}

float einasto_kappa_av( double r, double rs, int n, double rhos, double kappa_crit)
{
	float r_2=(rs*pow(2*n,n))/pow(d(n),n);
    	float rho_2=(rhos*exp(d(n)))/exp(2*n);
    	float h=r_2/pow(2*n,n);
    	float rho_0=rho_2*exp(2*n);
    	float x=r/h;
    	float rapport=pow(r/r_2,2);
    	float log_rapport=log10(rapport);
    	int val_a_chercher=floor((log_rapport+5)/0.1);

	float kappa_av=((kappa_crit)/(2*pow(2*M_PI,n-1)*sqrt(n)*exp(gammln(n))))*x*Tab1[val_a_chercher][2*n-1];
	return kappa_av;
}

double einasto_gamma( double r, double rs, int n, double rhos, double kappa_crit)
{
	return(einasto_kappa_av(r,rs,n,rhos,kappa_crit)-einasto_kappa(r,rs,n,rhos,kappa_crit));
}

double einasto_kappa_eps(double r, double rs, double n, double theta, double kappa_crit, double rhos, double eps)
{
	double kappa_eps;

	kappa_eps=einasto_kappa(r,rs,n,rhos, kappa_crit)+eps*cos(2.*theta)*einasto_gamma(r,rs,n,rhos,kappa_crit);

	return(kappa_eps);
}

struct point  einasto_gamma_eps(double r, double rs, double n, double theta, double kappa_crit, double rhos,double eps)
{
	    struct point gamma_eps;
	        double kappa, gamma;

		    kappa = einasto_kappa(r, rs, n,rhos, kappa_crit);
		        gamma = einasto_gamma(r,rs,n,rhos,kappa_crit);

//			    gamma_eps = sqrt(pow(einasto_gamma(r, rs, n,rhos, kappa_crit), 2.) + 2.*eps * cos(2.*theta) * einasto_kappa(r, rs, n,rhos, kappa_crit) * einasto_gamma(r, rs, n,rhos, kappa_crit) + eps * eps * (pow(einasto_kappa(r, rs, n,rhos, kappa_crit), 2.) - pow(sin(2.*theta) * einasto_gamma(r, rs, n,rhos, kappae), 2.)));
			        			    
			            gamma_eps.x = gamma * cos(2. * theta) + eps * kappa;  // gamma1 
			                gamma_eps.y = -sqrt(1. - eps * eps) * gamma * sin(2. * theta); // gamma2
			    
			                    return(gamma_eps);
			                    }
			    
float einasto_phi(double r, double rs, int n, double rhos, double kappa_crit)
{
	float r_2=(rs*pow(2*n,n))/pow(d(n),n);
    	float rho_2=(rhos*exp(d(n)))/exp(2*n);
    	float h=r_2/pow(2*n,n);
    	float rho_0=rho_2*exp(2*n);
    	float x=r/h;
    	float rapport=pow(r/r_2,2);
    	float log_rapport=log10(rapport);
    	int val_a_chercher=floor((log_rapport+5)/0.1);

	float phi=((kappa_crit)/(4*pow(2*M_PI,n-1)*sqrt(n)*exp(gammln(n))))*pow(x,3)*Tab3[val_a_chercher][2*n-1];
	return phi;
}

float einasto_alpha( double r, double rs, int n, double rhos, double kappa_crit)
{
	float r_2=(rs*pow(2*n,n))/pow(d(n),n);
    	float rho_2=(rhos*exp(d(n)))/exp(2*n);
    	float h=r_2/pow(2*n,n);
    	float rho_0=rho_2*exp(2*n);
    	float x=r/h;
    	float rapport=pow(r/r_2,2);
    	float log_rapport=log10(rapport);
    	int val_a_chercher=floor((log_rapport+5)/0.1);

	float alpha=((kappa_crit)/(2*pow(2*M_PI,n-1)*sqrt(n)*exp(gammln(n))))*pow(x,2)*Tab1[val_a_chercher][2*n-1];
	return alpha;
}
		

/****************** Partie numerical recipes in C *****************/

static double gammln(double xx)
{
	double x,y,tmp,ser;
	double cof[6]={ 76.18009172947146,
			   -86.50532032941677,
			   24.01409824083091,
			   -1.231739572450155,
			   0.1208650973866179e-2,
			   -0.5395239384953e-5};
	int j;
	
	y=x=xx;
	tmp =x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}


/* --------------------------------------------------------------*/
/* Pointers (Global variables?) used : 
 * - *gammcf, *gamser, *gln
 * Returns the incomplete gamma function
 * Taken from Numerical Recipes in C, second edition
 * */
static double gammp(double a, double x)
{
	double gamser,gammcf,gln;
	
	if( x < 0.0 || a <= 0.0 )
		fprintf(stderr, "ERROR: (Sersic:gammp) Invalid arguments\n");
	if( x <  a + 1.0 )
	{
		gser(&gamser,a,x,&gln);
		return gamser;
	} 
	else 
	{
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}

/* --------------------------------------------------------------*/
/*  Pointers (Global variables?) used :
 * *gamser, *gln
 *
 * Taken from Numerical Recipes in C, second edition
 * */
static void gser(double *gamser, double a, double x, double *gln)
{
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if( x <= 0.0 ) 
	{
		if( x < 0.0 ) 
			fprintf(stderr,"ERROR: (Sersic:gser) x less than 0\n");
			
		*gamser=0.0;
		return;
	} 
	else 
	{
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++)
		{
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) 
			{
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		fprintf(stderr,"ERROR: (Einasto:gser) a too large, ITMAX too small\n");
	}
}

/* --------------------------------------------------------------*/
/* Pointers used (global variables?) :
 * *gammcf, *gln
 *
 * Taken from Numerical Recipes in C, second edition
 * */
static void gcf(double *gammcf, double a, double x, double *gln)
{
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) 
	{
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) 
		fprintf(stderr,"ERROR: (Einasto:gcf) a too large, ITMAX too small\n");
		
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}


#undef ITMAX
#undef EPS
#undef FPMIN
