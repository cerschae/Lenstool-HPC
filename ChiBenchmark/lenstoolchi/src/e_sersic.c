#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

static double b(double n);
static double sersic_m(double r, double re, double n, double kappae);

// From Numerical recipes
static double gammln(double xx);
static double gammp(double a, double x);
static void gser(double *gamser, double a, double x, double *gln);
static void gcf(double *gammcf, double a, double x, double *gln);

/****************************************************************/
/*              name:           e_sersic                        */
/*              author:         Ardis Eliasdottir               */
/*              date:              01/2007                      */
/*              place:          Copenhagen                      */
/****************************************************************/

/* Global variables used :
 * - none
 * Given by re*alpha(x) where alpha(x) is the dimensionless deflection
 * angle (see eg. eq.8.3 in Schneider, Ehlers and Falco 1992)
 * */
double   sersic_dpl(double r, double re, double n, double kappae)
{
    double sersic_dpl;

    sersic_dpl = 0.;
    if ( r > 0 ) sersic_dpl = re * sersic_m(r, re, n, kappae) / (r / re);
    return( sersic_dpl );
}


/* --------------------------------------------------------------*/
/* Global variables used :
 * - none
 * */
double   sersic_kappa(double r, double re, double n, double kappae)
{
    return(kappae*exp(-b(n)*(pow(r / re, 1. / n) - 1.)));
}

/* --------------------------------------------------------------*/
/* Global variables used :
 * - none
 * Gamma(n)=exp(gammln(n))
 * Gamma(n,x)=Gamma(n)*gammq(n,x)
 * gamma(n,x)=Gamma(n)-Gamma(n,x)=Gamma(n)*gammp(n,x)
 * */
double   sersic_kappa_av(double r, double re, double n, double kappae)
{
    double sersic_kappa_av;
    double B = b(n);

    sersic_kappa_av = kappae * exp(B);
    if ( r > 0 ) sersic_kappa_av = kappae * 2.*pow(B, -2.*n) * n *
                                       exp(B) * exp(gammln(2.*n)) * gammp(2.*n, B * pow(r / re, 1. / n)) / r / r * re * re;

    return( sersic_kappa_av );
}

/* --------------------------------------------------------------*/
/* Global variables used :
 * - none
 * */
double   sersic_gamma(double r, double re, double n, double kappae)
{
    return(sersic_kappa_av(r, re, n, kappae) - sersic_kappa(r, re, n, kappae));
}

/* --------------------------------------------------------------*/
/* Global variables used :
 * - none
 * */
double   sersic_kappa_eps(double r, double re, double n, double theta, double kappae, double eps)
{
    double kappa_eps;

    kappa_eps = sersic_kappa(r, re, n, kappae) + eps * cos(2.*theta) * sersic_gamma(r, re, n, kappae);

    return(kappa_eps);
}

/* --------------------------------------------------------------*
 * Return gamma1 and gamma2 for a pseudo-elliptical potential
 * Global variables used :
 * - none
 */
struct point   sersic_gamma_eps(double r, double re, double n, double theta, double kappae, double eps)
{
    struct point gamma_eps;
    double kappa, gamma;

    kappa = sersic_kappa(r, re, n, kappae);
    gamma = sersic_gamma(r, re, n, kappae);

    //gamma_eps = sqrt(pow(sersic_gamma(r, re, n, kappae), 2.) + 2.*eps * cos(2.*theta) * sersic_kappa(r, re, n, kappae) * sersic_gamma(r, re, n, kappae) + eps * eps * (pow(sersic_kappa(r, re, n, kappae), 2.) - pow(sin(2.*theta) * sersic_gamma(r, re, n, kappae), 2.)));
    //

    gamma_eps.x = gamma * cos(2. * theta) + eps * kappa;  // gamma1 
    gamma_eps.y = -sqrt(1. - eps * eps) * gamma * sin(2. * theta); // gamma2

    return(gamma_eps);
}


/* --------------------------------------------------------------*/
/* Global variables used :
 * - none
 * This is an approximation (see Ciotti & Bertin 1999) to the
 * equation Gamma(2n,b(n))=Gamma(2n)/2
 * */
static double   b(double n)
{
    return( 2.*n - 1. / 3. + 4. / 405. / n + 46. / 25515. / n / n );
}

/* --------------------------------------------------------------*/
/* Global variables used :
 * - none
 * This calculates the (dimensionless) 2D mass within radius r
 * divided by pi (?), given by 2 Integrate(x' Kappa(x')_0^x)
 * */
static double sersic_m(double r, double re, double n, double kappae)
{
    double B = b(n);
    return 2.*kappae*pow(B, -2.*n)*exp(B)*n*exp(gammln(2.*n))*gammp(2.*n, B*pow(r / re, 1. / n));
}


/* --------------------------------------------------------------*/
/* Global variables used :
 * - none
 * This is the ln of the Gamma function
 * Taken from Numerical Recipes in C, second edition
 * */
static double gammln(double xx)
{
    double x, y, tmp, ser;
    double cof[6] = { 76.18009172947146,
                      -86.50532032941677,
                      24.01409824083091,
                      -1.231739572450155,
                      0.1208650973866179e-2,
                      -0.5395239384953e-5
                    };
    int j;

    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005*ser / x);
}


/* --------------------------------------------------------------*/
/* Pointers (Global variables?) used :
 * - *gammcf, *gamser, *gln
 * Returns the incomplete gamma function
 * Taken from Numerical Recipes in C, second edition
 * */
static double gammp(double a, double x)
{
    double gamser, gammcf, gln;

    if ( x < 0.0 || a <= 0.0 )
        fprintf(stderr, "ERROR: (Sersic:gammp) Invalid arguments\n");

    if ( x <  a + 1.0 )
    {
        gser(&gamser, a, x, &gln);
        return gamser;
    }
    else
    {
        gcf(&gammcf, a, x, &gln);
        return 1.0 - gammcf;
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
    double sum, del, ap;

    *gln = gammln(a);
    if ( x <= 0.0 )
    {
        if ( x < 0.0 )
            fprintf(stderr, "ERROR: (Sersic:gser) x less than 0\n");

        *gamser = 0.0;
        return;
    }
    else
    {
        ap = a;
        del = sum = 1.0 / a;
        for (n = 1; n <= ITMAX; n++)
        {
            ++ap;
            del *= x / ap;
            sum += del;
            if (fabs(del) < fabs(sum)*EPS)
            {
                *gamser = sum * exp(-x + a * log(x) - (*gln));
                return;
            }
        }
        fprintf(stderr, "ERROR: (Sersic:gser) a too large, ITMAX too small\n");
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
    double an, b, c, d, del, h;

    *gln = gammln(a);
    b = x + 1.0 - a;
    c = 1.0 / FPMIN;
    d = 1.0 / b;
    h = d;
    for (i = 1; i <= ITMAX; i++)
    {
        an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = b + an / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del - 1.0) < EPS) break;
    }
    if (i > ITMAX)
        fprintf(stderr, "ERROR: (Sersic:gcf) a too large, ITMAX too small\n");

    *gammcf = exp(-x + a * log(x) - (*gln)) * h;
}


#undef ITMAX
#undef EPS
#undef FPMIN
