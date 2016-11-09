#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

static double   hern_f(double x);
static double   hern_gamma(double r, double rs, double kappas);
static double hern_kappa(double r, double rs, double kappas);

/****************************************************************/
/*              nom:            e_hernquist                          */
/*              auteur:         Eric Jullo, Carlo Giocoli (MOKA) */
/*              date:             06/13                        */
/*              place:          Marseille                        */
/****************************************************************
 * Global variables used :
 * - none
 * */
double   hern_dpl(double r, double rs, double kappas)
{
    double x = r / rs;
    double dpl = kappas * rs * x * (1. - hern_f(x)) / (x * x - 1.);
    return( dpl );
}

/* --------------------------------------------------------------*/
/* Global variables used :
 * - none
 * */
static double hern_kappa(double r, double rs, double kappas)
{

  double x = r / rs;
  double kappa = kappas * ((2. + x * x) * hern_f(x) - 3.) / (x * x - 1.) / (x * x - 1.);
  return(kappa);
}

/* --------------------------------------------------------------*
 * Return kappa for elliptical NFW (cf. Golse & Kneib 2002, eq 17)
 * Global variables used :
 * - none
 * */
double   hern_kappa_eps(double r, double rs, double theta, double kappas, double eps)
{
    double kappa_eps;

    kappa_eps = hern_kappa(r, rs, kappas) + eps * cos(2 * theta) * hern_gamma(r, rs, kappas);

    return(kappa_eps);
}

/* --------------------------------------------------------------*
 * Return gamma1 and gamma2 for elliptical Hernquist. 
 * (cf. Drumet-Montoyat et al. 2012 Eq A10-A12, Golse & Kneib 2002, eq 18)
 * The formalism is valid independently of the profile, namely NFW or Hernquist
 * Global variables used :
 * - none
 */
struct point   hern_gamma_eps(double r, double rs, double theta, double kappas, double eps)
{
    struct point gamma_eps;
    double kappa, gamma;

    kappa = hern_kappa(r, rs, kappas);
    gamma = hern_gamma(r, rs, kappas);

//    gamma_eps = sqrt( pow(gamma, 2.) + 2 * eps * cos(2 * theta) * kappa * gamma + eps * eps * ( pow(kappa, 2.) - pow(sin(2 * theta) * gamma, 2.) ) );  //NOT USED ANYMORE

    gamma_eps.x = gamma * cos(2. * theta) + eps * kappa;  // gamma1 
    gamma_eps.y = -sqrt(1. - eps * eps) * gamma * sin(2. * theta); // gamma2

    return(gamma_eps);
}

/* --------------------------------------------------------------*/
/* Global variables used :
 * - none
 * */
static double   hern_gamma(double r, double rs, double kappas)
{
    return( hern_dpl(r, rs, kappas) / r - hern_kappa(r, rs, kappas) );

}

/* --------------------------------------------------------------*/
/* Global variables used :
 * Taken from Keeton 2002, Eq 48
 * - none
 * */
static double   hern_f(double x)
{
    double f = 0;

    if ( x > 1 )
        f = atan(sqrt((x - 1.) * (x + 1.))) / sqrt(x * x - 1);  
    else if ( x < 1 )
        f = atanh(sqrt((1. - x) * (1. + x))) / sqrt(1. - x * x);
    else
        f = 1.;

    return(f);
}
