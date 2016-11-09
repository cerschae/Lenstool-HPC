#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<structure.h>
#include "lt.h"
#include "gsl/gsl_integration.h"
#include "errors.h"

/****************************************************************/
/*      nom:        cosmoratio          */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

static double gg(double z);
static double sk(double x, double k);
static double chi1(double z);
static double chi2(double z1, double z2);
static double chiz_gsl(double z, void* param);
static double integral_chiz_ab(double a, double b);

/**************************************************************/

double distcosmo1(double z)
/* Return the angular distance DA(observator(z=0),object(z)) (no unit)
 * Multiply discosmo1(z) by c/H0 to get the true value in Mpc.
 * distcosmo1(z) * c/H0 = DA(0,z) = 1/(1+z) * Sk( integral( 0,z,c*dz/H(z) ) )
 *
 * Global variables used :
 * - C
 */
{
    const extern struct g_cosmo C;
    double g;

    if (C.omegaX == 0.)
    {
        g = gg(z);
        // Reformulation of the Mattig relation of OL = OK = 0 (De Sitter)
        return(2.*((1. - C.omegaM - g)*(1. - g)) / C.omegaM / C.omegaM / (1. + z) / (1. + z));
    }
    else
    {
        if (C.kcourb != 0.)
            return(sk(chi1(z)*sqrt(fabs(C.kcourb)), C.kcourb)
                   / (1 + z) / sqrt(fabs(C.kcourb)));
        else
            return(chi1(z) / (1 + z));
    }
}

/**************************************************************/
/* Return the angular distance DA(object(z1),object(z2)) divided by c/H0
 * distcosmo2(z1,z2) * c/H0 = DA(z1,z2) = 1/(1+z2) * Sk( integral( z1,z2,c*dz/H(z) ) )
 *
 * Global variables used :
 * - C
 * - in gg() : C
 * - in chi2() : C
 */
double distcosmo2(double z1, double z2)
{
    const extern struct   g_cosmo C;
    double  g1, g2;

    if ( z1 >= z2 )
        return 0.;

    if (C.omegaX == 0.)
    {
        g1 = gg(z1);
        g2 = gg(z2);
        // Mattig relation for a De Sitter Universe
        return(2.*((1. - C.omegaM - g1*g2)*(g1 - g2))
               / C.omegaM / C.omegaM / (1. + z1) / (1. + z2) / (1. + z2));
    }
    else
    {
        if ( C.kcourb != 0. )
            return(sk(chi2(z1, z2)*sqrt(fabs(C.kcourb)), C.kcourb)
                   / (1 + z2) / sqrt(fabs(C.kcourb)));
        else
            return(chi2(z1, z2) / (1 + z2));
    }
}

/****************************************************************/
/* Return the lens efficacity E=DA(LS) / DA(OS)
 *
 * If zl > zs, return 0.
 *
 * Global variables used :
 * - in distcosmo1() : C
 * - in distcosmo2() : C
 */
double dratio(double zl, double zs)
{
    if ( zl >= zs )
        return (0.);
    else
        return (distcosmo2(zl, zs) / distcosmo1(zs));
}

/* Same as dratio() but for arclet structures
 */
void dratio_gal(struct galaxie *arclet, double zl)
{
    if( zl < arclet->z )
        arclet->dl0s = distcosmo2(zl, arclet->z);
    else
        arclet->dl0s = 0.;
    arclet->dos = distcosmo1(arclet->z);
    arclet->dr = arclet->dl0s / arclet->dos;
}

/**************************************************************
 * Return the luminosity distance DL(0,z) = (1+z)^2 * DA(0,z)
 */
double dlumcosmo1(double z)
{
    return(distcosmo1(z)*(1. + z)*(1. + z));
}

/**************************************************************/

double distprime1(double z)
{
    const extern  struct  g_cosmo C;
    double   g;

    if (C.omegaX == 0.)
    {
        g = 1. / (1. + z);
        return(g*(g / gg(z) - distcosmo1(z)*2.*g));
    }
    else
    {
        return((distcosmo1(z + 0.1) - distcosmo1(z)) / 0.1);
    }
}

/**************************************************************
 * Global variables used :
 * - in gg() : C
 * - in chi2() : C
 */
double distprime2(double z1, double z2)
{
    const extern struct   g_cosmo C;
    double  g2; //g1

    if (C.omegaX == 0.)
    {
        g2 = 1. / (1. + z2);
        return(g2*(g2 / gg(z2) - distcosmo2(z1, z2)*2.*g2));
    }
    else
    {
        return((distcosmo2(z1, z2 + 0.1) - distcosmo2(z1, z2)) / 0.1);
    }
}

/**************************************************************/

double dratioprime(double zl, double zs)
{
    double  ds;

    ds = distcosmo1(zs);
    return( (distprime2(zl, zs) - distcosmo2(zl, zs)*distprime1(zs) / ds) / ds );
}

/**************************************************************/
/* Global variables used :
 * - in chiz() : C
 * */
static double gg(double z)
{
    const extern struct g_cosmo C;

    return(sqrt(1. + C.omegaM*z));
}
/**************************************************************
 * Global variables used :
 * - none
 */

static double sk(double x, double k)
{
    if (k > 0)
        return(sin(x));
    else if (k < 0)
        return(sinh(x));
    else
        return(x);
}

/**************************************************************
 * Return 1 / H(z) multiplied by H0
 *   chiz(z) / H0 = 1 / H(z) = 1/H0/(1+z)/sqrt( sum(i, Omega_i*(1+z)^(3*(w_i + 1) ) )
 *   with w_M = 0 and w_k = -1/3
 *
 * Global variables used :
 * - C
 * */
double chiz(double z)
{
    const extern struct g_cosmo C;
    double x;
    double yy, yyy, y4;
    double r0, e1, e2, frac;

    x = 1 + z;
    
    
    switch (C.model)       //TV  CPL Model
    {
        case(1):
            x = -x * x * C.kcourb + x * x * x * C.omegaM + C.omegaX * pow(x, 3 * (1 + C.wX + C.wa)) * exp(-3 * C.wa * z / x);
            break;
        case(2):        //TV Cardassian (wx is q, wa is n)
            yy = pow ( (1.+C.kcourb)/C.omegaM  ,C.wX);
            yyy = (yy - 1.) * pow( x,3.*C.wX*(C.wa-1.) );
            y4 = pow(1.+yyy,1./C.wX);
            x = -x * x * C.kcourb + x * x * x * C.omegaM*y4;
            break; 
        case(3):        //TV Interacting DE Model (wa is delta)
            yy = C.omegaX*pow( x,3.*(1.+C.wX) );
            yyy = ( C.omegaM/(C.wa+3.*C.wX) ) * (  C.wa*pow(x,3.*(1.+C.wX))   +  3.*C.wX*pow(x,(3.-C.wa))   );
            x = -x * x * C.kcourb + yy +yyy;
            break;
        case(4):        //TV Holographic Ricci Scale with CPL
            r0 = C.omegaM/(1.-C.omegaM);
            e1 =  (3./2.)*(  (1.+r0+C.wX+4*C.wa)/(1.+r0+3.*C.wa) );
            e2 =  (-1./2.)*(  (1.+r0-3.*C.wX)/(1.+r0+3.*C.wa) )  ;
            frac =  ( 1.+r0+3.*C.wa*(x-1.)/x )/(1.+r0)  ;
            yy = pow(x,e1);
            yyy = pow(frac,e2);
            x = (yy*yyy)*(yy*yyy);
            break;
        default:
            fprintf(stderr, "ERROR: Unknown cosmological model %d\n", C.model);
            exit(E_COSMO_MODEL);
            break;
    }
    
    if ( x <= 0 )
    {
        fprintf( stderr, "ERROR : H^2(z)<=0 produced (z,omegaM,omegaX,wX,wa) = (%.3lf,%.3lf,%.3lf,%.3lf,%.3lf)\n",
                 z, C.omegaM, C.omegaX, C.wX, C.wa );
        exit(-1);
    }
    return( 1. / sqrt(x) );
}
/**************************************************************
 * Return the proper distance Dp(0,z) divided by c/H0
 *  chi1(z) * c/H0 = Dp(0,z) = integral( 0, z, c*dz / H(z) )
 *
 * Global variables used :
 * - in chiz() : C
 * */
static double chi1(double z)
{
    double rez;
    rez = integral_chiz_ab(0., z);
    return rez; 
}
/**************************************************************
 * Return the proper distance Dp(z1,z2) divided by c/H0
 *  chi2(z) * c/H0 = Dp(z1,z2) = integral( z1, z2, c*dz / H(z) )
 *
 * Global variables used :
 * - in chiz() : C
 */
static double chi2(double z1, double z2)
{
   double rez;
   rez = integral_chiz_ab(z1, z2);
   return rez;
}

static double chiz_gsl(double z, void* param)
{
   return chiz(z);
}

static double integral_chiz_ab(double a, double b)
{
   gsl_function f;
   f.function = &chiz_gsl;
   const int limit = 10;
   gsl_integration_workspace * w1 = gsl_integration_workspace_alloc(limit);
   double rez, err;
   gsl_integration_qag(&f, a, b, 0, 1e-6, limit, GSL_INTEG_GAUSS15, w1 , &rez ,&err);    
   gsl_integration_workspace_free(w1);   
   return rez;
}
