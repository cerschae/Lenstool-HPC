/**
* @file   module_cosmodistances.cpp
* @Author Thomas Jalabert, EPFL (me@example.com)
* @date   July 2015
* @version 0,1
* @brief  Library for the computation of cosmological ratios
*
* compute the cosmological ratio of the distances between the lens and the source and the lens and the observer
*
*/




/// include header file
#include <string>
#include <stdio.h>
#include <math.h>
#include <cstring>
#include <stdlib.h>
#include "module_cosmodistances.h"





// Declare static functions that will only be used in this module
// The functions are defined further below
static double module_cosmodistances_cosmo_root(double z,cosmo_param C);
static double module_cosmodistances_sk(double x, double k);
static double module_cosmodistances_chi1(double z,cosmo_param C);
static double module_cosmodistances_chi2(double z1, double z2,cosmo_param C);
static double module_cosmodistances_chiz(double z,cosmo_param C);
static double module_cosmodistances_integral_chiz_ab(double a, double b,cosmo_param C);



// Function defintions
//==========================================================================================================

/** @brief Calculates ratio distance_(lens-source)/distance_(source)
*   Calculates ratio distance_(lens-source)/distance_(source)
* @param nsetofimages     number of set of images
* @param nImagesSet     set of images
* @param z_lens        redshift of lens
* @param source        sources
* @param cosmoratio     variable where result is stored
* @param cosmopar      cosmological parameter
*/
void module_cosmodistances_lensSourceToSource( const int nsetofimages, int nImagesSet[], double z_lens, galaxy source[], double cosmoratio[], cosmo_param cosmopar){


        //int imageCounter = 0;  // Count the total number of images up to now
  for(int i=0; i<nsetofimages; i++){
    //printf("z_lens %f , imag.redshift %f \n" , z_lens,source[0].redshift);
    cosmoratio[i]=module_cosmodistances_lensSourceToObserverSource(z_lens,source[i].redshift, cosmopar); /// lens efficiency (angular distance)
                //imageCounter += (nImagesSet[i] - 1);  // We add the number of images to skip to get to the next set
  }
}




/** @brief Return the angular distance DA(observator(z=0),object(z)) (no unit)
* Multiply observerObject(z) by c/H0 to get the true value in Mpc.
* observerObject(z) * c/H0 = DA(0,z) = 1/(1+z) * Sk( integral( 0,z,c*dz/H(z) ) )

* @param z          redshift of object
* @param cosmopar      cosmological parameter
*/
double module_cosmodistances_observerObject(double z, cosmo_param cosmopar)

{
    double g;

    if (cosmopar.omegaX == 0.)
    {
        g = module_cosmodistances_cosmo_root(z,cosmopar);
        // Reformulation of the Mattig relation of OL = OK = 0 (De Sitter)
        return(2.*((1. - cosmopar.omegaM - g)*(1. - g)) / cosmopar.omegaM / cosmopar.omegaM / (1. + z) / (1. + z));
    }
    else
    {
        if (cosmopar.curvature != 0.)
            return(module_cosmodistances_sk(module_cosmodistances_chi1(z,cosmopar)*sqrt(fabs(cosmopar.curvature)), cosmopar.curvature)
                   / (1 + z) / sqrt(fabs(cosmopar.curvature)));
        else
            return(module_cosmodistances_chi1(z,cosmopar) / (1 + z));
    }
}





/** @brief Return the angular distance DA(object(z1),object(z2)) divided by c/H0
* Return the angular distance DA(object(z1),object(z2)) divided by c/H0
* objectObject(z1,z2) * c/H0 = DA(z1,z2) = 1/(1+z2) * Sk( integral( z1,z2,c*dz/H(z) ) )
*
* @param z1          redshift of object 1
* @param z2          redshift of object 2
* @param cosmopar      cosmological parameter
*/
double module_cosmodistances_objectObject(double z1, double z2, cosmo_param cosmopar)
{
    double  g1, g2;

    if ( z1 >= z2 )
        return 0.;

    if (cosmopar.omegaX == 0.)
    {
        g1 = module_cosmodistances_cosmo_root(z1,cosmopar);
        g2 = module_cosmodistances_cosmo_root(z2,cosmopar);
        // Mattig relation for a De Sitter Universe
        return(2.*((1. - cosmopar.omegaM - g1*g2)*(g1 - g2))
               / cosmopar.omegaM / cosmopar.omegaM / (1. + z1) / (1. + z2) / (1. + z2));
    }
    else
    {
        if ( cosmopar.curvature != 0. )
            return(module_cosmodistances_sk(module_cosmodistances_chi2(z1, z2,cosmopar)*sqrt(fabs(cosmopar.curvature)), cosmopar.curvature)
                   / (1 + z2) / sqrt(fabs(cosmopar.curvature)));
        else
            return(module_cosmodistances_chi2(z1, z2,cosmopar) / (1 + z2));
    }
}




/****************************************************************/

 /** @brief Return the lens efficacity E=DA(LS) / DA(OS)
* If zl > zs, return 0.
*
* @param zl          redshift lens
* @param zs          redshift source
* @param cosmopar      cosmological parameter
*/

double module_cosmodistances_lensSourceToObserverSource(double zl, double zs, cosmo_param cosmopar)
{
    if ( zl >= zs ){
  printf("*******************\nWarning, a source is between the lens and the observer\n*******************\n");
  printf("Zl: %f , ZS: %f \n", zl, zs);
        return (0.);
    }
        return (module_cosmodistances_objectObject(zl, zs, cosmopar) / module_cosmodistances_observerObject(zs, cosmopar));
}



 /** @brief Calculate square root
*
* @param z      redshift
* @param C      cosmological parameter
*/

// Calculate square root
static double module_cosmodistances_cosmo_root(double z,cosmo_param C)
{
    return(sqrt(1. + C.omegaM*z));
}




// Calculate S_k for cosmology
 /** @brief Calculate S_k for cosmology
*
* @param x
* @param k      curvature
*/
static double module_cosmodistances_sk(double x, double k)
{
    if (k > 0)
        return(sin(x));
    else if (k < 0)
        return(sinh(x));
    else
        return(x);
}


 /** @brief Return 1 / H(z) multiplied by H0
* chiz(z) / H0 = 1 / H(z) = 1/H0/(1+z)/sqrt( sum(i, Omega_i*(1+z)^(3*(w_i + 1) ) )
*
* @param z      redshift
* @param C      cosmological parameter
*/

static double module_cosmodistances_chiz(double z,cosmo_param C)
{
    double x;
    double yy, yyy, y4;
    double r0, e1, e2, frac;

    x = 1 + z;


    switch (C.model)       //TV  CPL Model
    {
        case(1):
            x = -x * x * C.curvature + x * x * x * C.omegaM + C.omegaX * pow(x, 3 * (1 + C.wX + C.wa)) * exp(-3 * C.wa * z / x);
            break;
        case(2):        //TV Cardassian (wx is q, wa is n)
            yy = pow ( (1.+C.curvature)/C.omegaM  ,C.wX);
            yyy = (yy - 1.) * pow( x,3.*C.wX*(C.wa-1.) );
            y4 = pow(1.+yyy,1./C.wX);
            x = -x * x * C.curvature + x * x * x * C.omegaM*y4;
            break;
        case(3):        //TV Interacting DE Model (wa is delta)
            yy = C.omegaX*pow( x,3.*(1.+C.wX) );
            yyy = ( C.omegaM/(C.wa+3.*C.wX) ) * (  C.wa*pow(x,3.*(1.+C.wX))   +  3.*C.wX*pow(x,(3.-C.wa))   );
            x = -x * x * C.curvature + yy +yyy;
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
            printf("ERROR: Unknown cosmological model %d\n", C.model);
            exit(-1);

    }

    if ( x <= 0 )
    {
        printf("ERROR : H^2(z)<=0 produced (z,omegaM,omegaX,wX,wa) = (%.3lf,%.3lf,%.3lf,%.3lf,%.3lf)\n",
                 z, C.omegaM, C.omegaX, C.wX, C.wa );
         exit(-1);
    }
    return( 1. / sqrt(x) );
}



 /** @brief Return the proper distance Dp(0,z) divided by c/H0
* chi1(z) * c/H0 = Dp(0,z) = integral( 0, z, c*dz / H(z) )
*
* @param z      redshift
* @param C      cosmological parameter
*/
static double module_cosmodistances_chi1(double z,cosmo_param C)
{
    double rez;
    rez = module_cosmodistances_integral_chiz_ab(0., z,C);
    return rez;
}



 /** @brief Return the proper distance Dp(z1,z2) divided by c/H0
* chi2(z) * c/H0 = Dp(z1,z2) = integral( z1, z2, c*dz / H(z) )
*
* @param zl      redshift lens
* @param zs      redshift source
* @param C      cosmological parameter
*/
static double module_cosmodistances_chi2(double z1, double z2,cosmo_param C)
{
   double rez;
   rez = module_cosmodistances_integral_chiz_ab(z1, z2,C);
   return rez;
}


 /** @brief compute the integral of H0/H(z) by trapezes method
* compute the integral of H0/H(z) by trapezes method
* @param a,b, cosmologicalparameters cosmopar
*
* @param a
* @param b
* @param C    cosmological parameter
*/
static double module_cosmodistances_integral_chiz_ab(double a, double b,cosmo_param C)
{
   int i,nit;
   double res, epsilon=1.e-5; /// accuracy of the integration

   nit=(b-a)/epsilon;
   res=epsilon*(module_cosmodistances_chiz(a,C)+module_cosmodistances_chiz(b,C))/2.;

   for(i=1;i<nit;i++)
   {
       res=res+epsilon*module_cosmodistances_chiz(a+i*epsilon,C);
   }

   return res;
}






 /** @brief Debug function to print the calculated output of the functions
  * Debug function to print the calculated output of the functions
*/

// Debug function to print the calculated output of the functions

int module_cosmodistances_debug(int runmode[], double strongLensingRatios_lensSourceToSource[], double weakLensingRatios_lensSourceToSource[], double weakLensing_observerSource[], int numberCleanLens, double cleanlensRatios_lensSourceToSource[], double cleanlens_observerSource[], std::string DEBUG)
{
if (strcasecmp(DEBUG.c_str(), "True") == 0)  // If we are in debug mode
{

  if (runmode[0] == 1 or runmode[1] == 1)  // If we have strong lensing
  {
      printf("DEBUG: Strong lensing cosmo ratios D_LS/D_OS for first 2 sets: %lf, %lf\n\n", strongLensingRatios_lensSourceToSource[0], strongLensingRatios_lensSourceToSource[1]);
  };
  if (runmode[2] == 1)  // We have weak lensing
  {
      printf("DEBUG: Weak lensing D_LS/D_OS for first 2 arclets: %lf, %lf. Weak lensing D_OS for first 2 arclets: %lf, %lf.\n\n", weakLensingRatios_lensSourceToSource[0], weakLensingRatios_lensSourceToSource[1], weakLensing_observerSource[0], weakLensing_observerSource[1]);
  };
  if (numberCleanLens > 0)  // We have cleanlens mode
  {
      printf("DEBUG: Cleanlens D_LS/D_OS for first 2 sources: %lf, %lf. Cleanlens D_OS for first 2 sources: %lf, %lf.\n\n", cleanlensRatios_lensSourceToSource[0], cleanlensRatios_lensSourceToSource[1], cleanlens_observerSource[0], cleanlens_observerSource[1]);
  };


};

return 0;
}





