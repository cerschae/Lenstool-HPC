#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        o_chi_flux                  */
/*      auteur:     Ghislain Golse          */
/*      date:       10/99               */
/*      place:      Toulouse            */
/****************************************************************
 * Calcul de chi2 deduit de la difference de flux entre celui
 * donne par la source "moyenne" et celui de l'image
 *
 * mag = -2.5 log f
 * dmag = -2.5 df / f  -->  df = abs( dmag * f / 2.5 )
 *
 * Global variables used :
 * - I
 */
void o_chi_flux(struct galaxie *arclet, double fluxS, double *da, double *sig2flux)
{
    double Dflux;
    const extern struct g_image I;

    // If Dmag is very small compared to flux then 
    // the approximation is dfv = ln10/2.5*dM*fv
    Dflux = I.Dmag * arclet->flux / 2.5;

    *sig2flux = Dflux * Dflux;
    *da = arclet->flux - fluxS / arclet->A;

}


