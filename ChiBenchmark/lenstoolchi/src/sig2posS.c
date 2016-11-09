#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        sig2posS            */
/*      auteur:     Ghislain Golse          */
/*      date:       10/99               */
/*      place:      Toulouse            */
/****************************************************************/
/* Return the sigma in the source plan corresponding to sig2pos.
 *
 * The error image plane (1sigma) is converted in etendue (assuming
 * circular seeing AND FLAT brightness profile cf Kneib1993 Eq2.97, Eq2.91).
 * etendue_{image plane}  =  2. * sigma_pos^2
 * etendue_{src plane}    =  A^-1  etendue_img (cf Kneib1996 Eq24)
 * sigma_{src plane}^2    =  etendue_src / 2.
 *
 * Parameters :
 * - sig2pos : position error for the I.n_mult systems
 * - n : number of arclets for the current system
 * - n_famille:  source index  
 *
 * Global variables used :
 * - amplifi
 */
double   sig2posS(double sig2pos, int n, int n_famille, double *np_b0)
{
    const extern double   amplifi[NFMAX][NIMAX];
    extern struct galaxie   multi[NFMAX][NIMAX];

    double   sigma_src;  // error source plane
    double   ainv;  // averaged amplification^-1, ie averaged etendue in source plane
    // for a 1 arcsec^2 etendue in image plane.
    int i;

    ainv = 0.;

    for ( i = 0 ; i < n ; i++ )
        ainv += fabs(e_amp_gal(&multi[n_famille][i], np_b0)); //1. / fabs(amplifi[n_famille][i]);

    ainv = ainv / n;

    sigma_src = ainv * sig2pos;

    return(sigma_src);
}


