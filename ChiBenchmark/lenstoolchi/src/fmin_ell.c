#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        fmin_ell            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/11/93            */
/*      place:      ESO LaSilla         */
/****************************************************************
 *
 *
 * Global variables used :
 * - elix, SC
 * - in e_unmag() : G, lens, lens_table
 */
double  fmin_ell(double dl0s, double dos, double zs)
{
    const extern  double  elix;
    const extern  struct  point   SC;
    double  dp, tp, u, q;
    struct  ellipse ampli;

    ampli = e_unmag(&SC, dl0s, dos, zs);
    q = ampli.a / ampli.b;
    u = 1. / q;
    tp = q - u;
    dp = q + u;

    return(dp*elix - tp);
}


