#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        magfast             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Global variables used :
 * - in e_grad2() : G, lens, lens_table
 */
double  diff_mag(struct galaxie *arclet, struct point *guess)
{
    double  ampguess, A, B, C;
    struct  matrix  MA, MG;

    MA = e_grad2_gal(arclet, NULL);
    MA.a /= arclet->dos;
    MA.b /= arclet->dos;
    MA.c /= arclet->dos;

    A = 1. - MA.a;
    B = -MA.b;
    C = 1. - MA.c;
    arclet->A = fabs(A * C - B * B);

    MG = e_grad2(guess, arclet->dl0s, arclet->z);
    MG.a /= arclet->dos;
    MG.b /= arclet->dos;
    MG.c /= arclet->dos;

    A = 1. - MG.a;
    B = -MG.b;
    C = 1. - MG.c;
    ampguess = fabs(A * C - B * B);


    return( arclet->A - ampguess );
}
