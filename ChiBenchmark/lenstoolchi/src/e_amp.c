#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        e_amp               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            *
 ****************************************************************
 * Return the determinant of the amplification matrix ie the amplification
 * at the position position.
 *
 * - position is the position where to compute the amplification
 * - dl0s is the distance between lens[0] and zs
 * - dos is the normalised distance distcosmo1(zs)
 * - zs is the redshift of the source
 *
 * Global variables read :
 * - in e_grad2() : G, lens, lens_table
 */
double  e_amp(const struct point *position, double dl0s, double dos, double zs)
{
    struct  matrix M;

    M = e_grad2(position, dl0s, zs);
    M.a /= dos;
    M.b /= dos;
    M.c /= dos;

    return(((1. - M.a)*(1. - M.c) - M.b*M.b));
}

double  e_amp_gal(struct galaxie *image, double *np_b0)
{
    struct  matrix M;

    M = e_grad2_gal(image, np_b0);
    M.a /= image->dos;
    M.b /= image->dos;
    M.c /= image->dos;

    return(((1. - M.a)*(1. - M.c) - M.b*M.b));
}
