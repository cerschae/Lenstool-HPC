#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        e_mag               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Return an ellipse that contains the inverse of the eigenvalues of the
 * amplification matrix [ a=1/(1-k+gamma), b=1/(1-k-gamma) ] and the
 * orientation of the proper magnification axis.
 *
 * Used to convert source to image plane. (inverse of the lens equation)
 * Used in e_test.c
 */

struct  ellipse e_mag(struct point *position, double dl0s, double dos, double zs)
{
    struct  ellipse ampli;
    double  A, B, C;
    struct  matrix  M;

    M = e_grad2(position, dl0s, zs);
    M.a /= dos;
    M.b /= dos;
    M.c /= dos;
    A = 1. - M.a;
    B = -M.b;
    C = 1. - M.c;
    ampli = diag(A, B, C);

    return(ampli);
}

struct  ellipse e_mag_gal(struct galaxie *image)
{
    struct  ellipse ampli;
    double  A, B, C;
    struct  matrix  M;

    M = e_grad2_gal(image, NULL);
    M.a /= image->dos;
    M.b /= image->dos;
    M.c /= image->dos;
    A = 1. - M.a;
    B = -M.b;
    C = 1. - M.c;
    ampli = diag(A, B, C);

    return(ampli);
}

