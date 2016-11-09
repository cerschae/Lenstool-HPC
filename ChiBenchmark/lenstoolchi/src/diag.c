#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        diag                */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Same as formeli() but return the inverse of a and b, the proper
 * magnification axis.
 *
 * Used in e_test()/e_mag() to define the arclet size in the image plane
 * from the source size in the source plane.
 *
 * Parameters :
 * - a : 1 - dxx(Phi)
 * - b : - dxy(Phi)
 * - c : 1 - dyy(Phi)
 */
struct  ellipse diag(double a, double b, double c)
{
    struct  ellipse diago;

    diago = formeli(a, b, c);
    diago.a = 1 / diago.a;
    diago.b = 1 / diago.b;

    return (diago);
}
