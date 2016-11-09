#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        unmag               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse
 *****************************************************************
 * Return an ellipse that contains the eigenvalues of the amplification
 * matrix [ a=(1-k+gamma), b=(1-k-gamma) ] and the orientation of the proper
 * magnification axis.
 *
 * Used to convert image to source plane (using the lens equation)
 *
 * You have to multiply those eigenvalues by c2/2/DOL to get the true
 * eigenvalues.
 *
 * Here :
 * k = 0.5 * DLS/DS * ( d2Phixx + d2Phiyy )
 * gamma = DLS/DS * SQRT ( 0.25 * (d2Phiyy - d2Phixx)^2 + (d2Phixy)^2 )
 *
 * Parameters :
 * - position is the point where the amplification matrix is computed
 * - dl0s is distcosmo2 between lens[0] and zs
 * - dos is the normalised distance to the source
 * - zs is the redshift of the source
 *
 * Global variables used :
 * - in e_grad2() : G, lens, lens_table
 */
struct  ellipse e_unmag(const struct point *position, double dl0s, double dos, double zs)
{
    struct  ellipse ampli;
    double  A, B, C;
    struct  matrix  M;

    M = e_grad2(position, dl0s, zs);
    M.a /= dos;
    M.b /= dos;
    M.c /= dos;
    A = 1. - M.a;   // 1 - DLS/DS * d2phixx
    B = -M.b;   // - DLS/DS * d2phixy
    C = 1. - M.c;   // 1 - DLS/DS * d2phiyy
    ampli = formeli(A, B, C);


    return(ampli);
}

struct  ellipse e_unmag_gal(struct galaxie *image)
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
    ampli = formeli(A, B, C);


    return(ampli);
}
