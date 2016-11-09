#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        formeli             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/*****************************************************************
 * Diagonalisation d'une matrice de forme ou de deformation
 * Matrice de forme : M= | a  b |
 *                       | b  c |
 * lambda et mu : racines de det(M - XI)=0
 *
 * Return the diagonalized form of the input matrice M. Theta is the
 * shear orientation, a and b are the proper magnification axis.
 *
 * If :
 * 1) a = 1 - [g2].a = 1 - DLS/DS * d2phixx
 * 2) c = 1 - [g2].c = 1 - DLS/DS * d2phiyy
 * 3) b = - [g2].b = - DLS/DS * d2phixy
 * then :
 *     with gamma = DLS/DS * SQRT ( 0.25 * (d2Phiyy - d2Phixx)^2 + (d2Phixy)^2 )
 *     with k = 0.5 * DLS/DS * ( d2Phixx + d2Phiyy ) (cf phd_JPK eq 2.55)
 *
 *  delta = (a-c)^2+4b^2
 *        = [ DLS/DS * (d2Phiyy - d2Phixx) ]^2 + 4*[ DLS/DS * d2Phixy ]^2
 *        = 4*gamma^2 (cf phd_JPK eq 2.56)
 *
 *  lambda = 0.5*(a + c + SQRT(delta) )
 *         = .5*( 2 - DLS/DS * (d2Phixx + d2Phiyy) + 2*gamma)
 *         = ( 1 - k + gamma )
 *
 *      mu = ( 1 - k - gamma )
 *
 * Global variables used :
 * - none
 */
struct  ellipse formeli(double a, double b, double c)
{
    struct  ellipse eli;
    double  e, delta, lambda, mu;

    // eq carateristique : det(M-xI) = 0
    delta = (a - c)*(a - c) + 4*b*b;    //  4*gamma^2 (cf phd_JPK eq 2.56)
    e = sqrt(delta);  /*e is 2 * shear, ie 2*gamma*/
    lambda = .5*(a + c + e);  // 1 - k + gamma
    mu = .5*(a + c - e);      // 1 - k - gamma

    eli.a = lambda;
    eli.b = mu;
    if (lambda != mu && fabs(b) > 1e-5)
        eli.theta = atan2(lambda - a, b); // cf phd_JPK eq 2.58, and
    // tan(theta)= ( -cos(2theta) +- 1 ) / sin(2theta)
// ADDED by EJ 29/11/2007
    else if ( a >= c ) // ellipse aligned along the major axis of magnification
        eli.theta = 0.;
    else
        eli.theta = PI / 2.;    // ellipse aligned along the minor axis of magnification

    return(eli);
}
