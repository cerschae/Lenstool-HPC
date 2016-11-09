#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        e_mass              */
/*      auteur:     Eric Jullo          */
/*      date:       12/07/06            */
/*      place:      ESO Chile           */
/****************************************************************
 * Return the mass of 1 clump inside a given radius in Msol.
 *
 * Sum the surface mass given by e_grad2() over an ellitical area.
 * In case of elliptical clump, <radius> is the half major axis.
 *
 * Parameters:
 * - icl : number of the clump
 * - radius : max radius in arcsec. (-1 for total mass)
 */

double e_mass(long int icl, double radius)
{
    const extern struct   pot     lens[];
    const extern struct   g_cosmo C;
    const extern struct   g_source  S;

    struct point    pi;     // sampling point
    struct ellipse ampli;   // amplification matrix at that point
    struct matrix   grad2;  // second derivate of the projected lens potential

    double mass;            // total mass
    double conv;            //
    double x, y;            // sampling point coordinates
    double xmin, xmax;      // sampling area
    double ymin, ymax;
    double dos;             // distcosmo1 to S.zs
    double dl0s;            // distcosmo2 between lens[0] and zs
    double dl;              // luminosity distance to the lens
    double a, b, c;         // amplification matrix diagonalisation coefs

    dos = distcosmo1(S.zs);
    dl0s = lens[icl].dlsds * dos;
    dl = distcosmo1(lens[icl].z);

    // all the distances are in arcsec
    xmax = radius / sqrt(1 + lens[icl].emass) + lens[icl].C.x;
    ymax = radius / sqrt(1 - lens[icl].emass) + lens[icl].C.y;
    xmin = -xmax;
    ymin = -ymax;

    mass = 0.;
    for ( x = xmin; x < xmax; x += 0.05 )
        for ( y = ymin; y < ymax; y += 0.05 )
        {
            pi.x = x;
            pi.y = y;
            grad2 = e_grad2_pot(&pi, icl);
            grad2.a *= dl0s;
            grad2.b *= dl0s;
            grad2.c *= dl0s;

            a = 1. - grad2.a;
            b = 1. - grad2.c;
            c = - grad2.b;

            ampli = formeli(a, b, c);

            mass += 1. - (ampli.a + ampli.b) / 2.;    // the convergence k
        }

    conv = MCRIT12 / C.h * 0.05 * 0.05 / dl0s / dos * dl;
    mass *= conv;

    return mass;
}

