#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        o_mag_m             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * calcul de magnification pour les IMAGES MULTIPLES
 * Compute the deformation parameters thetaPot(thp), distortionPot(dp)
 * and tauPot(tp) for every arclet of the arclet list.
 *
 * arclet thetaPot (thp) = 2Gradxy(Phi)/(Gradyy(Phi) - Gradxx(Phi))
 * arclet distortion (dp) = (K^2 + G^2)/(K^2 - G^2) = trace(a^-2)/2det(a^-1)
 * arclet -tau (tp) = 2KG/(K^2 - G^2)
 *
 * Global variables used :
 * - in e_grad2() : G, lens, lens_table
 */
void    o_mag_m(int n, struct galaxie *arclet)
{
    register int    i;
    double  A, B, C;
    struct  matrix  MA;


    for (i = 0; i < n; i++)
    {
        MA = e_grad2_gal(&arclet[i], NULL);
        MA.a *= arclet[0].dr;
        MA.b *= arclet[0].dr;
        MA.c *= arclet[0].dr;
        MA.d *= arclet[0].dr;
        A = 1. - MA.a;  //Gradxx(Phi)
        B = -MA.b;  //Gradxy(Phi)
        C = 1. - MA.c;  //Gradyy(Phi)
        arclet[i].A = fabs(A * C - B * B);
        arclet[i].thp = PI / 2. + 0.5 * atan2(2.*B, (A - C));
        arclet[i].dp = (A * A + C * C + 2.*B * B) / arclet[i].A / 2.;
        arclet[i].tp = sqrt(arclet[i].dp * arclet[i].dp - 1.);
    }

    arclet[n] = arclet[0];
}
