#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        sig2posSe         */
/*      auteur:     Eric Jullo          */
/*      date:       01/07               */
/*      place:      Marseille           */
/****************************************************************
 * Assuming that the deflection angle is linear between the observed
 * and predicted image, it is possible to analytically relate the
 * image plane chi2 to the source plane chi2.
 *
 * chi2_img = (xi - xiobs) = (xs - bar(xs)) / (1 - d2xxPhi)
 *
 * The observational uncertainties are the image ellipticities a & b.
 *
 * Global variables used :
 *
 */
void sig2posSe(struct galaxie *multi, double *sigx2, double *sigy2)
{

    struct matrix grad2;

    grad2 = e_grad2_gal(multi, NULL);
    grad2.a /= multi->dos;
    grad2.c /= multi->dos;
    *sigx2 = (1. - grad2.a) * (1. - grad2.a) * multi->E.a * multi->E.a;
    *sigy2 = (1. - grad2.c) * (1. - grad2.c) * multi->E.b * multi->E.b;
}
