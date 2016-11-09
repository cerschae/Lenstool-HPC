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
/****************************************************************/


void    o_mag(int n, struct galaxie *arclet)
{
    register int    i;
    double  A, B, C;
    struct  matrix  MA;


    for (i = 0; i < n; i++)
    {
        MA = e_grad2_gal(&arclet[i], NULL);
        MA.a *= arclet[i].dr;
        MA.b *= arclet[i].dr;
        MA.c *= arclet[i].dr;
        MA.d *= arclet[i].dr;
        A = 1. - MA.a;
        B = -MA.b;
        C = 1. - MA.c;
        arclet[i].A = fabs(A * C - B * B);
        arclet[i].thp = PI / 2. + 0.5 * atan2(2.*B, (A - C));
        arclet[i].dp = (A * A + C * C + 2.*B * B) / arclet[i].A / 2.;
        arclet[i].tp = sqrt(arclet[i].dp * arclet[i].dp - 1.);
    }

}
