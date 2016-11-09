#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        zero (methode des tangentes)    */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/11/93            */
/*      place:      Toulouse            */
/****************************************************************
 * Global variables used :
 * - none
 */
double  zero_t(double c1, double c2, double (*f)(double))
{
    double  dc, c, fc, c0, fc0, tang, errc = .005;
    double  fc1, fc2;
    int end = 0, nb = 0;

    fc1 = (*f)(c1);
    fc2 = (*f)(c2);
    if (fabs(fc1) < fabs(fc2))
    {
        c = c2;
        c2 = c1;
        c1 = c;
        fc = fc2;
        fc2 = fc1;
        fc1 = fc;
    };

    do
    {
        nb++;
        tang = (c2 - c1) / (fc2 - fc1);
        c0 = c2 - fc2 * tang;

        if (c0 > 1.)
        {
            end = 1;
            c1 = c2 = 0.8;
        }
        else if (c0 < 0.)
        {
            end = 1;
            c1 = c2 = 0.;
        }
        else
        {
            fc0 = (*f)(c0);

            if ((fc0*fc2) < 0)
            {
                c1 = c0;
                fc1 = fc0;
            }
            else if (fabs(fc0) < fabs(fc2))
            {
                c1 = c2;
                fc1 = fc2;
                c2 = c0;
                fc2 = fc0;
            }
            else
                end = 1;
        };

        dc = fabs(c1 - c2);

    }
    while ((end == 0) && (dc > errc) && (nb < 50));

    return ((c1 + c2) / 2.);
}
