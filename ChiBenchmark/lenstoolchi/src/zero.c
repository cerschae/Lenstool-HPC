#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

static double   zerodich(double c1, double c2, double (*f)(double));

/****************************************************************/
/*      nom:        zero                */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse
 * Stop when *f(c1) and *f(c2) have the same signs.
 *
 * If f is the fz_fdlsds() function, return the solution z of
 * the equation f(z)=z_dlsds.
 *
 * - *f(c1) is always smaller than *f(c2)
 ****************************************************************/

double  zero(double c1, double c2, double (*f)(double))
{
    double  c;

    /*f(c1) is always smaller than *f(c2)*/
    if ((*f)(c1) > (*f)(c2))
    {
        c = c1;
        c1 = c2;
        c2 = c;
    }

    if ((*f)(c1)*(*f)(c2) >= 0.)
    {
        if ((*f)(c1) == 0.)
            return(c1);
        if ((*f)(c2) == 0.)
            return(c2);
        else
            return(c2);
    }
    else
        return(zerodich(c1, c2, f));
}

/****************************************************************/
/*If f is the fz_dlsds() function return */
static double   zerodich(double c1, double c2, double (*f)(double))
{
    double  c, z;

    c = (c1 + c2) / 2.;
    z = (*f)(c);
    if ((fabs(z) < PREC_ZERO) || (fabs(c1 - c2) < PREC_ZERO))
        return(c);
    else if (z*(*f)(c1) > 0)
        return(zerodich(c, c2, f));
    else
        return(zerodich(c1, c, f));
}

