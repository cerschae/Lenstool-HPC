#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        min_slope_par           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/*                              */
/* recherche d'un minimum dans le cas ou on a 3 points      */
/* tels que (x1<x0<x2 ou x1>x0>x2) et (y1<y0<y2 ou y1>y0>y2)    */
/****************************************************************/

double min_slope_par(double x0, double y0, double x1, double y1, double x2, double y2)
{
    double  D1, D2, a, xmin;

    D1 = slope(x0, y0, x1, y1);
    D2 = slope(x0, y0, x2, y2);
    a = slope(x1, D1, x2, D2);

    if (a < PREC_ZERO)
    {
        if (y1 < y2)
            return(x1);
        else
            return(x2);
    }
    else
    {
        xmin = .5 * (x0 + .5 * (x1 + x2 - (D1 + D2) / a));
        return(xmin);
    }
}
