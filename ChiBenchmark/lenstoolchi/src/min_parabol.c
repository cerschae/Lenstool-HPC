#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        min_parabol         */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/*                              */
/* recherche d'un minimum dans le cas ou on a 3 points      */
/* tels que (x1<x0<x2 ou x1>x0>x2) et (y0<y1<y2 ou y1>y2>y0)    */
/****************************************************************/

double  min_parabol(double x0, double y0, double x1, double y1, double x2, double y2)
{
    const extern  struct  g_mode          M;
    double  D1, D2, a, xmin;


    D1 = slope(x0, y0, x1, y1);
    D2 = slope(x0, y0, x2, y2);

    if (D1*D2 > 0.)
    {
        NPRINTF(stderr, "ERROR: x0 not an estimated minimum");
        NPRINTF(stderr, "\t%lf %lf %lf %lf %lf %lf %.3lf %.3lf\n",
                x0, y0, x1, y1, x2, y2, D1, D2);
        return(x0);
    };
    if (D1 == 0.)
        xmin = .5 * (x0 + x1);
    else if (D2 == 0.)
        xmin = .5 * (x0 + x2);
    else
    {
        a = slope(x1, D1, x2, D2);
        xmin = .5 * (x0 + .5 * (x1 + x2 - (D1 + D2) / a));
    };
    return(xmin);
}
