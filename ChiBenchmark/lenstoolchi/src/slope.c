#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        slope               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/*                              */
/* determination de la pente entre 2 points         */
/****************************************************************/

double slope(double x0, double y0, double x1, double y1)
{
    const extern struct   g_mode M;
    if (x0 != x1)
        return((y1 - y0) / (x1 - x0));
    else
    {
        NPRINTF(stderr, "WARNING: slope not determined [%lf,%lf], 0. assumed\n", x0, x1);
        return(x1 - x0);
    }
}
