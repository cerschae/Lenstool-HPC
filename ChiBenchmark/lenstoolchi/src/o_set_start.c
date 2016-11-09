#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_set_start         */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/
/* Set the starting values of each parameter for the lens i according
 * to the limits settings in the parameters file */

void  o_set_start(int i, int block[NLMAX][NPAMAX])
{
    extern struct ipot  ip;
    double a, b, lmin, lmax, x;

    int    ipx;

    for ( ipx = 0; ipx < ip.pmax ; ipx++ )
    {
        if ( block[i][ipx] > 0 )
        {
            lmin = o_get_lmin(i, ipx);
            lmax = o_get_lmax(i, ipx);
            a = ( lmin * 2. + lmax ) / 3.;
            b = ( lmin + 2. * lmax ) / 3.;
            x = o_get_lens(i, ipx);
            x = Min(b, x);
            x = Max(a, x);
            o_set_lens( i, ipx, x );
        }
    }
}
