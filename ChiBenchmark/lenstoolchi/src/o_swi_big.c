#include<stdio.h>
#include<math.h>
#include<float.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_swi_big           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            *
 ***************************************************************
 * Return 1 if the new parameter range is larger than the error bar
 * otherwise return 0.
 *
 * Parameters :
 * - i : index of the potential to optimize
 * - ipx : index of the parameter to optimize
 */
int  o_swi_big(int i, int ipx)
{
    extern double excu[][NPAMAX];
    extern double excd[][NPAMAX];
    extern int    block[][NPAMAX];

    double x0, x1, x2, y1, y2, lmin, lmax;

    x0 = x1 = x2 = y1 = y2 = DBL_MAX;

    lmin = o_get_lmin(i, ipx);
    lmax = o_get_lmax(i, ipx);
    x0 = o_get_lens(i, ipx);
    x1 = x0 + excu[i][ipx] * ( lmax - x0 );
    o_set_lens(i, ipx, x1);
    y1 = o_chi();
    x2 = x0 - excd[i][ipx] * ( x0 - lmin );
    o_set_lens(i, ipx, x2);
    y2 = o_chi();
    o_set_lens(i, ipx, x0);

    //TODO check if x0,x1,x2,y1,y2 have been modified

    /* if excd=excu (default value in o_set_ext() : 0.3) then
     *      x2-x1<err is equiv to 0.3(pmax-pmin)<err */
    if ( fabs(x2 - x1) <= o_get_err(i, ipx) )
    {
        block[i][ipx] = 0;
        return(0);
    }
    else
    {
        o_keep_min(x1, x2, y1, y2, i, ipx);
        return(1);
    }
}
