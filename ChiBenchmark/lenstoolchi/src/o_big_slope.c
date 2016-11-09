#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_big_slope         */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            *
 ***************************************************************
 * Return the number of parameters which variation range is outside
 * the error bars.
 * Variation range is defined by excd and excu. See o_set_exc()
 */
int  o_big_slope()
{
    extern  struct  galaxie multi[][NIMAX];
    extern  struct  z_lim   zlim[];

    const extern  struct  g_grille    G;
    const extern  struct  g_image     I;
    const extern  struct  ipot    ip;
    const extern  struct  pot lens[];
    const extern  int block[][NPAMAX];
//  const extern  double  x1min,x1max,y1min,y1max;
//  const extern  double  x2min,x2max,y2min,y2max;
//  const extern  int izmin,izmax;

    register int    i, j;
    int stop; /*number of parameters that are not in the error bars */
    double  x0, x1, x2, y1, y2;

    stop = 0;

    if (lens[0].type != 10)
    {
        for (i = 0; i < G.no_lens; i++)
            for (j = 0; j < ip.pmax; j++)
                if (block[i][j] > 0)
                    stop += o_swi_big(i, j);
        /*o_swi_big() return 1 if the parameter has to be optimized
          according to its optimisation limits*/
    }
    /*else  spline mapping
        stop+=o_slope_sp();*/


    /* For each z_m_limit in image keyword put stop=1 if
     * dr (the Efficiency radio DLS/DOS) is outside the range
     * 0.3(ddmax-ddmin) by default.
     * The min and max values are saved in the o_keepz_min() function.
     * See o_set_ext() for excu, excd definition
     * */
    for (i = 0; i < I.nzlim; i++)
    {
        // check if this parameter optimisation is not blocked
        if (zlim[i].bk > 0)
        {
            x0 = multi[i][0].dr;
            x1 = multi[i][0].dr = x0 + zlim[i].excu * (zlim[i].ddmax - x0);
            y1 = o_chi();
            x2 = multi[i][0].dr = x0 - zlim[i].excd * (x0 - zlim[i].ddmin);
            y2 = o_chi();
            multi[i][0].dr = x0;

            if (fabs((x2 - x1)) <= zlim[i].dderr)
                zlim[i].bk = 0; // block the optimisation
            else
            {
                o_keepz_min(x1, x2, y1, y2, i); /*check if the new z has modified
                                            the global chi2 min/max values*/
                stop++;
            };
        };
    };

    return(stop);
}
