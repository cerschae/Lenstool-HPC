#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        o_min_slope         */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/
/* facteur de blocage:  1: amin amax bloques    */
/*          2: amin amax mous   */
/*          3: amin mou, amax bloque    */
/*          4: amin bloque, amax mou    */
/* Parameter :
 * - y0 : current value of Khi2
 ****************************************************************/


double  o_min_slope(double y0)
{
    extern  struct  g_mode          M;
    extern  struct  g_grille    G;
    extern  struct  galaxie multi[][NIMAX];
    extern  struct  pot lens[];
    extern  struct  z_lim   zlim[];
    extern  int block[][NPAMAX];
    extern  double  excu[][NPAMAX];
    extern  double  excd[][NPAMAX];
    extern  double  x1min, x2min, y1min, y2min;
    extern  int izmin, ipmin, ilsmin;
    extern  int imapmin, jmapmin;
    extern  double  **map_p;
    extern  double  **map_axx;
    extern  double  **map_ayy;
    extern  double  **tmp_p;

    double  x0, xmin, ymin, amin, amax, aerr;
    int i;

    /* If the minimum chi2 comes from a lens limit parameter */
    if (ipmin >= 0)
    {
        x0 = o_get_lens(ilsmin, ipmin);
        if ((x0 != x1min) && (x0 != x2min))
            xmin = min_slope_par(x0, y0, x1min, y1min, x2min, y2min);
        else
            xmin = x1min;

        o_set_lens(ilsmin, ipmin, xmin);
        ymin = o_chi();

        if (ymin > y1min)
            for (i = 0; i < 10; i++)
            {
                xmin = min_slope_par(x1min, y1min, x0, y0, xmin, ymin);
                o_set_lens(ilsmin, ipmin, xmin);
                ymin = o_chi();
                if (ymin < y1min)
                    i += 10;
            };

        if (ymin < y1min)
        {
            amin = o_get_lmin(ilsmin, ipmin);
            amax = o_get_lmax(ilsmin, ipmin);
            aerr = o_get_err(ilsmin, ipmin);
            if (xmin < amin)
            {
                if ((block[ilsmin][ipmin] == 1) || (block[ilsmin][ipmin] == 4))
                {
                    o_set_lens(ilsmin, ipmin, amin + aerr*.4);
                    ymin = o_chi();
                    excd[ilsmin][ipmin] *= .9;
                    excu[ilsmin][ipmin] *= .9;
                }
                else if (block[ilsmin][ipmin] == 2)
                {
                    o_set_lmax(ilsmin, ipmin, amin);
                    o_set_lmin(ilsmin, ipmin, xmin*2. - amin);
                    excd[ilsmin][ipmin] *= .9;
                }
                else if (block[ilsmin][ipmin] == 3)
                {
                    o_set_lmin(ilsmin, ipmin, xmin*2. - amin);
                    excd[ilsmin][ipmin] *= .9;
                    excu[ilsmin][ipmin] *= .9;
                };
            }
            else if (xmin > amax)
            {
                if ((block[ilsmin][ipmin] == 1) || (block[ilsmin][ipmin] == 3))
                {
                    o_set_lens(ilsmin, ipmin, amax - aerr*.4);
                    ymin = o_chi();
                    excd[ilsmin][ipmin] *= .9;
                    excu[ilsmin][ipmin] *= .9;
                }
                else if (block[ilsmin][ipmin] == 2)
                {
                    o_set_lmin(ilsmin, ipmin, amax);
                    o_set_lmax(ilsmin, ipmin, xmin*2. - amax);
                    excu[ilsmin][ipmin] *= .9;
                }
                else if (block[ilsmin][ipmin] == 4)
                {
                    o_set_lmax(ilsmin, ipmin, xmin*2. - amax);
                    excu[ilsmin][ipmin] *= .9;
                    excd[ilsmin][ipmin] *= .9;
                };
            }

            return(ymin);
        }
        else
        {
            o_set_lens(ilsmin, ipmin, x1min);
            return(y1min);
        };
    } /* end if (ipmin>=0) */

    /* spline mapping minimization
     * If min Khi2 comes from a map parameter */
    else if (imapmin >= 0)
    {
        tmp_p = map_p;
        x0 = lens[G.nlens-1].b0 = 0.;

        NPRINTF(stderr, "min_slope: %lf (%lf,%lf)\n", x0, x1min, x2min);

        if ((x0 != x1min) && (x0 != x2min))
            xmin = min_slope_par(x0, y0, x1min, y1min, x2min, y2min);
        else
            xmin = x1min;

        add_pm(map_p, G.nx, G.ny, imapmin, jmapmin, xmin);

        ymin = o_chi();
        wr_mass("mass.tmp", map_axx, map_ayy);

        if (ymin > y1min)
        {
            NPRINTF(stderr,
                    "Not a minimum: Min:%.3lf:%.3lf C:%.3lf:%.3lf M:%.3lf:%.3lf\n",
                    xmin, ymin, x0, y0, x1min, y1min);
            add_pm(map_p, G.nx, G.ny, imapmin, jmapmin, x1min);
            return(o_chi());
        };
        if (ymin < y1min)
        {
            return(ymin);
        }
        else
        {
            add_pm(map_p, G.nx, G.ny, imapmin, jmapmin, x1min);
            return(y1min);
        };
    } /* end if (imapmin>=0)*/

    /* steepest descent redshift */
    /* If min chi2 comes from a z_m_limit parameter */
    else if (izmin >= 0)
    {
        x0 = multi[izmin][0].dr;
        if ((x0 != x1min) && (x0 != x2min))
            xmin = min_slope_par(x0, y0, x1min, y1min, x2min, y2min);
        else
            xmin = x1min;

        multi[izmin][0].dr = xmin;
        ymin = o_chi();
        if (ymin > y1min)
            for (i = 0; i < 10; i++)
            {
                xmin = min_slope_par(x1min, y1min, x0, y0, xmin, ymin);
                multi[izmin][0].dr = xmin;
                ymin = o_chi();
                if (ymin < y1min)
                    i += 10;
            };
        if (ymin < y1min)
        {
            amin = zlim[izmin].ddmin;
            amax = zlim[izmin].ddmax;
            aerr = zlim[izmin].dderr;
            if (xmin < amin)
            {
                if ((zlim[izmin].bk == 1) || (zlim[izmin].bk == 4))
                {
                    multi[izmin][0].dr = amin + aerr * .4;
                    ymin = o_chi();
                    zlim[izmin].excd *= .9;
                    zlim[izmin].excu *= .9;
                }
                else if (zlim[izmin].bk == 2)
                {
                    zlim[izmin].ddmax = amin;
                    zlim[izmin].ddmin = xmin * 2. - amin;
                    zlim[izmin].excd *= .9;
                }
                else if (zlim[izmin].bk == 3)
                {
                    zlim[izmin].ddmin = xmin * 2. - amin;
                    zlim[izmin].excd *= .9;
                    zlim[izmin].excu *= .9;
                };
            }
            else if (xmin > amax)
            {
                if ((zlim[izmin].bk == 1) || (zlim[izmin].bk == 3))
                {
                    multi[izmin][0].dr = amax - aerr * .4;
                    ymin = o_chi();
                    zlim[izmin].excd *= .9;
                    zlim[izmin].excu *= .9;
                }
                else if (zlim[izmin].bk == 2)
                {
                    zlim[izmin].ddmin = amax;
                    zlim[izmin].ddmax = xmin * 2. - amax;
                    zlim[izmin].excu *= .9;
                }
                else if (zlim[izmin].bk == 4)
                {
                    zlim[izmin].ddmax = xmin * 2. - amax;
                    zlim[izmin].excd *= .9;
                    zlim[izmin].excu *= .9;
                };
            }

            return(ymin);
        }
        else
        {
            multi[izmin][0].dr = x1min;
            return(y1min);
        };

    } /* end if (izmin>=0)*/

    /*********** erreur ************/
    else
    {
        fprintf(stderr, "FATAL ERROR: arguments confusion\n");
        exit(-1);
    };
} /*end*/

