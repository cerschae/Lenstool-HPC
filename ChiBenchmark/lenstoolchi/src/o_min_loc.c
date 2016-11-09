#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        o_min_loc           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

double  o_min_loc(double y0)
{
    extern  struct  g_mode          M;
    extern  struct  g_grille    G;
    extern  struct  g_image I;
    extern  struct  ipot    ip;
    extern  struct  pot lens[]; //,prec[];
    extern  double  excu[][NPAMAX];
    extern  double  excd[][NPAMAX];
    extern  struct  galaxie multi[][NIMAX];
    extern  struct  z_lim   zlim[];
    extern  double  x1max, x2max, y1max, y2max;
    extern  int ipmax, ilsmax, izmax;
    extern  int imapmax, jmapmax;
    extern  double  **map_p;
    extern  double  **tmp_p;

    double  xmin, ymin;
    double  x0;
    register int    i, j;


    if (lens[0].type != 10)
    {
        if (ip.extend == 1)
        {
            NPRINTF(stderr, "\n\t\textend 1\n");
            ip.extend = 2;
            for (i = 0; i < G.no_lens; i++)
                for (j = 0; j < ip.pmax; j++)
                {
                    excu[i][j] *= 1.5;
                    excd[i][j] *= 1.5;
                };

            for (i = 0; i < I.nzlim; i++)
            {
                zlim[i].excu *= 1.5;
                zlim[i].excd *= 1.5;
            };
            return(y0);
        }
        if (ip.extend == 2)
        {
            NPRINTF(stderr, "\n\t\textend 2\n");
            ip.extend = 0;
            for (i = 0; i < G.no_lens; i++)
                for (j = 0; j < ip.pmax; j++)
                {
                    excu[i][j] *= 1.5;
                    excd[i][j] *= 1.5;
                };

            for (i = 0; i < I.nzlim; i++)
            {
                zlim[i].excu *= 1.5;
                zlim[i].excd *= 1.5;
            };
            return(y0);
        }
    }

    /* retrieve the value of x0 */

    if (ipmax >= 0)
        x0 = o_get_lens(ilsmax, ipmax);
    else if (imapmax >= 0)
        x0 = lens[G.nlens-1].b0 = 0.;
    else if (izmax >= 0)
    {
        x0 = multi[izmax][0].dr;
        NPRINTF(stderr, "a(%d) z:[%lf]  [%lf]:%lf [%lf]%lf\n",
                izmax, x0, x1max, y1max, x2max, y2max);
    }
    else
    {
        fprintf(stderr, "ERROR: arguments confusion\n");
        exit(-1);
    };

    /* find the minimum with min_parabol */

    if (((x0 > x1max + 0.001) && (x0 < x2max - 0.001)) || ((x0 < x1max - 0.001) && (x0 > x2max + .001)))
    {
        xmin = min_parabol(x0, y0, x1max, y1max, x2max, y2max);
    }
    else
    {
        NPRINTF(stderr,
                "WARNING: x0:%lf not between x1:%lf-x2:%lf l:%d p:%d z:%d\n",
                x0, x1max, x2max, ilsmax, ipmax, izmax);

        xmin = (x1max + x2max) / 2.;
    };

    if (ipmax >= 0)
    {
        o_set_lens(ilsmax, ipmax, xmin);
        excu[ilsmax][ipmax] *= .6;
        excd[ilsmax][ipmax] *= .6;
    }
    else if (imapmax >= 0)
    {
        tmp_p = map_p;
        add_pm(map_p, G.nx, G.ny, imapmax, jmapmax, xmin);
        G.exc *= .8;
    }
    else if (izmax >= 0)
    {
        multi[izmax][0].dr = xmin;
        zlim[izmax].excu *= .6;
        zlim[izmax].excd *= .6;
    };

    ymin = o_chi();

    if (ymin < y0)
        return(ymin);
    else
    {
        if (ipmax >= 0)
            o_set_lens(ilsmax, ipmax, x0);
        else if (imapmax >= 0)
            map_p = tmp_p;
        else if (izmax >= 0)
            multi[izmax][0].dr = x0;
        return(y0);
    };
}
