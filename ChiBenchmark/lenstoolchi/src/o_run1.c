#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"

/****************************************************************/
/*      nom:        o_run1              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void  o_run1()
{
    extern  struct  g_mode      M;
    extern  struct  g_image     I;
    extern  struct  g_grille    G;
    extern  struct  ipot    ip;
    extern  struct  pot lens[], lmin[], lmax[], prec[];
    extern  struct  z_lim   zlim[];
    extern  struct  galaxie multi[][NIMAX];
    extern  int block[][NPAMAX];
    extern  double   chip, chis, chil, chia;//,chix,chiy
    extern  double   **map_p;

    register int    i, j, k;
    int ils = -1, ipx = -1, iz = -1, n1, stop = 0, ic = 0;
    int sblock[NLMAX][NPAMAX];
    FILE    *OUT, *MAP;
    double chi0;
    double  savchi, minchi;
    double  pp, pmin, pmax, dp;
    struct  pot best[NLMAX], save[NLMAX], slmin[NLMAX], slmax[NLMAX], sprec[NLMAX];
    struct  z_lim   szlim[NFMAX];
    double  bestz[NIMAX], z0[NIMAX];

    double  **mbest_p, **smap_p;

    OUT = fopen("map.res", "w");
    MAP = fopen("map.iso", "w");

    /* sauvegarde des parametres des lentilles */
    if (lens[0].type != 10)
        for (j = 0; j < G.no_lens; j++)
        {
            best[j] = lens[j];
            save[j] = lens[j];
            slmin[j] = lmin[j];
            slmax[j] = lmax[j];
            sprec[j] = prec[j];
            for (k = 0; k < ip.pmax; k++)
                sblock[j][k] = block[j][k];
        }
    else
    {
        mbest_p = (double **) alloc_square_double (G.nx, G.ny);
        mbest_p = map_p;
        smap_p = (double **) alloc_square_double (G.nx, G.ny);
        smap_p = map_p;
    }

    for (j = 0; j < I.nzlim; j++)
        szlim[j] = zlim[j];

    if (ip.lens[0] >= 0)
    {
        ils = ip.lens[0];
        ipx = ip.para[0];
        n1 = abs(block[ils][ipx]);
        pmin = o_get_lmin(ils, ipx);
        pmax = o_get_lmax(ils, ipx);
    }
    else if (ip.zlim[0] >= 0)
    {
        iz = ip.zlim[0];
        n1 = abs(zlim[iz].bk);
        pmin = zlim[iz].min;
        pmax = zlim[iz].max;
    }
    else
    {
        fprintf(stderr, "ERROR: in o_run1\n");
        exit(-1);
    };
    dp = (pmax - pmin) / ((double)(n1 - 1));

    savchi = minchi = o_chi();
    NPRINTF(stderr, "minchi %lf \n", minchi);
    NPRINTF(stderr, "boucle %d %lf %lf %lf \n", n1, pmin, pmax, dp);

    for (j = 0; j < n1; j++)
    {
        pp = pmin + ((double) j) * dp;
        if (ipx >= 0)
            o_set_lens(ils, ipx, pp);
        else
            multi[iz][0].dr = dratio(lens[0].z, pp);

        chi0 = o_prep();
        NPRINTF(stderr, "\n\n%d %d %.3lf par:%.3lf\n", j, ic = 0, chi0, pp);
        stop = 0;
        do
        {
            chi0 = o_step(chi0);
            if (chi0 < 0.)
            {
                chi0 = -chi0;
                stop = 1;
            }
            ic++;
            FPRINTF(stderr, "%d/%d %.3lf p:%.3lf s:%.3lf l:%.3lf a:%.3lf\r", ic, M.itmax, chi0,
                    chip, chis, chil, chia);
        }
        while ((chi0 > M.minchi0) && (stop != 1) && (ic < M.itmax));

        NPRINTF(stderr, "%d/%d %.3lf p:%.3lf s:%.3lf l:%.3lf a:%.3lf\n",
                ic, M.itmax, chi0, chip, chis, chil, chia);

        if (chi0 < minchi)
        {
            minchi = chi0;
            if (lens[0].type != 10)
                for (i = 0; i < G.no_lens; i++)
                    best[i] = lens[i];
            else
                mbest_p = map_p;

            for (i = 0; i < I.nzlim; i++)
            {
                bestz[i] = multi[i][0].dr;
                z0[i] = pp;
            }
        };

        fprintf(MAP, "%lf %lf\n", pp, chi0);
        o_print(OUT, chi0);

        if (lens[0].type != 10)
            for (i = 0; i < G.no_lens; i++)
            {
                if (savchi < chi0)
                    lens[i] = save[i];
                lmin[i] = slmin[i];
                lmax[i] = slmax[i];
                prec[i] = sprec[i];
                for (k = 0; k < ip.pmax; k++)
                    block[i][k] = sblock[i][k];
            }
        else
            map_p = smap_p;

        for (i = 0; i < I.nzlim; i++)
            zlim[i] = szlim[i];
    };

    fclose(MAP);
    fclose(OUT);

    if (minchi < savchi)
    {
        if (lens[0].type != 10)
            for (i = 0; i < G.no_lens; i++)
                lens[i] = best[i];
        else
            map_p = mbest_p;

        for (i = 0; i < I.nzlim; i++)
        {
            multi[i][0].dr = bestz[i];
            multi[i][0].z = z0[i];
        }
    };
    /*
    *  free maps
    */

    free_square_double((double **) mbest_p, G.nx);
    free_square_double((double **) smap_p, G.nx);
}
