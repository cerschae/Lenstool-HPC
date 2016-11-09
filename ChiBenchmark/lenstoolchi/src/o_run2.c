#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"

/****************************************************************/
/*      nom:        o_run2              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void  o_run2()
{
    extern  struct  g_grille    G;
    extern  struct  g_image     I;
    extern  struct  g_mode      M;
    extern  struct  ipot    ip;
    extern  struct  pot lens[], lmin[], lmax[], prec[];
    extern  int block[][NPAMAX];
    extern  struct  z_lim   zlim[];
    extern  struct  galaxie multi[][NIMAX];
    extern  double   chip, chis, chil, chia;//chix,chiy,
    extern  double   **map_p;


    register int    i, j, k, l;
    int ils1 = -1, ipx1 = -1, n1, stop = 0, ic = 0;
    int ils2 = -1, ipx2 = -1, n2;
    int iz1 = -1, iz2 = -1;
    int sblock[NLMAX][NPAMAX];
    FILE    *OUT, *MAP;
    double chi0;
    double  savchi, minchi;
    double  p1, pmin1, pmax1, dp1;
    double  p2, pmin2, pmax2, dp2;
    struct  pot best[NLMAX], save[NLMAX], slmin[NLMAX], slmax[NLMAX], sprec[NLMAX];
    struct  z_lim   szlim[NFMAX];
    double  bestz[NIMAX];

    double   **mbest_p, **smap_p;


    OUT = fopen("map.res", "w");
    MAP = fopen("map.iso", "w");

    /* sauvegarde des parametres des lentilles */
    if (lens[0].type != 10)
        for (j = 0; j < G.no_lens; j++)
        {
            save[j] = best[j] = lens[j];
            slmin[j] = lmin[j];
            slmax[j] = lmax[j];
            sprec[j] = prec[j];
            for (k = 0; k < 10; k++)
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
        ils1 = ip.lens[0];
        ipx1 = ip.para[0];
        n1 = abs(block[ils1][ipx1]);
        pmin1 = o_get_lmin(ils1, ipx1);
        pmax1 = o_get_lmax(ils1, ipx1);
    }
    else if (ip.zlim[0] >= 0)
    {
        iz1 = ip.zlim[0];
        n1 = abs(zlim[iz1].bk);
        pmin1 = zlim[iz1].min;
        pmax1 = zlim[iz1].max;
    }
    else
    {
        fprintf(stderr, "ERROR: in o_run2\n");
        exit(-1);
    };
    if (ip.lens[1] >= 0)
    {
        ils2 = ip.lens[1];
        ipx2 = ip.para[1];
        n2 = abs(block[ils2][ipx2]);
        pmin2 = o_get_lmin(ils2, ipx2);
        pmax2 = o_get_lmax(ils2, ipx2);
    }
    else if (ip.zlim[1] >= 0)
    {
        iz2 = ip.zlim[1];
        n2 = abs(zlim[iz2].bk);
        pmin2 = zlim[iz2].min;
        pmax2 = zlim[iz2].max;
    }
    else
    {
        fprintf(stderr, "ERROR: in o_run2\n");
        exit(-1);
    };


    dp1 = (pmax1 - pmin1) / ((double)(n1 - 1));
    dp2 = (pmax2 - pmin2) / ((double)(n2 - 1));
    NPRINTF(stderr, "%d %d %d %.3lf %.3lf \n", ils1, ipx1, n1, pmin1, pmax1);
    NPRINTF(stderr, "%d %d %d %.3lf %.3lf \n", ils2, ipx2, n2, pmin2, pmax2);
    if (ipx1 >= 0)
        o_set_lens(ils1, ipx1, pmin1);
    else
        multi[iz1][0].dr = dratio(lens[0].z, pmin1);
    if (ipx2 >= 0)
        o_set_lens(ils2, ipx2, pmin2);
    else
        multi[iz2][0].dr = dratio(lens[0].z, pmin2);

    savchi = minchi = o_chi();

    fprintf(MAP, "2 %lf %lf %lf %lf\n\n%d\n%d\n", pmin2, pmax2, pmin1, pmax1, n1, n2);
    fprintf(MAP, "double txt real\nmap_optim\n");

    for (l = 0; l < n1; l++)
        for (j = 0; j < n2; j++)
        {
            ip.extend = 1;
            p1 = pmin1 + ((double) l) * dp1;
            p2 = pmin2 + ((double) j) * dp2;

            if (ipx1 >= 0)
                o_set_lens(ils1, ipx1, p1);
            else
                multi[iz1][0].dr = dratio(lens[0].z, p1);
            if (ipx2 >= 0)
                o_set_lens(ils2, ipx2, p2);
            else
                multi[iz2][0].dr = dratio(lens[0].z, p2);

            chi0 = o_prep();
            ic = 0;
            NPRINTF(stderr, "start %d %d %d %.3lf\tp1: %.3lf p2:%.3lf\n",
                    l, j, ic, chi0, p1, p2);
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
                FPRINTF(stderr, "%d/%d %.3lf p:%.3lf s:%.3lf l:%.3lf\r", ic, M.itmax, chi0,
                        chip, chis, chil);
            }
            while ((stop != 1) && (ic < M.itmax));
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
                    bestz[i] = multi[i][0].dr;
            };

            fprintf(MAP, "%lf\n", chi0);
            o_print(OUT, chi0);

            if (lens[0].type != 10)
                for (i = 0; i < G.no_lens; i++)
                {
                    if (M.minchi0 < chi0)
                        lens[i] = save[i];
                    lmin[i] = slmin[i];
                    lmax[i] = slmax[i];
                    prec[i] = sprec[i];
                    for (k = 0; k < 10; k++)
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
            multi[i][0].dr = bestz[i];
    };

    /*
    * free maps
    */

    free_square_double((double **) mbest_p, G.nx);
    free_square_double((double **) smap_p, G.nx);

}
