#include<stdio.h>
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
/*
 * Look for the rcut or the sigma value for the potfile that give the
 * best chi2.
 * Run the standard optimisation o_chires() function to terminate.
 *
 * This function establish a map of 1 parameter rcut or sigma that
 * is increased by the sampling value P.ircut or P.isigma at each
 * iteration.
 * At each iteration, a full standard optimisation process is run.
 *
 * Parameter :
 * - flag = 1 optimise rcut
 * - flag = 2 optimise sigma (default)
 *
 * Display more debug info if M.verbose>1.
 * Write 2 files on the disk :
 * - map.iso contains the best parameter value and the associated chi2
 * - map.res contains the chi2 and the potential parameters.
 *
 ****************************************************************/

void  o_runpot1(int flag)
{
    extern  struct  g_pot       P;
    extern  struct  g_mode      M;
    extern  struct  g_image     I;
    extern  struct  g_grille    G;
    extern struct   g_cosmo     C;
    extern  struct  ipot    ip;
    extern  struct  pot lens[], lmin[], lmax[], prec[];
    extern  struct  z_lim   zlim[];
    extern  struct  galaxie multi[][NIMAX];
    extern  int block[][NPAMAX];
    extern  double   chip, chis, chil;//,chix,chiy,
    extern  double   **map_p;

    register int    i, j, k;
    int n1, stop = 0, ic = 0; //ils=-1,ipx=-1,iz=-1
    int sblock[NLMAX][NPAMAX];
    FILE    *OUT, *MAP;
    double chi0;
    double  savchi, minchi;
    double  pp, pmin, pmax, dp;
    struct  pot best[NLMAX], save[NLMAX], slmin[NLMAX], slmax[NLMAX], sprec[NLMAX];
    struct  z_lim   szlim[NFMAX];
    double  bestz[NIMAX];

    double  **mbest_p, **smap_p;
    double  scale;  // scaling kpc/arcsec

    OUT = fopen("map.res", "w");
    MAP = fopen("map.iso", "w");

    /* sauvegarde des parametres des lentilles */
    if (lens[0].type != 10)
    {
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

    scale = d0 / C.h * distcosmo1(P.zlens);

    if (flag == 1)
    {
        if (P.ircut < 2)
            P.ircut = 2;

        n1 = P.ircut;
        pmin = P.cut1;
        pmax = P.cut2;
        NPRINTF(stderr, "INFO: rcut (kpc) : %d %.3lf %.3lf\n",
                n1, pmin*scale, pmax*scale);
    }
    else //if (flag==2)
    {
        if (P.isigma < 2)
            P.isigma = 2;

        n1 = P.isigma;
        pmin = P.sigma1;
        pmax = P.sigma2;
        NPRINTF(stderr, "INFO: sigma (km/s) : %d %.3lf %.3lf\n",
                n1, pmin, pmax);
    }

    dp = (pmax - pmin) / ((double)(n1 - 1));

    savchi = minchi = o_chi();
    NPRINTF(stderr, "INFO: minchi %lf \n", minchi);

    // sample the parameter range
    for (j = 0; j < n1; j++)
    {
        pp = pmin + ((double) j) * dp;

        if (flag == 1)
            P.cut1 = pp;
        else if (flag == 2)
            P.sigma1 = pp;

        //set_potfile(G.nplens);
        //set_lens();
        o_scale_pot();

        chi0 = o_prep();
        ic = 0;
        if (flag == 1)
        {
            NPRINTF(stderr, "\n%d/%d %.3lf %.3lf\n", j, ic, chi0, pp*scale);
        }
        else
        {
            NPRINTF(stderr, "\n%d/%d %.3lf %.3lf\n", j, ic, chi0, pp);
        }

        stop = 0;
        do
        {
            /*iterate with the pp(rcut or sigma) parameter for the potfile
             * galaxies until chi0<0 or ic>M.itMax.
             * Display more debug info if M.verbose>1. */
            chi0 = o_step(chi0);
            if (chi0 < 0.)
            {
                chi0 = -chi0;
                stop = 1;
            }

            ic++;
            NPRINTF(stderr, "%d/%d %.3lf p:%.3lf s:%.3lf l:%.3lf\r",
                    ic, M.itmax, chi0, chip, chis, chil);
        }
        while ((chi0 > M.minchi0) && (stop != 1) && (ic < M.itmax));

        if (chi0 < minchi)
        {
            minchi = chi0;

            /*Save the map parameters (rcut,sigma) that gave the best chi2*/
            if (lens[0].type != 10)
                for (i = 0; i < G.no_lens; i++)
                    best[i] = lens[i];
            else
                /*if there is a map defined in the runmode, save the best one*/
                mbest_p = map_p;

            /*for all the zm_limit sources, keep the z value that
             * gave the best chi2 */
            for (i = 0; i < I.nzlim; i++)
                bestz[i] = multi[i][0].dr;
        };

        /*print the parameter value and the associated chi2 in map.iso file*/
        if (flag == 1)
            fprintf(MAP, "%lf %lf\n", pp*scale, chi0);
        else
            fprintf(MAP, "%lf %lf\n", pp, chi0);

        /*print the chi2 and the potential parameters in the map.res file*/
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

    }; // end of sampling of the parameter range

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
    *  free maps
    */
    if (lens[0].type == 10)
    {
        free_square_double((double **) mbest_p, G.nx);
        free_square_double((double **) smap_p, G.nx);
    }
}
