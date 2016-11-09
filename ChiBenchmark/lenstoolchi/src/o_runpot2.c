#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"

/****************************************************************/
/*      nom:        o_runpot2               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/*
 * Look for the rcut and sigma values for the potfile that give the
 * best chi2.
 * Run the standard optimisation o_chires() function to terminate.
 *
 * This function establish a map of 2 parameters rcut and sigma that
 * are increased by the sampling value P.ircut or P.isigma at each
 * iteration.
 * At each iteration, a full standard optimisation process is run.
 *
 * Default behavior :
 * If o_runpot2() is called but ircut or isigma are < 2 then the default
 * behavior is to set ircut=2 and isigma=2.
 *
 * Display more debug info if M.verbose>1.
 * Write 2 files on the disk :
 * - map.iso contains the best chi2.
 * - map.res contains the chi2 and the potential parameters.
 *
 ****************************************************************/

void  o_runpot2()
{
    extern struct g_grille G;
    extern struct g_image  I;
    extern struct g_mode   M;
    extern struct g_pot    P;
    extern struct g_cosmo  C;
    extern struct ipot     ip;
    extern struct pot        lens[], lmin[], lmax[], prec[];
    extern int    block[][NPAMAX];
    extern struct z_lim    zlim[];
    extern struct galaxie  multi[][NIMAX];
    extern double chip, chis, chil; //,chix,chiy,
    extern double **map_p;

    register int  i, j, k, l;
    int    n1, n2;
    int    stop = 0, ic = 0;
//  int  iz1=-1,iz2=-1;
    int    sblock[NLMAX][NPAMAX];
    FILE     *OUT, *MAP;
    double chi0;
    double savchi, minchi;
    double p1, pmin1, pmax1, dp1;
    double p2, pmin2, pmax2, dp2;
    struct pot   best[NLMAX], save[NLMAX], slmin[NLMAX], slmax[NLMAX], sprec[NLMAX];
    struct z_lim szlim[NFMAX];
    double bestz[NIMAX];

    double **mbest_p, **smap_p;
    double scale;   // scaling kpc/arcsec

    OUT = fopen("map.res", "w");
    MAP = fopen("map.iso", "w");

    /* sauvegarde des parametres des lentilles */
    if (lens[0].type != 10)
    {
        for (j = 0; j < G.no_lens; j++)
        {
            save[j] = best[j] = lens[j];
            slmin[j] = lmin[j];
            slmax[j] = lmax[j];
            sprec[j] = prec[j];
            for (k = 0; k < 10; k++)
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


    if (P.ircut < 2)
        P.ircut = 2;

    if (P.isigma < 2)
        P.isigma = 2;

    n1 = P.ircut;
    pmin1 = P.cut1;
    pmax1 = P.cut2;
    scale = d0 / C.h * distcosmo1(P.zlens);

    n2 = P.isigma;
    pmin2 = P.sigma1;
    pmax2 = P.sigma2;

    dp1 = (pmax1 - pmin1) / ((double)(n1 - 1));
    dp2 = (pmax2 - pmin2) / ((double)(n2 - 1));
    NPRINTF(stderr, "INFO: rcut (kpc) : %d %.3lf %.3lf \n", n1, pmin1*scale, pmax1*scale);
    NPRINTF(stderr, "INFO: sigma (km/s) : %d %.3lf %.3lf \n", n2, pmin2, pmax2);

    savchi = minchi = o_chi();

    fprintf(MAP, "2 sigmin:%lf sigmax:%lf rcutmin:%lf rcutmax:%lf\n\n",
            pmin2, pmax2, pmin1*scale, pmax1*scale);
    fprintf(MAP, "Nrcut:%d\nNsig:%d\n", n1, n2);
    fprintf(MAP, "double txt real\nmap_optim\n");

    for (l = 0; l < n1; l++)
        for (j = 0; j < n2; j++)
        {
            ip.extend = 1;
            p1 = pmin1 + ((double) l) * dp1;
            p2 = pmin2 + ((double) j) * dp2;

            // reset the clumps parameters
            P.cut = p1;
            P.sigma = p2;

            //set_potfile(G.nplens);
            //set_lens();
            o_scale_pot();

            chi0 = o_prep();
            ic = 0;
            NPRINTF(stderr, "start %d %d %d %.3lf p1:%.3lf p2:%.3lf\n",
                    l, j, ic, chi0, p1*scale, p2);
            stop = 0;

            /*iterate with the p1(rcut) and p2(sigma) parameters for the
             * potfile galaxies until chi0<0 or ic>M.itMax.
             * Display more debug info if M.verbose>1. */
            do
            {
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
            while ((stop != 1) && (ic < M.itmax));

            // Print last iteration ==> best chi2 for p1 and p2
            NPRINTF(stderr, "%d/%d %.3lf p:%.3lf s:%.3lf l:%.3lf\n",
                    ic, M.itmax, chi0, chip, chis, chil);

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
            fprintf(MAP, "%lf %lf %lf\n", p1*scale, p2, chi0);
            /*print the chi2 and the potential parameters in the map.res file*/
            fprintf(OUT, "%d/%d : %lf %lf\n", ic, M.itmax, p1*scale, p2);
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

        }; // end of loop over p1 and p2

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
