#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_run_mc            */
/*      auteur:     Ghislain Golse          */
/*      date:       10/00               */
/*      place:      Toulouse            */
/****************************************************************/
/* generation de parametres avec la methode de Monte Carlo
pour optimiser les parametres du potentiel */


void  o_run_mc()
{
    extern  struct  g_mode      M;
//  extern  struct  g_image     I;
    extern  struct  g_grille    G;
    extern  struct  ipot    ip;
    extern  struct  MCarlo      mc;
    extern  struct  pot lens[];//,lmin[],lmax[],prec[];
//  extern  struct  galaxie multi[NFMAX][NIMAX];
//  extern  int block[][NPAMAX];
    extern  double   chip, chis, chil, chia; //,chix,chiy,
//  extern  double   **map_p;
    extern  int      optim_z;

    register int    i, k, ipar, i_percent;
    int ils = -1, ipx = -1, stop = 0, fini, fini0, ic = 0, n_MC;//iz=-1,n1,
    long    j, n_cases; //,nMC_Max;
    int     n_over, npar_MC, n_it, nOK;
    int     n_effectif, iboucle, iboucle0;
    char    ichar;
//  int sblock[NLMAX][NPAMAX];
    FILE    *OUT, *MAP;
    double   chi0, para, para_rand, parametre[NPAMAX], para_save[NPAMAX];
    double   chi_min, chi_max;
    double  savchi, minchi, chi_total, chi_moyen, chi_case;
    double  nsqrt(double, int);


    double  par_min[NPARMAX], par_max[NPARMAX], dpar[NPARMAX];
    struct  pot best[NLMAX], save[NLMAX]; //,slmin[NLMAX],slmax[NLMAX],sprec[NLMAX];

    FILE **chisq, *nMC;
    char concat[15], dat[] = ".dat";

    struct Chi_MC
    {
        int n;
        int n100;
        int OK;
        double *chi2;
    } *chi_MC;

    int position[NPARMAX], pos;

    OUT = fopen("map.res", "w");
    MAP = fopen("map.iso", "w");
    nMC = fopen("nMC.dat", "w");

    chisq = (FILE**)malloc(mc.iterations * sizeof(*chisq));

    optim_z = 0;

    /* determination des parametres de la grille et des bornes */

    o_set_map_mc();

    /* determination des bornes de chaque parametre */

    n_MC = mc.n_MonteCarlo;
    npar_MC = mc.squares_par;
    n_cases = (long int)pow(((double)npar_MC), ((double)ip.map));
    n_MC = n_cases * mc.tosses_sq;
    chi_MC = (struct Chi_MC*)malloc(n_cases * sizeof(*chi_MC));
    fini0 = n_cases;

    NPRINTF(stderr, "n squares=%d^%d=%li n_MC=%d\n", npar_MC, ip.map, n_cases, n_MC);

    /* for(i=0;i<10;i++)
    {
    n_cases=1700000+100000*i;
    NPRINTF(stderr,"n cases=%d\n",n_cases);
    chi_MC=malloc(n_cases*sizeof(*chi_MC)); */

    for (j = 0; j < n_cases; j++)
    {
        chi_MC[j].OK = 1;
    };

    /* free(chi_MC);
      }; */

    for (j = 0; j < ip.map; j++)
    {
        ils = ip.lens[j];
        ipx = ip.para[j];
        par_min[j] = o_get_lmin(ils, ipx);
        par_max[j] = o_get_lmax(ils, ipx);
        dpar[j] = (par_max[j] - par_min[j]) / ((double)npar_MC);
    }

    savchi = minchi = o_chi();
    /* NPRINTF(stderr,"dr=%.3lf z0=%.3lf chi0=%.3lf",multi[iz][0].dr,multi[iz][0].z,minchi); */

    NPRINTF(stderr, "minchi %lf \n", minchi);

    if (lens[0].type != 10)
    {
        for (i = 0; i < G.no_lens; i++)
            best[i] = lens[i];

        for (i = 0; i < ip.map; i++)
        {
            ils = ip.lens[i];
            ipx = ip.para[i];
            para_save[i] = o_get_lens(ils, ipx);
            /* NPRINTF(stderr,"parameter %d : %.3lf\n",ip.para[i],para_save[i]); */
        };
    }

    chi_moyen = 1000000.;

    for (iboucle = 0; iboucle < mc.iterations; iboucle++)
    {
        strcpy(concat, "chisq");
        ichar = iboucle + 48;
        strcat(concat, &ichar);
        strcat(concat, dat);
        chisq[iboucle] = fopen(concat, "w");
        n_effectif = 0;
        chi_total = 0.;
        n_over = 0;
        fprintf(nMC, "%d ", n_MC);
        for (j = 0; j < n_cases; j++)
        {
            chi_MC[j].n100 = 0;
            chi_MC[j].n = 0;
            chi_MC[j].chi2 = (double*)malloc(5 * mc.tosses_sq * sizeof(double));
        };
        i_percent = 1;
        chi_min = 1000000.;
        chi_max = 0.;

        for (ipar = 0; ipar < n_MC; ipar++) /*** boucle sur les parametres ***/
        {
            n_effectif++;
            /*** determination des valeurs des parametres, au hasard ***/

            for (j = 0; j < ip.map; j++)
            {
                ils = ip.lens[j];
                ipx = ip.para[j];
                para_rand = (rand() % 1000) / 1000.*(par_max[j] - par_min[j]);
                para = par_min[j] + para_rand;
                position[j] = floor(para_rand / dpar[j]);
                /* srand(position[j]); */
                parametre[j] = para;
                o_set_lens(ils, ipx, para);
                /* NPRINTF(stderr,"parametre %d valeur %.3lf\n",ipx,para); */
            }

            pos = 0;

            for (i = 0; i < ip.map; i++)
                pos = pos + position[i] * pow(((double)npar_MC), ((double)i));

            if (chi_MC[pos].n >= 5*mc.tosses_sq - 1)
            {
                n_over++;
                /* NPRINTF(stderr,"chi_MC[%d].n=%d, chi_MC[%d].n100=%d \n",pos,chi_MC[pos].n,pos,chi_MC[pos].n100); */
                chi_MC[pos].n100++;
                if (chi_MC[pos].n100 >= 1)
                    srand(chi_MC[pos].n100);

                ipar--;
            }

            /* calcul du chi si on est dans la bonne zone et s'il n'y a pas eu trop de
             * tirages dans cette case**/
            if ((chi_MC[pos].OK == 1) && (chi_MC[pos].n100 == 0))
            {
                chi_MC[pos].n++;
                /* NPRINTF(stderr,"chi_MC[%d].n=%d\n",pos,chi_MC[pos].n); */
                chi0 = o_chi();
                /* NPRINTF(stderr," chi0=%.3lf\n",chi0); */
                stop = 0;

                if (chi0 < minchi)
                {
                    minchi = chi0;
                    if (lens[0].type != 10)
                    {
                        for (i = 0; i < G.no_lens; i++)
                            best[i] = lens[i];

                        for (i = 0; i < ip.map; i++)
                        {
                            para_save[i] = parametre[i];
                            /* NPRINTF(stderr,"parameter %d : %.3lf\n",ip.para[i],para_save[i]); */
                        };
                    }
                };

                if (chi0 < chi_min)
                    chi_min = chi0;

                if (chi0 > chi_max)
                    chi_max = chi0;

                o_print(OUT, chi0);
                fprintf(chisq[iboucle], "%.3lf\n", chi0);
                chi_MC[pos].chi2[chi_MC[pos].n] = chi0;
                chi_total = chi_total + chi0 ;
            }
            else
                ipar--;

            if (((double)ipar / n_MC*10) >= i_percent)
            {
                NPRINTF(stderr, "%d  ", 10*i_percent);
                i_percent++;
            }

            n_it = n_MC;
            /* NPRINTF(stderr,"iboucle=%d ipar=%d n effectif=%d \n ",iboucle,ipar, n_effectif); */

            if (((double)n_over / n_MC) > 0.2)
            {
                n_it = ipar;
                ipar = n_MC;
            }

        };

        chi_total = chi_total / n_it;

        NPRINTF(stderr, "\n#%d chi2: min=%.3lf mean=%.3lf max=%.3lf\n", iboucle, chi_min, chi_total, chi_max);

        fprintf(nMC, "%.3lf\n", chi_total);
        if (chi_total > chi_moyen)
            iboucle = mc.iterations;

        chi_moyen = chi_total;

        /*** calcul de la zone de chi minimum ***/

        fini = 0;
        nOK = 0;
        for (j = 0; j < n_cases; j++)
        {
            chi_case = 0.;
            /* NPRINTF(stderr,"case %d nombre de chi : %d\n",j,chi_MC[j].n); */

            if (chi_MC[j].n != 0 && chi_MC[j].OK == 1)
            {
                for (k = 0; k < chi_MC[j].n; k++)
                    chi_case = chi_case + chi_MC[j].chi2[k];

                chi_case = chi_case / chi_MC[j].n;
                if (chi_case > 0.5*chi_total)
                {
                    chi_MC[j].OK = 0;
                    nOK++;
                }
                else
                    fini++;
            }
        };

        if (fini != 0)
        {
            mc.tosses_sq = (double)fini0 / fini * mc.tosses_sq;
            n_MC = (double)fini0 / fini * n_MC;
            fini0 = fini;

            NPRINTF(stderr, "effective n=%d n squares=%d n_MC=%d tosses per square=%d\n\n", n_effectif, fini, n_MC, mc.tosses_sq);
            /* NPRINTF(stderr,"fini=%d n_over=%d n_over/n_MC=%lf\n",fini,n_over,((double)n_over/n_MC)); */
        }
        else
        {
            NPRINTF(stderr, "effective n=%d n squares=%d\n\n", n_effectif, fini);
            /* NPRINTF(stderr,"fini=%d n_over=%d n_over/n_MC=%lf\n",fini,n_over,((double)n_over/n_MC)); */
        }

        iboucle0 = iboucle;

        if (fini <= 2 || ((double)n_over / n_MC) > 0.2 || ((double)fini / n_cases) < 0.1  )
            iboucle = mc.iterations;

        fclose(chisq[iboucle0]);
    }

    fclose(MAP);
    fclose(OUT);
    fclose(nMC);
    free(chi_MC);

    if (minchi < savchi)
    {
        if (lens[0].type != 10)
            for (i = 0; i < G.no_lens; i++)
                lens[i] = best[i];
    };


    NPRINTF(stderr, "\n*** Best fit : chi= %.3lf ***\n", minchi);
    for (j = 0; j < ip.map; j++)
        NPRINTF(stderr, "parameter %d : %.3lf\n", ip.para[j], para_save[j]);


    /**** optimisation des parametres ****/

    savchi = minchi;
    if (lens[0].type != 10)
        for (i = 0; i < G.no_lens; i++)
            save[i] = lens[i];

    if (lens[0].type != 10)
        for (i = 0; i < G.no_lens; i++)
            lens[i] = save[i];

    chi0 = o_prep();
    stop = 0;
    ic = 0;
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

    NPRINTF(stderr, "%d/%d chi=%.3lf p:%.3lf s:%.3lf l:%.3lf a:%.3lf\n",
            ic, M.itmax, chi0, chip, chis, chil, chia);

    if (chi0 < minchi)
    {
        minchi = chi0;
        for (i = 0; i < G.no_lens; i++)
            best[i] = lens[i];
    };


    if (lens[0].type != 10)
        for (i = 0; i < G.no_lens; i++)
            lens[i] = best[i];

    NPRINTF(stderr, "\n****  optimisation ****\n");
    NPRINTF(stderr, "chi2= %.3lf\n", minchi);

    /*
    *  free maps
    */

}



