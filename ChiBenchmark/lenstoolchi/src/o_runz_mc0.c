#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_runz_mc0          */
/*      auteur:     Ghislain Golse          */
/*      date:       09/99               */
/*      place:      Toulouse            */
/****************************************************************/
/* generation de parametres avec la methode de Monte Carlo */

void  o_runz_mc0()
{
    extern  struct  g_mode      M;
    extern  struct  g_image     I;
    extern  struct  g_grille    G;
    extern  struct  ipot    ip;
    extern  struct  MCarlo      mc;
    extern  struct  pot lens[]; //,lmin[],lmax[],prec[];
    extern  struct  z_lim   zlim[];
    extern  struct  galaxie multi[NFMAX][NIMAX];
//extern    int block[][NPAMAX];
    extern  double   chip, chis, chil, chia; //chix,chiy,
//extern  double   **map_p;
    extern  int      optim_z;

    register int    i, j, k, ipar, i_percent;
    int ils = -1, ipx = -1, iz = -1, n1, stop = 0, fini, ic = 0, n_cases, n_MC;
    int     npar_MC, nOK, n_over;//,nMC_Max;
    int     iboucle, iboucle0, n_effectif;
    char    ichar;
//int   sblock[NLMAX][NPAMAX];
    FILE    *OUT, *MAP;
    double   chi0, chiz, z0, z_boucle, para, para_rand, parametre[NPAMAX], para_save[NPAMAX];
    double   chi_min, chi_max;
    double  savchi, minchi, chi_total, chi_moyen, chi_case;

    double  pmin, pmax, dp, zmin, zmax, percent;
    double  par_min[NPARMAX], par_max[NPARMAX], dpar[NPARMAX];
    struct  pot best[NLMAX], save[NLMAX];//,slmin[NLMAX],slmax[NLMAX],sprec[NLMAX];
//struct    z_lim   szlim[NFMAX];
    double  bestz[NIMAX], z_save[NIMAX];



    FILE **chisq, *nMC;
    char concat[15], dat[] = ".dat";

    struct Chi_MC
    {
        int n;
        int n100;
        int OK;
        double *chi2;
    }
    *chi_MC;

    int position[NPARMAX], pos;

    OUT = fopen("map.res", "w");
    MAP = fopen("map.iso", "w");
    nMC = fopen("nMC.dat", "w");

    chisq = (FILE**)malloc(mc.iterations * sizeof(*chisq));

    /* determination des parametres de la grille et des bornes */

    o_set_map_z();

    /* determination des bornes de chaque parametre */

    npar_MC = mc.squares_par;
    n_cases = (long int)pow(((double)npar_MC), ((double)ip.map_z));
    n_MC = n_cases * mc.tosses_sq;
    chi_MC = (struct Chi_MC*)malloc(n_cases * sizeof(*chi_MC));

    NPRINTF(stderr, "n squares=%d^%d=%d n_MC=%d\n", npar_MC, ip.map_z, n_cases, n_MC);

    for (j = 0; j < n_cases; j++)
    {
        chi_MC[j].OK = 1;
        chi_MC[j].chi2 = (double*)malloc(5 * mc.tosses_sq * sizeof(double));
    };

    for (j = 0; j < ip.map_z; j++)
    {
        ils = ip.lens_z[j];
        ipx = ip.para_z[j];
        par_min[j] = o_get_lmin(ils, ipx);
        par_max[j] = o_get_lmax(ils, ipx);
        dpar[j] = (par_max[j] - par_min[j]) / ((double)npar_MC);
    }

    /* bornes en z */

    if (ip.zlim[0] >= 0)
    {
        iz = ip.zlim[0];
        n1 = abs(zlim[iz].bk);
        pmin = zlim[iz].min;
        pmax = zlim[iz].max;
        dp = (pmax - pmin) / ((double)(n1 - 1));
    }
    else
    {
        fprintf(stderr, "ERROR: in o_runz\n");
        exit(-1);
    };

    savchi = minchi = o_chi();
    /* NPRINTF(stderr,"dr=%.3lf z0=%.3lf chi0=%.3lf",multi[iz][0].dr,multi[iz][0].z,minchi); */

    NPRINTF(stderr, "minchi %lf \n", minchi);

    if (lens[0].type != 10)
    {
        for (i = 0; i < G.no_lens; i++)
            best[i] = lens[i];

        for (i = 0; i < ip.map_z; i++)
        {
            ils = ip.lens_z[i];
            ipx = ip.para_z[i];
            para_save[i] = o_get_lens(ils, ipx);
            /* NPRINTF(stderr,"parameter %d : %.3lf\n",ip.para[i],para_save[i]); */
        };
    }

    NPRINTF(stderr, "boucle %d %lf %lf %lf \n", n1, pmin, pmax, dp);

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
        };
        i_percent = 1;
        chi_min = 1000000.;
        chi_max = 0.;

        for (ipar = 0; ipar < n_MC; ipar++) /*** boucle sur les parametres ***/
        {
            n_effectif++;
            /*** determination des valeurs des parametres, au hasard ***/
            fprintf(nMC, "%.3lf\n", chi_total);

            for (j = 0; j < ip.map_z; j++)
            {
                ils = ip.lens_z[j];
                ipx = ip.para_z[j];
                para_rand = (rand() % 1000) / 1000.*(par_max[j] - par_min[j]);
                para = par_min[j] + para_rand;
                position[j] = floor(para_rand / dpar[j]);
                parametre[j] = para;
                o_set_lens(ils, ipx, para);
                /* NPRINTF(stderr,"parametre %d valeur %.3lf\n",ipx,para); */
            }

            pos = 0;
            for (i = 0; i < ip.map_z; i++)
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


            if ((chi_MC[pos].OK == 1) && (chi_MC[pos].n100 == 0)) /** calcul du chi si on est dans la bonne zone **/
            {
                chi_MC[pos].n++;
                /* NPRINTF(stderr,"position= %d n=%d\n",pos,chi_MC[pos].n); */

                /* z0 fixe pour trouver le gradient*/
                optim_z = 0;

                z0 = pmin;

                multi[iz][0].dr = dratio(lens[0].z, z0);

                chi0 = o_chi();
                /* NPRINTF(stderr," chi0=%.3lf z0=%.3lf\n",chi0,z0); */
                stop = 0;

                /******* boucle en z ****************/
                optim_z = 1;
                chiz = 100000.;

                /* for(i=0;i<I.n_mult;i++)
                   for(j=0;j<I.mult[i];j++)
                   NPRINTF(stderr,"GradPot[%d][%d]=(%.3lf,%.3lf)\n",i,j,multi[i][j].Grad.x,multi[i][j].Grad.y); */


                for (j = 0; j < n1; j++)
                {
                    z0 = pmin + ((double) j) * dp;
                    multi[iz][0].dr = dratio(lens[0].z, z0);

                    chi0 = o_chi();
                    /* NPRINTF(stderr,"dr=%.3lf z0=%.3lf chi0=%.3lf\n",multi[iz][0].dr,z0,chi0); */
                    /* NPRINTF(stderr,"%d/z0=%.3lf chi0=%.3lf p:%.3lf s:%.3lf l:%.3lf a:%.3lf\n",j,z0,chi0,chip,chis,chil,chia); */

                    /* for(i=0;i<I.n_mult;i++)
                       for(j=0;j<I.mult[i];j++)
                       NPRINTF(stderr,"GradPot[%d][%d]=(%.3lf,%.3lf)\n",i,j,multi[i][j].Grad.x,multi[i][j].Grad.y); */

                    if (chi0 < minchi)
                    {
                        minchi = chi0;
                        if (lens[0].type != 10)
                        {
                            for (i = 0; i < G.no_lens; i++)
                                best[i] = lens[i];

                            for (i = 0; i < ip.map_z; i++)
                                para_save[i] = parametre[i];
                        }
                        for (i = 0; i < I.nzlim; i++)
                        {
                            bestz[i] = multi[i][0].dr;
                            z_save[i] = z0;
                            /* for(i=0;i<I.nzlim;i++)
                               NPRINTF(stderr,"chi=%.3lf  z =%.3lf\n",chi0,z_save[i]); */
                        }
                    };

                    if (chi0 < chiz)  /* meilleur z de la boucle */
                    {
                        chiz = chi0;
                        z_boucle = z0;
                    };

                    fprintf(MAP, "%lf %lf\n", z0, chi0);
                    o_print(OUT, chi0);

                }
                if (chiz < chi_min)
                    chi_min = chiz;
                if (chiz > chi_max)
                    chi_max = chiz;

                fprintf(chisq[iboucle], "%.3lf\n", chiz);

                chi_MC[pos].chi2[chi_MC[pos].n] = chiz;
                chi_total = chi_total + chiz;
                /* NPRINTF(stderr,"#%d meilleur chi= %.3lf pour z= %.3lf\n\n",ipar,chiz,z_boucle); */
            }
            else
                ipar--;

            if (((double)ipar / n_MC*10) >= i_percent)
            {
                NPRINTF(stderr, "%d \n", 10*i_percent);
                i_percent++;
            }

            if (((double)n_over / n_MC) > 0.2)
            {
                ipar = n_MC;
            }
        };

        chi_total = chi_total / n_MC;

        NPRINTF(stderr, "#%d chi2: min= %.3lf mean=%.3lf max=%.3lf\n", iboucle, chi_min, chi_total, chi_max);

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

        n_MC = fini * mc.tosses_sq;

        NPRINTF(stderr, "effective n=%d n squares=%d n_MC=%d\n\n", n_effectif, fini, n_MC);
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

        for (i = 0; i < I.nzlim; i++)
        {
            multi[i][0].dr = bestz[i];
            multi[i][0].z = z_save[i];
        }
    };


    NPRINTF(stderr, "\n*** Best fit : chi= %.3lf ***\n", minchi);
    for (j = 0; j < ip.map_z; j++)
        NPRINTF(stderr, "parameter %d : %.3lf\n", ip.para_z[j], para_save[j]);

    for (i = 0; i < I.nzlim; i++)
        NPRINTF(stderr, "#%s, z=%.3lf\n\n", multi[i][0].n, multi[i][0].z);

    /**** z fixe a +/- percent%, optimisation des parametres ****/
    savchi = minchi;
    if (lens[0].type != 10)
        for (i = 0; i < G.no_lens; i++)
            save[i] = lens[i];

    optim_z = 0;
    z0 = multi[iz][0].z;
    percent = zlim[iz].percent / 100.;
    zmin = (1. - percent) * z0;
    zmax = (1. + percent) * z0;
    n1 = zlim[iz].bk0;
    dp = (zmax - zmin) / (n1 - 1);

    for (j = 0; j < n1; j++)
    {
        z0 = zmin + ((double) j) * dp;
        multi[iz][0].z = z0;
        multi[iz][0].dr = dratio(lens[0].z, z0);
        if (lens[0].type != 10)
            for (i = 0; i < G.no_lens; i++)
                lens[i] = save[i];
        chi0 = o_prep();
        /* NPRINTF(stderr,"chi0= %.3lf z= %.3lf\n",chi0,multi[iz][0].z); */
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

        NPRINTF(stderr, "%d/%d z=%.3lf chi=%.3lf p:%.3lf s:%.3lf l:%.3lf a:%.3lf\n",
                ic, M.itmax, z0, chi0, chip, chis, chil, chia);

        if (chi0 < minchi)
        {
            minchi = chi0;
            for (i = 0; i < G.no_lens; i++)
                best[i] = lens[i];
            for (i = 0; i < I.nzlim; i++)
            {
                bestz[i] = multi[i][0].dr;
                z_save[i] = multi[i][0].z;
            }
        };
    };

    if (lens[0].type != 10)
        for (i = 0; i < G.no_lens; i++)
            lens[i] = best[i];

    for (i = 0; i < I.nzlim; i++)
    {
        multi[i][0].dr = bestz[i];
        multi[i][0].z = z_save[i];
    }


    NPRINTF(stderr, "\n**** z optimisation ****\n");
    NPRINTF(stderr, "chi2= %.3lf\n", minchi);
    for (i = 0; i < I.nzlim; i++)
        NPRINTF(stderr, "#%s, z=%.3lf\n", multi[i][0].n, multi[i][0].z);

    /*
    *  free maps
    */

}

