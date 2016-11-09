#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_runz0             */
/*      auteur:     Ghislain Golse          */
/*      date:       09/99               */
/*      place:      Toulouse            */
/****************************************************************/
/* parametres choisis sur une grille */

int optim_z;

void  o_runz0()
{
    extern  struct  g_mode      M;
    extern  struct  g_image     I;
    extern  struct  g_grille    G;
    extern  struct  ipot    ip;
    extern  struct  pot lens[]; //,lmin[],lmax[],prec[];
    extern  struct  z_lim   zlim[];
    extern  struct  galaxie multi[][NIMAX];
//extern    int block[][NPAMAX];
    extern  double   chip, chis, chil, chia; //,chix,chiy,
//extern  double   **map_p;

    register int    i, j, ipar; //k,
    int ils = -1, ipx = -1, iz = -1, n1, stop = 0, ic = 0, n_boucle, indice[NPARMAX];
    int     par_courant, fini, pas_para;
//int   sblock[NLMAX][NPAMAX];
    FILE    *OUT, *MAP;
    double   chi0, chiz, z0, z_boucle, para, parametre[NPAMAX], para_save[NPAMAX];
    double  savchi, minchi;

    double  pmin, pmax, dp;
    double  par_min[NPARMAX], par_max[NPARMAX], dpar[NPARMAX];
    struct  pot best[NLMAX];//,save[NLMAX],slmin[NLMAX],slmax[NLMAX],sprec[NLMAX];
//struct    z_lim   szlim[NFMAX];
    double  bestz[NIMAX], z_save[NIMAX];

    OUT = fopen("map.res", "w");
    MAP = fopen("map.iso", "w");


    /* determination des parametres de la grille et des bornes */
    pas_para = 10;
    o_set_map_z();

    /* determination des bornes de chaque parametre */

    for (j = 0; j < ip.map_z; j++)
    {
        ils = ip.lens_z[j];
        ipx = ip.para_z[j];
        par_min[j] = o_get_lmin(ils, ipx);
        par_max[j] = o_get_lmax(ils, ipx);
        dpar[j] = (par_max[j] - par_min[j]) / (pas_para - 1);
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
    NPRINTF(stderr, "minchi %lf \n", minchi);
    NPRINTF(stderr, "boucle %d %lf %lf %lf \n", n1, pmin, pmax, dp);

    n_boucle = 1;
    for (i = 0; i < ip.map_z; i++)
    {
        n_boucle = n_boucle * pas_para;
        indice[i] = 0;
    }
    par_courant = ip.map_z - 1;
    indice[ip.map_z-1] = -1;

    for (ipar = 0; ipar < n_boucle; ipar++) /*** boucle sur les parametres ***/
    {
        indice[par_courant]++;

        /* determination des valeurs des parametres dans l'iteration courante */

        do
        {
            if (indice[par_courant] >= pas_para)
            {
                indice[par_courant] = 0;
                par_courant--;
                /* if(par_courant==ip.map_z)
                  par_courant=0; */
                indice[par_courant]++;
                fini = 1;
            }
            else
            {
                fini = 0;
                par_courant = ip.map_z - 1;
            }
        }
        while (fini != 0);

        for (j = 0; j < ip.map_z; j++)
        {
            ils = ip.lens_z[j];
            ipx = ip.para_z[j];
            para = par_min[j] + ((double) indice[j]) * dpar[j];
            parametre[j] = para;
            o_set_lens(ils, ipx, para);
            /* NPRINTF(stderr,"parametre %d valeur %.3lf\n",ipx,para); */
        }

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
        NPRINTF(stderr,"GradPot[%d][%d]=(%.3lf,%.3lf)\n",i,j,multi[i][j].Grad.x,multi[i][j].Grad.y);
        */

        for (j = 0; j < n1; j++)
        {
            z0 = pmin + ((double) j) * dp;
            multi[iz][0].dr = dratio(lens[0].z, z0);


            chi0 = o_chi();

            /* NPRINTF(stderr,"%d/z0=%.3lf chi0=%.3lf p:%.3lf s:%.3lf l:%.3lf a:%.3lf\n",
                j,z0,chi0,chip,chis,chil,chia);
            */
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
        /* NPRINTF(stderr,"meilleur chi= %.3lf pour z= %.3lf\n\n",chiz,z_boucle); */

    };

    fclose(MAP);
    fclose(OUT);

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


    NPRINTF(stderr, "*** Meilleur ajustement : chi= %.3lf ***\n", minchi);
    for (j = 0; j < ip.map_z; j++)
        NPRINTF(stderr, "parametre %d : %.3lf\n", ip.para_z[j], para_save[j]);

    for (i = 0; i < I.nzlim; i++)
        NPRINTF(stderr, "#%d, z=%.3lf", i, multi[i][0].z);

    /**** z fixe, optimisation des parametres ****/

    optim_z = 0;
    chi0 = o_prep();
    NPRINTF(stderr, "\n\nchi0= %.3lf z= %.3lf\n", chi0, multi[i][0].z);
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


    /*
     * reset lens parameters
     */

    o_chires("chires.dat");
    set_res_par();

    /*
    *  print results
    */

    o_print_res(chi0, 0.);

    /*
    *  free maps
    */
}
