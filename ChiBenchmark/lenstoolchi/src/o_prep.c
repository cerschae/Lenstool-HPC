#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_prep              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

double  o_prep()
{
    extern  struct  g_mode          M;
    extern  struct  g_grille    G;
    extern  struct  g_image     I;
    extern  struct  z_lim   zlim[];
    extern  struct  pot lens[]; //,lmin[],lmax[],prec[];
    extern  int     block[NLMAX][NPAMAX];
    extern  double  excu[NLMAX][NPAMAX];
    extern  double  excd[NLMAX][NPAMAX];
    extern  struct sigposStr sigposAs;

    double  chi0;
    register int    i, j;

    /* chargement des parametres initiaux de l'excursion */
    // for the optimised clumps
    if (lens[0].type != 10)
    {
        for (i = 0; i < G.no_lens; i++)
        {
            o_set_exc(i, excu, excd, block);
            o_set_start(i, block);
        }
    }
    else
    {
        NPRINTF(stderr, "exc:%lf min:%lf\n", G.exc, G.excmin);
        G.no_lens = 1;
        G.nlens++;
        lens[G.nlens-1].type = 11;
        lens[G.nlens-1].b0 = 0.;
        lens[G.nlens-1].rc = 2.*G.echant * G.dx;
    }

    // for the images redshift
    for (i = 0; i < I.nzlim; i++)
    {
        zlim[i].excu = .3;
        zlim[i].excd = .3;
    }

    // for the sigposAs
    sigposAs.excu = .3;
    sigposAs.excd = .3;
    for ( i = 0; i < I.n_mult; i++)
        for ( j = 0; j < I.mult[i]; j++)
            I.sig2pos[i][j] = sigposAs.min * sigposAs.min;

    chi0 = o_chi();

    return(chi0);
}
