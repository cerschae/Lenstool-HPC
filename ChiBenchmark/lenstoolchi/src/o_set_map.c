#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_set_map           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/* determination des parametres qui vont varier par step        */
/* pour obtenir une carte de chi2               */
/****************************************************************/

void    o_set_map()
{
    extern  struct  g_mode          M;
    extern  struct  g_grille    G;
    extern  struct  g_image     I;
    extern  struct  pot lens[];//,lmin[],lmax[],prec[];
    extern  struct  z_lim   zlim[];
    extern  struct  ipot    ip;
    extern  int block[NLMAX][NPAMAX];
    register int    i, ils, ipx;

    ip.lens[0] = -1;
    ip.lens[1] = -1;
    ip.para[0] = -1;
    ip.para[1] = -1;
    ip.zlim[0] = -1;
    ip.zlim[1] = -1;
    if (lens[0].type != 10)
    {
        for (ils = 0; ils < G.no_lens; ils++)
            for (ipx = 0; ipx < ip.pmax; ipx++)
            {
                if (block[ils][ipx] < 0)
                {
                    if (ip.map < 2)
                    {
                        NPRINTF(stderr, "map %d %d\n", ils, ipx);
                        ip.lens[ip.map] = ils;
                        ip.para[ip.map] = ipx;
                        o_set_err(ils, ipx, 0.);
                        ip.map++;
                    }
                    else
                    {
                        fprintf(stderr, "ERROR: More than 2 mapping directions detected\n");
                        exit(-1);
                    };
                };
            };
    }
    else    /* spline map */
    {
    }

    for (i = 0; i < I.nzlim; i++)
    {
        if (zlim[i].bk < 0)
        {
            if (ip.map < 2)
            {
                NPRINTF(stderr, "map_z %d \n", i);
                ip.zlim[ip.map] = i;
                zlim[i].dderr = 0.;
                ip.map++;
            }
            else
            {
                fprintf(stderr, "ERROR: More than 2 mapping directions detected\n");
                exit(-1);
            };
        };
    };

}
