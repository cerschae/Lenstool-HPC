#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_set_map_z         */
/*      auteur:     Ghislain Golse          */
/*      date:       09/99               */
/*      place:      Toulouse            */
/* determination des parametres qui vont varier par o_runz_mc   */
/* pour obtenir une carte de chi2 dans le cas de                */
/* l'optimisation rapide en z                   */
/****************************************************************/

void    o_set_map_z()
{
    extern  struct  g_mode          M;
    extern  struct  g_grille    G;
    //extern    struct  g_image     I;
    extern  struct  pot lens[];//,lmin[],lmax[],prec[];
    //extern  struct    z_lim   zlim[];
    extern  struct  ipot    ip;
    extern  int block[NLMAX][NPAMAX];
    register int    ils, ipx; //i,

    ip.map_z = 0;

    if (lens[0].type != 10)
    {
        for (ils = 0; ils < G.no_lens; ils++)
            for (ipx = 0; ipx < ip.pmax; ipx++)
            {
                if (block[ils][ipx] == 1)
                {
                    NPRINTF(stderr, "map_z %d %d\n", ils, ipx);
                    ip.lens_z[ip.map_z] = ils;
                    ip.para_z[ip.map_z] = ipx;
                    ip.map_z++;
                }
            };
    }

    else    /* spline map */
    {
    }

    NPRINTF(stderr, "grid for %d parameter(s)\n", ip.map_z);

}
