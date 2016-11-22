#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"
/****************************************************************/
/*      nom:        g_mass              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    wr_pot(char *name, double **map)
{
    const extern  struct  g_grille    G;
//int   size[2];
    register int    i, j;
    double  **mp;

    mp = (double **) alloc_square_double(G.nx, G.ny);

    for (i = 0; i < G.nx; i++)
        for (j = 0; j < G.ny; j++)
            mp[j][i] = map[i][j];

    wrf_fits(name, mp, G.nx, G.ny, G.xmin, G.xmax, G.ymin, G.ymax);
    /*
    size[0]=G.nx;
    size[1]=G.ny;
    wr_ipx(name,mp,2,size,"double","bin","real","poten_map",
            G.xmin,G.xmax,G.ymin,G.ymax);
    */

    free_square_double(mp, G.nx);
}