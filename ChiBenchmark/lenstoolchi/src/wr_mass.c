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

void    wr_mass(char *name, double **map_xx, double **map_yy)
{
    const extern  struct  g_grille    G;

    register int    i, j;
//int   size[2];
    double  **mass;

    mass = (double **) alloc_square_double(G.nx, G.ny);

    for (i = 0; i < G.nx; i++)
        for (j = 0; j < G.ny; j++)
            mass[j][i] = map_xx[i][j] + map_yy[i][j];

    wrf_fits(name, mass, G.nx, G.ny, G.xmin, G.xmax, G.ymin, G.ymax);
    /*
    size[0]=G.nx;
    size[1]=G.ny;
    wr_ipx(name,mass,2,size,"double","bin","real","mass_map",
            G.xmin,G.xmax,G.ymin,G.ymax);
    */

    free_square_double(mass, G.nx);
}
