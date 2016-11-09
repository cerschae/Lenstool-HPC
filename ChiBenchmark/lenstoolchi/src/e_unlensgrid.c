#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        e_unlensgrid            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse
 ****************************************************************
 * Fill the gsource global variable at the redshift corresponding to
 * dlsds so that it establishes a bijection with the grid gimage
 * global variable in the image plane.
 *
 * gsource is a 2D map of G.ngrid^2 points in the source plane defined
 * by dlsds. Each point in gsource is linked to a single point in the
 * gimage global variable. Those 2 maps define a kind of bijection
 * between the source and the image planes.
 */
void    e_unlensgrid(struct point gsource[][NGGMAX], double dlsds)
{
    const extern  struct  g_grille    G;
    const extern  struct  point   gimage[NGGMAX][NGGMAX];

    int    i, j;

    for (i = 0; i < G.ngrid; i++)
        for (j = 0; j < G.ngrid; j++)
            e_dpl(&gimage[i][j], dlsds, &gsource[i][j]);
}


