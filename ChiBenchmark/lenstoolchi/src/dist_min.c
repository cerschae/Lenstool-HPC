#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        dist_min                */
/*      auteur:     Ghislain Golse          */
/*      date:       10/99               */
/*      place:      Toulouse            */
/****************************************************************/
/* Calcule la distance la plus courte entre 2 images d'une meme famille
 *
 * Pour chaque famille de n images, calcule les n-1 distances euclidiennes
 * entre 2 images et sauvegarde la distance la plus courte dans la liste
 * globale distmin[NFMAX].
 *
 * Call by o_global() if I.n_mult>0.
 * */

void   dist_min()
{
    extern  double  distmin[NFMAX];

    const extern  struct  g_image I;
    const extern  struct  g_mode  M;
    const extern  struct  galaxie     multi[NFMAX][NIMAX];

    int i, j, k;
    double d, d_min;

    NPRINTF(stderr, "INFO: Minimal distance between arcs for each familly\n");
    for (i = 0; i < I.n_mult; i++)
    {
        d_min = 1000.;
        for (j = 0; j < I.mult[i]; j++)
        {
            for (k = j + 1; k < I.mult[i]; k++)
            {
                d = dist(multi[i][j].C, multi[i][k].C);
                NPRINTF(stderr, "\tdist(%d:%d,%d):%.3lf\n", i, j, k, d);

                if (d < d_min)
                    d_min = d;
            };
        };
        if (d_min > 15.)
            d_min = 15.;

        distmin[i] = d_min;

        NPRINTF(stderr, "\tdistmin[%d]=%.3lf\n", i, distmin[i]);
    };

}

