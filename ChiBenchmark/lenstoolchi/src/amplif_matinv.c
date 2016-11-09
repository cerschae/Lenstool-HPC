#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        amplif_matinv           */
/*      auteur:     Ghislain Golse          */
/*      date:       08/01               */
/*      place:      Toulouse            */
/****************************************************************
 * Calcul de la matrice inverse d'amplification pour chaque image (i,j)  *
 *
 * Global variables used :
 * - I, multi, amplifi_matinv
 * - in e_grad2() : G, lens, lens_table
 * - in distcosmo1() : C
 * - in e_amp() : G, lens, lens_table
 */
void   amplif_matinv()
{
    const extern  struct g_image  I;
    extern  struct galaxie  multi[NFMAX][NIMAX];
    extern  struct matrix   amplifi_matinv[NFMAX][NIMAX];

    struct  matrix Mat;
    double  dlsds;
    register int i, j;
    int   n;

    /* NPRINTF(stderr,"Coucou.......\n"); */

    for (i = 0; i < I.n_mult; i++)
    {
        n = I.mult[i];
        dlsds = multi[i][0].dr;
        /* NPRINTF(stderr,"n=%d dlsds=%.3lf\n",n,dlsds);*/
        for (j = 0; j < n; j++)
        {
            Mat = e_grad2_gal(&multi[i][j], NULL);
            Mat.a /= multi[i][j].dos;
            Mat.b /= multi[i][j].dos;
            Mat.c /= multi[i][j].dos;

            amplifi_matinv[i][j].a = 1. - Mat.a;
            amplifi_matinv[i][j].b = Mat.b;
            amplifi_matinv[i][j].c = 1. - Mat.c;
            amplifi_matinv[i][j].d = Mat.b;
        };
    };

}
