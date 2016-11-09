#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        amplif              */
/*      auteur:     Ghislain Golse          */
/*      date:       10/00               */
/*      place:      Toulouse            */
/****************************************************************/
/* Calcul de l'amplification pour chaque image (i,j)  */

/* Compute the amplification factor for each arclet and save the
 * computed values in the global array amplifi[idSource][idArclet]
 * (declared in o_global).
 *
 * Global variables used :
 * - multi, I, amplifi
 * - in e_amp() : G, lens, lens_table
 * - in distcosmo1() : C
 */
void   amplif(double *np_b0, double **amplifi)
{
    const extern struct g_image   I;
    extern struct galaxie   multi[NFMAX][NIMAX];
    //extern double amplifi[NFMAX][NIMAX];

    double A;
    int i, j;
    int   n;

    /*For each image*/
    for (i = 0; i < I.n_mult; i++)
    {
        n = I.mult[i];
        /* NPRINTF(stderr,"n=%d dlsds=%.3lf\n",n,dlsds);*/
        for (j = 0; j < n; j++)
        {
            A = 1. / e_amp_gal(&multi[i][j], np_b0);
            amplifi[i][j] = A;
            /* NPRINTF(stderr,"  A[%d][%d]=%.3lf\n ",i,j,amplifi[i][j]); */
        };
    };
}
