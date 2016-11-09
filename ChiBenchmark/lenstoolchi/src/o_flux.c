#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        o_flux                  */
/*      auteur:     Ghislain Golse          */
/*      date:       10/99               */
/*      place:      Toulouse            */
/****************************************************************/
/* calcul de magnification pour les IMAGES MULTIPLES        */
/* calcul du flux source moyen pour une famille                 */
/****************************************************************
 * Global variables used :
 * - multi
 * - in e_grad2() : G, lens, lens_table
 */
void    o_flux(int n, double *fluxS, int n_famille, double *np_b0)
{
//  extern  struct  g_mode  M;
    extern  struct  galaxie multi[NFMAX][NIMAX];

    register int    i;
    double  A, B, C, magn0 = 0.;
    struct  matrix  MA;
    double   flux;

    flux = 0.;

    for (i = 0; i < n; i++)
    {
        multi[n_famille][i].flux = pow(10., (magn0 - multi[n_famille][i].mag) / 2.5);
        /* NPRINTF(stderr,"mag[%d]=%.3lf flux[%d]=%.3lf\n",i,multi[n_famille][i].mag,i,multi[n_famille][i].flux); */
    };
    /* NPRINTF(stderr,"\n"); */

    for (i = 0; i < n; i++)
    {
        MA = e_grad2_gal(&multi[n_famille][i], np_b0);
        MA.a /= multi[n_famille][0].dos;
        MA.b /= multi[n_famille][0].dos;
        MA.c /= multi[n_famille][0].dos;
        MA.d /= multi[n_famille][0].dos;
        A = 1. - MA.a;
        B = -MA.b;
        C = 1. - MA.c;
        multi[n_famille][i].A = fabs(A * C - B * B);
        /* NPRINTF(stderr,"A[%d]=%.3lf ",i,1./multi[n_famille][i].A); */
        flux += multi[n_famille][i].flux * multi[n_famille][i].A;
        /* NPRINTF(stderr,"flux[%d]=%.3lf ",i,multi[n_famille][i].flux*multi[n_famille][i].A); */
    }

    *fluxS = flux / n;
    /* NPRINTF(stderr,"\n fluxS=%.3lf\n",*fluxS); */

    multi[n_famille][n] = multi[n_famille][0];
}
