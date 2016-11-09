#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        amplif_mat              */
/*      auteur:     Ghislain Golse          */
/*      date:       08/01               */
/*      place:      Toulouse            */
/****************************************************************
 * Calcul de la matrice d'amplification pour chaque image (i,j)  *
 *
 * Global variables used :
 * - I, multi, amplifi_mat
 * - in e_grad2() : G, lens, lens_table
 * - in distcosmo1() : C
 * - in e_amp() : G, lens, lens_table
 */
void   amplif_mat()
{
    const extern struct   g_image I;
    extern struct   galaxie multi[NFMAX][NIMAX];
    extern struct   matrix  amplifi_mat[NFMAX][NIMAX];

    struct  matrix Mat;
    double  A;
    register int i, j;

    for (i = 0; i < I.n_mult; i++)
        for (j = 0; j < I.mult[i]; j++)
        {
            Mat = e_grad2_gal(&multi[i][j], NULL);
            Mat.a /= multi[i][j].dos;
            Mat.b /= multi[i][j].dos;
            Mat.c /= multi[i][j].dos;

            // Prevent infinite amplification issues
            // and makes sure amplification is bounded by MAXAMPINV or -MAXAMPINV
#define MAXAMPINV 0.001
            A = (1. - Mat.a) * (1. - Mat.c) - Mat.b * Mat.b;
//            if( fabs(A) < MAXAMPINV )
//            {
//                double k, g2, koffset, d;
//                k = 1. - 0.5 * (Mat.a + Mat.c); // 1-k
//                g2 = 0.5 * (Mat.a - Mat.c);
//                g2 = g2 * g2 + Mat.b * Mat.b;
//                if( A < 0. && g2 >= MAXAMPINV ) 
//                    d = A + MAXAMPINV;
//                else
//                    d = A - MAXAMPINV;
//  
//                if( k * k - d < 0. ) exit(1);
//  
//                koffset = k + sqrt(k * k - d);
//                Mat.a += koffset;
//                Mat.c += koffset;
//                A = (1. - Mat.a) * (1. - Mat.c) - Mat.b * Mat.b;
//            }

            amplifi_mat[i][j].a = (1. - Mat.c) / A;
            amplifi_mat[i][j].b = -Mat.b / A;
            amplifi_mat[i][j].c = (1. - Mat.a) / A;
            amplifi_mat[i][j].d = amplifi_mat[i][j].b;
        }
}
