#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
/* Append a Poisson noise to the pixels of image **z
 */
void d_bruiter(double **z, int nx, int ny)
{
    extern struct g_observ O;

    int i, j;
    double res;


    for (i = 0; i < ny; i++)
        for (j = 0; j < nx; j++)
        {
            res = z[i][j];
            res = (res + O.SKY) * O.gain;
            res = floor(d_poisson(res, &O.idum) / O.gain);
            z[i][j] = res;
        }
}
