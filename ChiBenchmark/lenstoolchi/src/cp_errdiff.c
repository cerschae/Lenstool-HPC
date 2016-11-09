#include<stdio.h>

/* function compute  norme of (test image-image)    ----------- */

double cp_errdiff(double **im1, double **im2, int ni, int nj, double sigim)
{
    register int    i, j;
    double ndiff = 0;

    for (i = 0; i < ni; i++)
        for (j = 0; j < nj; j++)
            ndiff += (im1[i][j] - im2[i][j]) * (im1[i][j] - im2[i][j]);

    ndiff /= ni * nj;
    ndiff /= sigim * sigim;

    return(ndiff);
}


