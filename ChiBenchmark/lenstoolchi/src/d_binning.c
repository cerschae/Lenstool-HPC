#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"

void d_binning(double **im, int *nx, int *ny, int bin)
{
    const extern struct g_mode   M;
    register int i, j, ii, jj;
    int    nxb, nyb;
    double s, bb;

    bb = bin * bin;

    NPRINTF(stderr, "Bin image by a factor %d\n", bin);

    nxb = (int) ((double) * nx) / ((double)bin);
    nyb = (int) ((double) * ny) / ((double)bin);

    for (i = 0; i < *ny; i += bin)
        for (j = 0; j < *nx; j += bin)
        {
            s = 0.;
            for (ii = 0; ii < bin; ii++)  // y
                for (jj = 0; jj < bin; jj++)  //x
                    s += im[i+ii][j+jj];

            im[i/bin][j/bin] = s / bb;
        }

    *nx = nxb;
    *ny = nyb;
}
