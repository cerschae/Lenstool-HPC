#include<stdio.h>
#include<math.h>
#include "fonction.h"
#include "constant.h"
#include"dimension.h"
#include "structure.h"
#include "lt.h"

long int im_bin(double **im, int nx, int ny, int bin)
{
    register int    i, j, ii, jj;
    double s, B;
    double **imbin;

    B = bin;

    imbin = (double **)alloc_square_double(nx / bin, ny / bin);

    for (i = 0; i < nx; i += bin)
        for (j = 0; j < ny; j += bin)
        {
            s = 0.;
            for (ii = 0; ii < bin; ii++)
                for (jj = 0; jj < bin; jj++)
                    s += im[i+ii][j+jj];

            imbin[i/bin][j/bin] = s / B / B;
        };

    return((long int ) imbin);
}
