#include<stdio.h>
#include<math.h>
#include "fonction.h"
#include "constant.h"
#include"dimension.h"
#include "structure.h"

void im_stat(double **im, int nx, int ny, double *mean, double *disp)
{
    register int    i, j, k;
    double  s, s2, moy, sig;


    s = s2 = 0.;
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
        {
            s += im[i][j];
            s2 += im[i][j] * im[i][j];
        };
    moy = s / nx / ny;
    sig = sqrt(s2 / nx / ny - moy * moy);

    k = 0;
    s = s2 = 0.;
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            if (fabs(im[i][j] - moy) < 3.*sig)
            {
                s += im[i][j];
                s2 += im[i][j] * im[i][j];
                k++;
            };
    moy = s / nx / ny;
    sig = sqrt(s2 / nx / ny - moy * moy);

    k = 0;
    s = s2 = 0.;
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            if (fabs(im[i][j] - moy) < 3.*sig)
            {
                s += im[i][j];
                s2 += im[i][j] * im[i][j];
                k++;
            };
    *mean = s / k;
    *disp = sqrt(s2 / k - (*mean) * (*mean));
}
