#include<stdio.h>
#include<math.h>
#include "fonction.h"
#include "constant.h"
#include"dimension.h"
#include "structure.h"
#include "lt.h"

long int im_convolve(double **ima, int nx, int ny, double **filt, int nf)
{
    register int    i, j, ii, jj;
    int i_im, j_im, nfd;
    double s;
    double **im_conv;

    nfd = (nf + 1) / 2;
    nfd = (nf - 1) / 2;

    im_conv = (double **)alloc_square_double(nx, ny);

    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
        {
            s = 0.;
            for (ii = 0; ii < nf; ii++)
            {
                i_im = i - nfd + ii;
                for (jj = 0; jj < nf; jj++)
                {
                    j_im = j - nfd + jj;
                    if ((i_im >= 0) && (i_im < nx) && (j_im >= 0) && (j_im < ny))
                        s += ima[i_im][j_im] * filt[ii][jj];
                };
            };

            im_conv[i][j] = s;
        };

    return( (long int) im_conv);
}
