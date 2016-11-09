#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"fonction.h"
#include "lt.h"
/*
JPK - OMP - 20 Mai 97

Convolve an image by a circular gaussien PSF of FWHM=seeing

calls crea_filtre(seeing,pixel-size,filter,n)

*/

void    cv_cpsf(double **im,
                int nx, int ny,
                double xmin, double xmax, double ymin, double ymax,
                double seeing)
{
    register int    i, j;
    int n, N, err;
    double  scalex, scaley;

    double **z, **tfr_z, **tfi_z, **filter, **tfr_f, **tfi_f;
    double **r_conv, **i_conv, **tfr_conv, **tfi_conv;

    /* compute pixel scale ---------------------------------------*/
    scalex = (xmax - xmin) / (nx - 1);
    scaley = (ymax - ymin) / (ny - 1);

    if ( fabs(scalex - scaley) > 1e-5)
    {
        fprintf(stderr, "WARNING: pixel is not a square\n");
    }

    /* include the image in a larger 2^N image ---------------------*/

    N = (int) ceil(log(((double)Max(nx, ny))) / log(2.));
    n = (int) pow(2., ((double)N));
    z = (double **) alloc_square_double(n, n);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            if ((i < nx) && (j < ny))
                z[i][j] = im[i][j];
            else
                z[i][j] = 0.;
        }

    /* compute  filter  -------------------------------------------*/

    filter = (double **) alloc_square_double(n, n);
    crea_filtre(seeing, (scalex + scaley) / 2., filter, n);

    /* Compute Fourrier Transform  -------------------------------*/
    tfr_z = (double **) alloc_square_double(n, n);
    tfi_z = (double **) alloc_square_double(n, n);
    fftc_im(z, tfr_z, tfi_z, n, 1);
    free_square_double(z, n);

    tfr_f = (double **) alloc_square_double(n, n);
    tfi_f = (double **) alloc_square_double(n, n);
    fftc_im(filter, tfr_f, tfi_f, n, 1);
    free_square_double(filter, n);

    tfr_conv = (double **) alloc_square_double(n, n);
    tfi_conv = (double **) alloc_square_double(n, n);
    ic_product(tfr_z, tfi_z, tfr_f, tfi_f, tfr_conv, tfi_conv, n, n);
    free_square_double(tfr_z, n);
    free_square_double(tfi_z, n);
    free_square_double(tfr_f, n);
    free_square_double(tfi_f, n);

    r_conv = (double **) alloc_square_double(n, n);
    i_conv = (double **) alloc_square_double(n, n);
    fftcc_im(tfr_conv, tfi_conv, r_conv, i_conv, n, -1);
    free_square_double(tfr_conv, n);
    free_square_double(tfi_conv, n);

    err = (int) split_image(r_conv, n, n);
    if (err != 0)
    {
        fprintf(stderr, "\nFATAL ERROR: when splitting input image\n");
        exit(-1);
    }

    /* copy results into the original image --------------------*/

    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            im[i][j] = r_conv[i][j];

    free_square_double(r_conv, n);
    free_square_double(i_conv, n);

}
