#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"

/* Append seeing effect to the image im.
 * Parameters :
 * im : image in [ny][nx] format
 * scale : pixel size in arcsec
 */
void    d_seeing(double **im, int nx, int ny, double scale)
{
    extern struct g_mode    M;
    const extern struct g_observ  O;

    register int i, j;
    int    n, N, err;

    double **z, **tfr_z, **tfi_z, **tfr_f, **tfi_f; //f
    double **tfr_conv, **tfi_conv; //**r_conv,**i_conv,

    NPRINTF(stderr, "ADD: seeing\n");

    /* on va inclure l'image dans une image plus grande de dimension n=2^N */

    N = (int) ceil(log(((double)Max(nx, ny))) / log(2.));
    n = (int) pow(2., ((double)N));

    NPRINTF(stderr, "\tResize image in 2^N for FFT [%d,%d] --> [%d,%d]\n", nx, ny, n, n);

    z = (double **) alloc_square_double(n, n);  // z[ny][nx]

    NPRINTF(stderr, "\tFill [%d,%d] image with original image\n", n, n);

    for (i = 0; i < n; i++) // ny
        for (j = 0; j < n; j++) // nx
        {
            if ((i < ny) && (j < nx))
                z[i][j] = im[i][j];
            else
                z[i][j] = 0.;
        }

    NPRINTF(stderr, "\tAllocate images for Real & Imaginary part of the FT\n");

    tfr_z = (double **) alloc_square_double(n, n);
    tfi_z = (double **) alloc_square_double(n, n);

    NPRINTF(stderr, "\tPerform the image FFT\n");

    fftc_im(z, tfr_z, tfi_z, n, 1); // z is not used anymore for image

    // By default O.filtre is null (see set_default.c) so do it!!!
    if ( !O.filtre )
    {
        NPRINTF(stderr, "\tCreate the seeing filter\n");
        printf("\tCreate the seeing filter\n");
        if(O.setseeing==1) crea_filtre(O.seeing, scale, z, n);  // z contains the seeing filter
        if(O.setseeing==2) crea_filtre_e(O.seeing_a, O.seeing_b, O.seeing_angle,scale, z, n);  // z contains the seeing filter
    }
   
   //additional normalization of seeing filter
   double psf_sum = 0;
   for (i = 0 ; i < n ; i++)
     for (j = 0 ; j < n ; j++)
       psf_sum += z[i][j];
   for (i = 0 ; i < n ; i++)
     for (j = 0 ; j < n ; j++)
       z[i][j] /= psf_sum;
   
   
    tfr_f = (double **) alloc_square_double(n, n);
    tfi_f = (double **) alloc_square_double(n, n);

    NPRINTF(stderr, "\tPerform the filter FFT\n");

    fftc_im(z, tfr_f, tfi_f, n, 1); // z is not used anymore for filter

    free_square_double(z,n);
    tfr_conv = (double **) alloc_square_double(n, n);
    tfi_conv = (double **) alloc_square_double(n, n);

    NPRINTF(stderr, "\tCompute the (image X filter) product in Fourier space\n");

    ic_product(tfr_z, tfi_z, tfr_f, tfi_f, tfr_conv, tfi_conv, n, n);
//  free_square_double(tfr_z,n);
//  free_square_double(tfi_z,n);
    free_square_double(tfr_f, n);
    free_square_double(tfi_f, n);

    //r_conv= (double **) alloc_square_double(n,n);
    //i_conv= (double **) alloc_square_double(n,n);

    NPRINTF(stderr, "\tCompute the IFFT of the product\n");

    fftcc_im(tfr_conv, tfi_conv, tfr_z, tfi_z, n, -1);
    free_square_double(tfr_conv, n);
    free_square_double(tfi_conv, n);

    NPRINTF(stderr, "\tSplit the image product of the convolution\n");

    err = (int) split_image(tfr_z, n, n);
    if (err != 0)
    {
        fprintf(stderr, "ERROR: when splitting the image\n");
        exit(-1);
    }

    for (i = 0; i < ny; i++)
        for (j = 0; j < nx; j++)
            im[i][j] = tfr_z[i][j]; // [ny][nx]

    free_square_double(tfr_z, n);
    free_square_double(tfi_z, n);

//  free_square_double(r_conv,n);
//  free_square_double(i_conv,n);

}
