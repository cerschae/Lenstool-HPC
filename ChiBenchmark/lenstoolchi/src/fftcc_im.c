#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

void    fftcc_im(double **r_im, double **i_im, double **tfr_im, double **tfi_im, int n, int flag)
{
    const extern  struct  g_mode          M;
    register int    i, j, k;
    int err, isign, nbr_lin, nbr_col;
    int ndim, nn[2];

    double  *image_FT;

    nbr_lin = n;
    nbr_col = n;

    if (flag == -1)
    {
        err = (int)split_image (r_im, nbr_lin, nbr_col);
        if (err != 0)
        {
            fprintf(stderr, "\nFATAL ERROR: when splitting input image\n\n");
            exit(-1);
        }
        err = (int)split_image (i_im, nbr_lin, nbr_col);
        if (err != 0)
        {
            NPRINTF(stderr, "\nFATAL ERROR: when splitting input image\n\n");
            exit(-1);
        }
    }

    /*    FAST FOURIER TRANSFORM OF THE IMAGE */

    image_FT = (double *) malloc (2 * nbr_lin * nbr_col * sizeof(double));
    for (i = 0; i < nbr_lin; i++)
        for (j = 0; j < nbr_col; j++)
        {
            image_FT[2*(j+i*nbr_col)] = r_im[i][j];
            image_FT[2*(j+i*nbr_col)+1] = i_im[i][j];
        }
    ndim = 2;
    isign = -flag;
    nn[0] = nbr_lin;
    nn[1] = nbr_col;

    fft(image_FT - 1, nn - 1, ndim, isign);

    for (k = 0; k < nbr_lin*nbr_col; k++)
    {
        i = k / nbr_col;
        j = k - i * nbr_col;
        if (flag == 1)
        {
            tfr_im[i][j] = image_FT[2*k];
            tfi_im[i][j] = image_FT[2*k+1];
        }
        else if (flag == -1)
        {
            tfr_im[i][j] = image_FT[2*k] / ((double)(nbr_lin * nbr_col));
            tfi_im[i][j] = image_FT[2*k+1] / ((double)(nbr_lin * nbr_col));
        }
    }
    free(image_FT);


    /*    SPLITS RESULTS TO HAVE ORIGIN AT CENTER OF IMAGE
    */
    if (flag == 1)
    {
        err = (int)split_image(tfr_im, nbr_lin, nbr_col);
        if (err != 0)
        {
            NPRINTF(stderr, "\a\nFATAL ERROR: problem to split results\n\n");
            exit(-1);
        }
        err = (int)split_image(tfi_im, nbr_lin, nbr_col);
        if (err != 0)
        {
            NPRINTF(stderr, "\a\nFATAL ERROR: problem to split results\n\n");
            exit(-1);
        }
    }


}
