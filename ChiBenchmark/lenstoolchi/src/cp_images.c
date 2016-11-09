#include<stdio.h>
#include<math.h>

/* function compute  difference of two images    ----------- */

void cp_diffim(double **im1, double **im2, int ni, int nj, double **resim)
{
    register int    i, j;

    for (i = 0; i < ni; i++)
        for (j = 0; j < nj; j++)
            resim[i][j] = (im1[i][j] - im2[i][j]);
}


/* function compute ratio of two images    ----------- */

void cp_divim(double **im1, double **im2,
              int ni, int nj,
              double **resim, double dbz)
{
    register int    i, j;

    for (i = 0; i < ni; i++)
        for (j = 0; j < nj; j++)
            if (im2[i][j] != 0.)
                resim[i][j] = (im1[i][j] / im2[i][j]);
            else
                resim[i][j] = dbz;
}


/* function compute  exponetiel of an image    ----------- */

void cp_expim(double **im, int ni, int nj, double **resim)
{
    register int    i, j;

    for (i = 0; i < ni; i++)
        for (j = 0; j < nj; j++)
            resim[i][j] = exp(im[i][j]);

}

/* function compute  product of two images    ----------- */

void cp_mulim(double **im1, double **im2,
              int ni, int nj, double **resim)
{
    register int    i, j;

    for (i = 0; i < ni; i++)
        for (j = 0; j < nj; j++)
            resim[i][j] = (im1[i][j] * im2[i][j]);
}


/* function compute  the square of an image   ----------- */

void cp_sqrim(double **im, int ni, int nj, double **resim)
{
    register int    i, j;

    for (i = 0; i < ni; i++)
        for (j = 0; j < nj; j++)
            resim[i][j] = im[i][j] * im[i][j];
}


/* function compute  the sum of two images    ----------- */

void cp_sumim(double **im1, double **im2,
              int ni, int nj,
              double **resim)
{
    register int    i, j;

    for (i = 0; i < ni; i++)
        for (j = 0; j < nj; j++)
            resim[i][j] = (im1[i][j] + im2[i][j]);
}

/* function threshold an image    ----------- */

void cp_thim(double **im, int ni, int nj,
             double min, double max,
             double **resim)
{
    register int    i, j;

    for (i = 0; i < ni; i++)
        for (j = 0; j < nj; j++)
        {
            if (im[i][j] > max)
                resim[i][j] = max;
            else if (im[i][j] < min)
                resim[i][j] = min;
            else
                resim[i][j] = max;
        }

}


