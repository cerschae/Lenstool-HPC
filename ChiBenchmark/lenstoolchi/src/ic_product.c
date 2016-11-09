#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

void    ic_product( double **r_im, double **i_im, double **r_im2, double **i_im2,
                    double **r_prod, double **i_prod,
                    int nbr_lin, int nbr_col)
{
    register int    i, j;
    double  prodr, prodi;

    for (i = 0; i < nbr_lin; i++)
        for (j = 0; j < nbr_col; j++)
        {
            prodr = r_im[i][j] * r_im2[i][j] - i_im[i][j] * i_im2[i][j];
            prodi = r_im[i][j] * i_im2[i][j] + i_im[i][j] * r_im2[i][j];
            r_prod[i][j] = prodr;
            i_prod[i][j] = prodi;
        }
}
