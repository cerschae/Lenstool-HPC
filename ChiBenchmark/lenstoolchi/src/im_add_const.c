#include<stdio.h>
#include<math.h>
#include "fonction.h"
#include "constant.h"
#include"dimension.h"
#include "structure.h"
#include "lt.h"

long int im_add_const(double **im, int nx, int ny, double c)
{
    register int    i, j;
    double **im_res;

    im_res = (double **)alloc_square_double(nx, ny);

    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            im_res[i][j] = im[i][j] + c;

    return((long int ) im_res);
}
