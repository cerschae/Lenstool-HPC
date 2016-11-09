#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/*
JPK - OMP - 20 Mai 1995

comp the derivatives of chi2 of images using a source object

*/

extern double   **im, **tim;
extern int  nx, ny;
extern double   dlsds, xmin, xmax, ymin, ymax;

void    comp_dchi_osv(double *vect, double *dvect)
{
    double  step = .001;
    register int    i, j;
    int n = 6;
    double  chi1, chi2, vect1[6], vect2[6];

    for (j = 0; j < n; j++)
        vect2[j] = vect1[j] = vect[j];

    for (i = 0; i < n; i++)
    {
        vect1[i] -= step;
        vect2[i] += step;
        chi1 = comp_chi_osv(vect1);
        chi2 = comp_chi_osv(vect2);
        dvect[i] = (chi2 - chi1) / (vect2[i] - vect[i]);
        vect1[i] = vect2[i] = vect[i];
    }

}
