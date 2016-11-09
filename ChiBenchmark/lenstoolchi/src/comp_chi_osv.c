#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/*
JPK - OMP - 20 Mai 1995

comp the chi2 of images using a source object

*/
extern double   **im, **tim;
extern int     nx, ny, nbs;
extern double   dlsds, xmin, xmax, ymin, ymax;
extern  struct  galaxie osource[3];

double  comp_chi_osv(double *vect)
{
    extern struct g_observ O;

    register    int i, j;
    double  sigim = 0.001;
    double  merrit;

    osource[0].C.x = vect[0];
    osource[0].C.y = vect[1];
    osource[0].E.a = vect[2];
    osource[0].E.b = vect[3];
    osource[0].E.theta = vect[4];
    osource[0].I0 = vect[5];
    osource[0].dr = dlsds;


    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            tim[i][j] = 0.;
    cp_im(tim, nx, ny, xmin, xmax, ymin, ymax, osource, nbs);
    cv_cpsf(tim, nx, ny, xmin, xmax, ymin, ymax, O.seeing);

    /* compute  merrit function of the difference of the test and observed images*/
    merrit = cp_errdiff(tim, im, nx, ny, sigim);    /* sigim is the image noise */

    return(merrit);
}
