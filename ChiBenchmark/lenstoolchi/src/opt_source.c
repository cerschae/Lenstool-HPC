#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"

/*
JPK - OMP - 20 Mai 1995

optimize the source shape of a multiple image not well resolved

*/

double   **im, **tim;
int     nx, ny;
long int nbs;
double   dlsds, xmin, xmax, ymin, ymax;
struct  galaxie osource[3];


void    opt_source()
{
    extern  struct  g_mode  M;
    extern struct g_observ O;
    extern  struct  g_pixel imFrame; //PSF
    extern  struct  pot lens[];

    register    int i, k, l, m, n, o, p; //ii,j,
    //int   nimage;
    //double    err2;
    double  **dim;
    double  **imTrue;
    struct  galaxie best0, best1;
    double  bm = 100, merrit, sigim = 0.001;
    double  x0, y0, a0, b0, t0, i0;//,in0;
    double  x1, y1, a1, b1, t1, i1;//,in1;
    double  vect[12];//,dvect[6];
    FILE    *out;
    int nv, niter;
    double  ftol, fret;


    NPRINTF(stderr, "DO: optimize source of a multiple image not well resolved\n");

    /* read source file  and compute the distance ratio -----------------*/
//    f_source(M.sourfile, osource, &nbs);
    f_shape(&nbs, osource, M.sourfile, 0);
    if (nbs > 3) nbs = 3;

    for (i = 0; i < nbs; i++)
    {
        osource[i].dr = dlsds = dratio(lens[0].z, osource[i].z);
        osource[i].I0 = osource[i].mag;
    }

    NPRINTF(stderr, "nbs:%ld dlsds:: %.3lf\n", nbs, dlsds);

    /* read observed image -----------------------------------------*/
    NPRINTF(stderr, "READ: image frame\n");
    im = (double **) readimage(&imFrame);

    /* define test image with the definition of the image read ----------*/
    nx = imFrame.nx;
    ny = imFrame.ny;
    xmin = imFrame.xmin;
    xmax = imFrame.xmax;
    ymin = imFrame.ymin;
    ymax = imFrame.ymax;

    tim = (double **) alloc_square_double(nx, ny);
    dim = (double **) alloc_square_double(nx, ny);
    imTrue = (double **) alloc_square_double(5 * nx, 5 * ny);

    cp_im(tim, nx, ny, xmin, xmax, ymin, ymax, osource, nbs);
    cv_cpsf(tim, nx, ny, xmin, xmax, ymin, ymax, O.seeing);
    cp_im(imTrue, 5*nx, 5*ny, xmin, xmax, ymin, ymax, osource, nbs);
    cp_diffim(tim, im, nx, ny, dim);
    merrit = cp_errdiff(tim, im, nx, ny, sigim);
    wrf_fits("im.fits", im, nx, ny, xmin, xmax, ymin, ymax);
    wrf_fits("tim.fits", tim, nx, ny, xmin, xmax, ymin, ymax);
    wrf_fits("true.fits", imTrue, 5*nx, 5*ny, xmin, xmax, ymin, ymax);
    wrf_fits("dim.fits", dim, nx, ny, xmin, xmax, ymin, ymax);


    vect[0] = osource[0].C.x;
    vect[1] = osource[0].C.y;
    vect[2] = osource[0].E.a;
    vect[3] = osource[0].E.b;
    vect[4] = osource[0].E.theta;
    vect[5] = osource[0].I0;
    NPRINTF(stderr, "m: %.4lf %.4lf %.4lf %.4lf %.1lf %.3lf %.4lf %.3lf \n",
            vect[0], vect[1], vect[2], vect[3], vect[4]*RTD, osource[0].z, vect[5], merrit);

    x0 = osource[0].C.x;
    y0 = osource[0].C.y;
    a0 = osource[0].E.a;
    b0 = osource[0].E.b;
    t0 = osource[0].E.theta;
    i0 = osource[0].I0;
    best0 = osource[0];

    x1 = osource[1].C.x;
    y1 = osource[1].C.y;
    a1 = osource[1].E.a;
    b1 = osource[1].E.b;
    t1 = osource[1].E.theta;
    i1 = osource[1].I0;
    best1 = osource[1];

    k = l = m = n = o = p = 0;
    /*
    for(k=-1;k<2;k++)
    for(l=-1;l<2;l++)
    */
    for (m = -1; m < 2; m++)
        for (n = -1; n < 2; n++)
            for (o = -1; o < 2; o++)
                for (p = -1; p < 2; p++)
                {
                    vect[0] = x0 + k * 0.01;
                    vect[1] = y0 + l * 0.01;
                    vect[2] = a0 * (1 + m * 0.1);
                    vect[3] = b0 * (1 + n * 0.1);
                    vect[4] = t0 + p * 0.1;
                    vect[5] = i0 * (1 + o * 0.05);

                    vect[6] = x1 + k * 0.01;
                    vect[7] = y1 + l * 0.01;
                    vect[8] = a1 * (1 + m * 0.1);
                    vect[9] = b1 * (1 + n * 0.1);
                    vect[10] = t1 + p * 0.1;
                    vect[11] = i1 * (1 + o * 0.05);

                    nv = 12;
                    ftol = 0.01;
                    fret = comp_chi_osv(vect);
                    if (fret < 5.)
                    {
                        frprmn(vect, nv, ftol, &niter, &fret, comp_chi_osv, comp_dchi_osv);
                        NPRINTF(stderr, "niter:%d fret:%.3lf\n", niter, fret);
                    }

                    NPRINTF(stderr, "m: %.4lf %.4lf %.4lf %.4lf %.1lf %.3lf %.4lf %.3lf \n",
                            vect[0], vect[1], vect[2], vect[3], vect[4]*RTD, osource[0].z, vect[5], fret);

                    if (fret < bm)
                    {
                        best0.C.x = vect[0];
                        best0.C.y = vect[1];
                        best0.E.a = vect[2];
                        best0.E.b = vect[3];
                        best0.E.theta = vect[4];
                        best0.I0 = vect[5];

                        best1.C.x = vect[6];
                        best1.C.y = vect[7];
                        best1.E.a = vect[8];
                        best1.E.b = vect[9];
                        best1.E.theta = vect[10];
                        best1.I0 = vect[11];
 
                        bm = fret;
                     }
                 }
 
     out = fopen("best.dat", "w");
     fprintf(out, "#REFERENCE 3 %lf %lf\n", M.ref_ra, M.ref_dec);
     fprintf(out, "1 %.4lf %.4lf  %.4lf %.4lf %.1lf %.3lf %.2lf %.3lf\n",
            best0.C.x, best0.C.y,
            best0.E.a, best0.E.b, best0.E.theta*RTD, osource[0].z, best0.I0, bm);

    NPRINTF(stderr, "1 %.4lf %.4lf  %.4lf %.4lf %.1lf %.3lf %.2lf %.3lf\n",
            best0.C.x, best0.C.y,
            best0.E.a, best0.E.b, best0.E.theta*RTD, osource[0].z, best0.I0, bm);

     fprintf(out, "2 %.4lf %.4lf  %.4lf %.4lf %.1lf %.3lf %.2lf %.3lf\n",
            best1.C.x, best1.C.y,
            best1.E.a, best1.E.b, best1.E.theta*RTD, osource[1].z, best1.I0, bm);

    NPRINTF(stderr, "2 %.4lf %.4lf  %.4lf %.4lf %.1lf %.3lf %.2lf %.3lf\n",
            best1.C.x, best1.C.y,
            best1.E.a, best1.E.b, best1.E.theta*RTD, osource[1].z, best1.I0, bm);

    fclose(out);

    osource[0] = best0;
    osource[1] = best1;

    cp_im(tim, nx, ny, xmin, xmax, ymin, ymax, osource, nbs);
    cv_cpsf(tim, nx, ny, xmin, xmax, ymin, ymax, O.seeing);

    cp_im(imTrue, 5*nx, 5*ny, xmin, xmax, ymin, ymax, osource, nbs);

    /* compute  differnce of test and observed image  --------------*/
    cp_diffim(tim, im, nx, ny, dim);

    if ( M.iref > 0 )
    {
        wrf_fits_abs("im.fits", im, nx, ny, xmin, xmax, ymin, ymax, M.ref_ra, M.ref_dec);
        wrf_fits_abs("tim.fits", tim, nx, ny, xmin, xmax, ymin, ymax, M.ref_ra, M.ref_dec);
        wrf_fits_abs("true.fits", imTrue, 5*nx, 5*ny, xmin, xmax, ymin, ymax, M.ref_ra, M.ref_dec);
        wrf_fits_abs("dim.fits", dim, nx, ny, xmin, xmax, ymin, ymax, M.ref_ra, M.ref_dec);
    }
    else
    {
        wrf_fits("im.fits", im, nx, ny, xmin, xmax, ymin, ymax);
        wrf_fits("tim.fits", tim, nx, ny, xmin, xmax, ymin, ymax);
        wrf_fits("true.fits", imTrue, 5*nx, 5*ny, xmin, xmax, ymin, ymax);
        wrf_fits("dim.fits", dim, nx, ny, xmin, xmax, ymin, ymax);
    }

    free_square_double(im, nx);
    free_square_double(tim, nx);
    free_square_double(dim, nx);
    free_square_double(imTrue, 5*nx);


}

