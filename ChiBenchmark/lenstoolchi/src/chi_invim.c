#include <fonction.h>

/*
*       nom:        chi_invim
*       auteur:     Jean-Paul Kneib
*       date:       10/03/94
*       place:      ESO
*
* major changes:
*       date: oct 94
*       place: IoA Cambridge
*
* Return the error associated to the transformation image -> source
*
* Global variables used :
* - imFrame, ps, M
* - in err_invim() : iline, radial, tangent, nrline, ntline
* - in criticinv() : CL, radial, tangent, lens, nrline, ntline, flagr, flagt, G, lens_table
* - in s_sourcebox() : M, G, lens, lens_table
* - in do_itos() : M, imFrame, ps, G, lens, lens_table, iline, radial, tangent, nrline, ntline
*/

double chi_invim(double **im, struct pixlist pl[], int npl, double dlsds,
                 double **so, double **er, int **imu)
{
    extern struct g_pixel ps;
    const extern struct g_pixel imFrame;
    const extern struct g_mode  M;

    double err;

    /*
    * compute the new critical and caustics line - quick version
    */
    criticinv(imFrame.xmin, imFrame.xmax, imFrame.ymin, imFrame.ymax);

    /*
    * define the source box as the barycenter of the source image-center points
    */
    s_sourcebox(&ps, M.centerfile, dlsds);

    /*
    *  do the transform image -> source + compute the error
    */
    do_itos(im, pl, npl, dlsds, so, er, imu);

    /*
    * estimate the error
    */
    err = err_invim(er, imu);

    return(err);
}
