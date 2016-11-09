#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"

void o_pixel(double **zz, int nx, int ny, double scale, double xmin, 
        double xmax, double ymin, double ymax, struct galaxie *source, double **dpl_x, double **dpl_y);

/****************************************************************/
/*      nom:        e-pixel             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Create a FITS file containing the computed images of the arcs
 * and arclets in the source plane
 */
void    e_pixel(int np, char *iname, char *sname, struct galaxie *source)
{
    const extern struct g_frame   F;
    const extern struct g_mode    M;
    const extern struct g_observ  O;
    //const extern    struct  g_large L;

    double xmin, ymin, xmax, ymax, dx, dy;
    double f, scale;
    int nx, ny, ny_sav;
    double **zz;

    xmin = F.xmin;
    xmax = F.xmax;
    ymin = F.ymin;
    ymax = F.ymax;

    dx = xmax - xmin;
    dy = ymax - ymin;
    nx = np;    // default: assign np to nx and adjust ny
    scale = dx / (nx-1);
    ny = dy / scale + 1;  // warning: trunc(ny) behavior 
//    f = dx / scale;
//    nx = (int) f;  // BUG in C
//    f = dy / scale;
//    ny = (int) f;
    if (O.setbin)
    {
        if ( O.bin * nx <= NTMAX && O.bin * ny <= NTMAX )
        {
            nx = O.bin * nx;
            ny = O.bin * ny;
            scale = scale / ((double) O.bin);
        }
        else
        {
            fprintf(stderr, "ERROR: reaching maximal size of array\n");
            exit(-1);
        }
    }

    dx = (nx - 1) * scale;
    dy = (ny - 1) * scale;
    xmax = xmin + dx;
    ymax = ymin + dy;

    NPRINTF(stderr, "\timage (%d,%d) s=%.3lf [%.3lf,%.3lf] [%.3lf,%.3lf]\n",
            nx, ny, scale, xmin, xmax, ymin, ymax);

    zz = (double **) alloc_square_double(ny, nx);

    o_pixel(zz, nx, ny, scale, xmin, xmax, ymin, ymax, source, NULL, NULL);


    if (O.setseeing)
        d_seeing(zz, nx, ny, scale);

    ny_sav = ny;
    if (O.setbin)
        d_binning(zz, &nx, &ny, O.bin);

    if (O.bruit)
        d_bruiter(zz, nx, ny);

    NPRINTF(stderr, "COMP: image of arcs and arclets --> %s\n", iname);

    if ( M.iref > 0 )
        wrf_fits_abs(iname, zz, nx, ny, xmin, xmax, ymin, ymax, M.ref_ra, M.ref_dec);
    else
        wrf_fits(iname, zz, nx, ny, xmin, xmax, ymin, ymax);

    free_square_double(zz, ny_sav);

    // generate a simulated image of the source
    const extern struct g_source  S;
    extern struct g_pixel ps;

    NPRINTF(stdout, "Simulated %ld sources\n", S.ns);
    // define a source box
    double dlsds = source[0].dl0s / source[0].dos;

    if( strcmp(M.centerfile, "") )
    {
        ps.pixelx = scale / ps.ech;
        ps.pixely = scale / ps.ech;
        s_sourcebox(&ps, M.centerfile, dlsds);
    }
    else
    {
        dx = ps.xmax - ps.xmin;
        dy = ps.ymax - ps.ymin;
        ps.pixelx = dx / (ps.nx-1);
        ps.pixely = dy / (ps.ny-1);
    }


    zz = alloc_square_double(ps.ny, ps.nx);
    struct point A; 
    int j, k, l;
    for ( j = 0; j < ps.ny; j++ )
    {
        A.y = ps.ymin + ((double)j + 1) * ps.pixely;
        for ( k = 0; k < ps.nx; k++ )
        {
            A.x = ps.xmin + ((double)k + 1) * ps.pixelx;
            zz[j][k] = 0.;
            for( l = 0; l < S.ns; l++ )
                zz[j][k] += d_profil(A.x, A.y, &source[l]);
        }
    }

    if ( M.iref > 0 )
        wrf_fits_abs(sname, zz, ps.nx, ps.ny, ps.xmin, ps.xmax, ps.ymin, ps.ymax, M.ref_ra, M.ref_dec);
    else
        wrf_fits(sname, zz, ps.nx, ps.ny, ps.xmin, ps.xmax, ps.ymin, ps.ymax);

    free_square_double(zz, ps.ny);
}

/* Function to simulate an image from a source
 * Called by o_chi() and e_pixel()
 * zz is filled in a field defined by champ section with xmin, ymin, xmax, ymax
 *
 */ 
void o_pixel(double **zz, int nx, int ny, double scale, double xmin, 
        double xmax, double ymin, double ymax, struct galaxie *source, double **dpl_x, double **dpl_y)
{

    const extern struct g_mode    M;
    const extern struct g_source  S;
    


    int j, k, l;
    struct point pi, ps;

    if (dpl_x == NULL)
    {
        dpl_x = (double **) alloc_square_double(ny, nx);
        dpl_y = (double **) alloc_square_double(ny, nx);
        for (j = 0; j < ny; j++)
        {
            pi.y = ymin + j * scale;
            for (k = 0; k < nx; k++)
            {
                pi.x = xmin + k * scale;
                e_dpl(&pi, 1., &ps);
                dpl_x[j][k] = pi.x - ps.x;
                dpl_y[j][k] = pi.y - ps.y;
            }
        }
    }
   
#pragma omp parallel for private(ps,pi,j,k,l)
    for (j = 0; j < ny; j++)
    {
        pi.y = ymin + j * scale;
        for (k = 0; k < nx; k++)
        {
            pi.x = xmin + k * scale;
            zz[j][k] = 0.;
            for (l = 0 ; l < S.ns ; l++ )
            {
                ps.x = pi.x - dpl_x[j][k] * source[l].dr;
                ps.y = pi.y - dpl_y[j][k] * source[l].dr;
                double f = d_profil(ps.x, ps.y, &source[l]);
                /*                if ((f>sqrt(O.SKY/O.gain)) || (f>source[l].I0/100.))
                                    zz[j][k]+=d_integrer(source[l],pi.x,pi.y,scale/2.,f,0);
                        else
                */
                zz[j][k] += f;
            }
        }
    }

}
