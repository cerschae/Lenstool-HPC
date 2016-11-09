#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*      nom:        g_curv              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    g_curv(int icurv, int np, double z, char *file1, char *file2, char *file3)
{
    const extern  struct  g_mode          M;
    const extern  struct  g_frame     F;
    const extern  struct  pot lens[];

    register int    i, j;
    double  dx, dy;
    double  dl0s, dos, dlsds;
    struct  point   pi;
    struct  matrix  gg;
    double  **cxx, **cxy, **cyy;

    NPRINTF(stderr, "COMP: curvature map for zs=%.3lf =>%s %s %s\n", z, file1, file2, file3);

    dl0s = distcosmo2( lens[0].z, z);
    dos = distcosmo1( z );
    dlsds = dl0s / dos;
    dx = (F.xmax - F.xmin) / (np - 1);
    dy = (F.ymax - F.ymin) / (np - 1);

    cxx = (double **) alloc_square_double(np, np);
    cxy = (double **) alloc_square_double(np, np);
    cyy = (double **) alloc_square_double(np, np);

    for (j = 0; j < np; j++)
    {
        pi.y = j * dy + F.ymin;
        for (i = 0; i < np; i++)
        {
            pi.x = i * dx + F.xmin;
            gg = e_grad2(&pi, dl0s, z);
            gg.a /= dos;
            gg.b /= dos;
            gg.c /= dos;

            cxx[j][i] = gg.a;
            cxy[j][i] = gg.b;
            cyy[j][i] = gg.c;
        };
    };

    if (M.iref > 0)
    {
        wrf_fits_abs(file1, cxx, np, np, F.xmin, F.xmax, F.ymin, F.ymax, M.ref_ra, M.ref_dec);
        wrf_fits_abs(file2, cxy, np, np, F.xmin, F.xmax, F.ymin, F.ymax, M.ref_ra, M.ref_dec);
        wrf_fits_abs(file3, cyy, np, np, F.xmin, F.xmax, F.ymin, F.ymax, M.ref_ra, M.ref_dec);
    }
    else
    {
        wrf_fits(file1, cxx, np, np, F.xmin, F.xmax, F.ymin, F.ymax);
        wrf_fits(file2, cxy, np, np, F.xmin, F.xmax, F.ymin, F.ymax);
        wrf_fits(file3, cyy, np, np, F.xmin, F.xmax, F.ymin, F.ymax);
    }

    free_square_double(cxx, np);
    free_square_double(cxy, np);
    free_square_double(cyy, np);
}
