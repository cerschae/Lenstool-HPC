#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*      nom:        g_dpl               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    g_dpl(int idpl, int np, double z, char *filex, char *filey)
{
    const extern  struct  g_mode          M;
    const extern  struct  g_frame     F;
    const extern  struct  pot lens[];

    register int    i, j;
    double  dx, dy;
    double  dlsds;
    struct  point   pi;
    struct  point   dpl;
    double  **dplx, **dply;

    NPRINTF(stderr, "COMP: displacment map for z_s=%.3lf =>%s %s\n", z, filex, filey);

    dlsds = dratio(lens[0].z, z);
    dx = (F.xmax - F.xmin) / (np - 1);
    dy = (F.ymax - F.ymin) / (np - 1);

    dplx = (double **) alloc_square_double(np, np);
    dply = (double **) alloc_square_double(np, np);

    for (j = 0; j < np; j++)
    {
        pi.y = j * dy + F.ymin;
        for (i = 0; i < np; i++)
        {
            pi.x = i * dx + F.xmin;
            dpl = e_grad(&pi);
            dpl.x *= dlsds;
            dpl.y *= dlsds;

            dplx[j][i] = dpl.x;
            dply[j][i] = dpl.y;
        };
    };

    if (M.iref > 0)
    {
        wrf_fits_abs(filex, dplx, np, np, F.xmin, F.xmax, F.ymin, F.ymax, M.ref_ra, M.ref_dec);
        wrf_fits_abs(filey, dply, np, np, F.xmin, F.xmax, F.ymin, F.ymax, M.ref_ra, M.ref_dec);
    }
    else
    {
        wrf_fits(filex, dplx, np, np, F.xmin, F.xmax, F.ymin, F.ymax);
        wrf_fits(filey, dply, np, np, F.xmin, F.xmax, F.ymin, F.ymax);
    }

    free_square_double(dplx, np);
    free_square_double(dply, np);
}
