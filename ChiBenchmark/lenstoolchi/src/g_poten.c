#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*      nom:        g_poten             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    g_poten(int ipoten, int np, double z, char *file)
{
    const extern  struct  g_mode          M;
    const extern  struct  g_frame     F;
    const extern  struct  g_cosmo     C;
    const extern  struct  pot lens[];

    register int    i, j;
    double  dx, dy;
    double  dcrit, dlsds, dl;
    struct  point   pi;
    double  **poten;
    //int   size[4];

    if (ipoten == 1)
    {
        NPRINTF(stderr,
                "COMP: projected relative potential map for z_s=%.3lf =>%s\n", z, file);
    }
    else
    {
        NPRINTF(stderr, "COMP: projected absolute potential map =>%s\n", file);
    }

    dlsds = dratio(lens[0].z, z);
    dl = distcosmo1(lens[0].z);
    dcrit = cH2piG * C.h / dl / dlsds;
    dx = (F.xmax - F.xmin) / (np - 1);
    dy = (F.ymax - F.ymin) / (np - 1);

    poten = (double **) alloc_square_double(np, np);

    for (j = 0; j < np; j++)
    {
        pi.y = j * dy + F.ymin;
        for (i = 0; i < np; i++)
        {
            pi.x = i * dx + F.xmin;
            poten[j][i] = e_pot(pi, dlsds);
            if (ipoten == 2)
                poten[j][i] *= dcrit;
            else if (ipoten == 3)
                poten[j][i] /= dlsds;
        };
    };

    wrf_fits(file, poten, np, np, F.xmin, F.xmax, F.ymin, F.ymax);

    free_square_double(poten, np);
}
