#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*      nom:        g_time              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    g_time(int flag, int np, double z, char *file)
{
    const extern  struct  g_mode          M;
    const extern  struct  g_frame     F;
    const extern  struct  g_cosmo     C;
    const extern  struct  pot lens[];

    register int    i, j;
    double  cst, dlsds, dl;
    struct  point   pi, gr;
    double  **time;
    //int   size[4];


    dlsds = dratio(lens[0].z, z);
    dl = distcosmo1(lens[0].z);

    if (flag == 1)
    {
        cst = (1 + lens[0].z) * dl / dlsds * th_a2_day / C.h;
        NPRINTF(stderr, "COMP: time-delay (Image Plane) for z_s=%.3lf in days =>%s\n",
                z, file);
    }
    else
    {
        cst = (1 + lens[0].z) * dl / dlsds * tyaa / C.h;
        NPRINTF(stderr, "COMP: time-delay (Image Plane) for z_s=%.3lf in years =>%s\n",
                z, file);
    }

    time = (double **) alloc_square_double(np, np);

    for (j = 0; j < np; j++)
    {
        pi.y = j * (F.ymax - F.ymin) / (np - 1) + F.ymin;
        for (i = 0; i < np; i++)
        {
            pi.x = i * (F.xmax - F.xmin) / (np - 1) + F.xmin;
            gr = e_grad(&pi);
            gr.x *= dlsds;
            gr.y *= dlsds;


            /* time is calculated in year */
            time[j][i] = cst * ((gr.x * gr.x + gr.y * gr.y) / 2. - e_pot(pi, dlsds));

        }
    }

    if (M.iref > 0)
    {
        wrf_fits_abs(file, time, np, np, F.xmin, F.xmax, F.ymin, F.ymax, M.ref_ra, M.ref_dec);
    }
    else
    {
        wrf_fits(file, time, np, np, F.xmin, F.xmax, F.ymin, F.ymax);
    }

    free_square_double(time, np);
}
