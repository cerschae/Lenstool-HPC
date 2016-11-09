#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        g_ampf          */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    g_amplif(int iampf, double z, char *file)
{
    const extern  struct  g_mode          M;
    const extern  struct  g_frame     F;
//const extern    struct  g_cosmo     C;
    const extern  struct  pot lens[];

    register int    i, j;
    int nampf = 25;
    double  e, amp;
    double  dl0s, dos, dlsds;
    struct  point   pi;
    struct  pointell    ampf[25][25];
    FILE    *OUT;

    NPRINTF(stderr, "COMP: Amp. in the Image Plane for z_s=%.3lf =>%s\n", z, file);

    dl0s = distcosmo2(lens[0].z, z);
    dos = distcosmo1(z);
    dlsds = dl0s / dos;


    for (j = 0; j < nampf; j++)
    {
        pi.y = j * (F.ymax - F.ymin) / (nampf - 1) + F.ymin;
        for (i = 0; i < nampf; i++)
        {
            pi.x = i * (F.xmax - F.xmin) / (nampf - 1) + F.xmin;
            ampf[i][j].E = e_unmag(&pi, dl0s, dos, z);
            ampf[i][j].C = pi;
        };
    };

    OUT = fopen(file, "w");
    if (iampf == 1)
    {
        for (j = 0; j < nampf; j++)
            for (i = 0; i < nampf; i++)
            {
                amp = (ampf[i][j].E.a * ampf[i][j].E.b);
                e = 2 * fabs(amp) / (amp * amp + 1.);

                fprintf(OUT, "%d %lf %lf\t%lf %lf %lf %.1lf\n", i, ampf[i][j].C.x,
                        ampf[i][j].C.y, e, e, 0., 0.);
            };
    }
    else
    {
        for (j = 0; j < nampf; j++)
            for (i = 0; i < nampf; i++)
            {
                amp = (ampf[i][j].E.a * ampf[i][j].E.b);
                e = 2 * fabs(amp) / (amp * amp + 1.);
                fprintf(OUT, "%d %lf %lf\t%lf %lf %lf %.1lf\n", i, ampf[i][j].C.x,
                        ampf[i][j].C.y, e, e / 2., RTD*ampf[i][j].E.theta, 0.);
            };
    }

    fclose(OUT);
}
