#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

struct propertie gprop[NTMAX][NTMAX];

/****************************************************************/
/*      nom:        g_prop              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    g_prop(int nprop, double z)
{

    const extern  struct  g_cosmo     C;
    const extern  struct  g_frame     F;
    extern  struct  g_mode      M;
    const extern  struct  pot lens[];
    extern struct point gsource_global[NGGMAX][NGGMAX];
   
    register int    i, j;
    double  dl0s, dos, dlsds;
    double  qpot;
    struct  point   Centre, pi, ps;
    struct  ellipse ampli;
    double  mass;

    NPRINTF(stderr, "COMP: properties grid\n");

    Centre.x = 0.;
    Centre.y = 0.;
    dl0s = distcosmo2( lens[0].z, z);
    dos = distcosmo1( z );
    dlsds = dl0s / dos;
    e_unlensgrid(gsource_global, dlsds);
    M.masse = 0.;

    for (i = 0; i < nprop; i++)
    {
        pi.x = i * (F.xmax - F.xmin) / (nprop - 1) + F.xmin;
        for (j = 0; j < nprop; j++)
        {
            pi.y = j * (F.ymax - F.ymin) / (nprop - 1) + F.ymin;
            e_dpl(&pi, dlsds, &ps);
            ampli = e_unmag(&pi, dl0s, dos, z);

            gprop[i][j].P.x = pi.x;
            gprop[i][j].P.y = pi.y;
            gprop[i][j].k = 1. - (ampli.a + ampli.b) / 2.;
            gprop[i][j].s = (ampli.a - ampli.b) / 2.;
            gprop[i][j].theta = ampli.theta;
            qpot = fabs(ampli.b / ampli.a);
            gprop[i][j].q = 1. / qpot;
            gprop[i][j].e = (1. - qpot * qpot) / (1. + qpot * qpot);
            gprop[i][j].d = (1. + qpot * qpot) / 2. / qpot;
            gprop[i][j].t = fabs(1. / qpot - qpot) / 2. / qpot;
            gprop[i][j].A = 1. / fabs(ampli.a * ampli.b);
            gprop[i][j].g = dist(ps, pi);

            if (dist(Centre, pi) <= M.radius)
                M.masse += 1. - 0.5 * (ampli.a + ampli.b);

        };
    };

    mass = MCRIT / C.h * M.radius * M.radius *
           distcosmo1(z) * distcosmo1(lens[0].z) / distcosmo2(lens[0].z, z);

    NPRINTF(stderr, "\n masse dans le rayon %lf : %lf\n", M.radius, M.masse);
    NPRINTF(stderr, "\n D_ratio : %lf\n",
            distcosmo1(z)*distcosmo1(lens[0].z) / distcosmo2(lens[0].z, z));
    NPRINTF(stderr, "\n masse_crit dans le rayon (M_sol) : %e\n\n", mass);
}
