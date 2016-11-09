#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        g_radial            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    g_radial(int iradial, double zradial, double theta)
{
    const extern  struct  g_mode          M;
    const extern  struct  g_frame F;
    const extern  struct  pot lens[];

    struct  polar   P;
    struct  point   Q, QS;
    struct  ellipse ampli;
    struct  matrix  MA;
    double  alpha;
    double  dl0s, dos, dlsds;
    double  N, A;
    double  qpot, epot, dpot, tpot;
    double  kappa, gam, ga1, ga2, gp;
    double  kappa0 = 1.;
    FILE    *OUT, *OUT2;
    register    int i;

    NPRINTF(stderr, "COMP: radialprop %d\n", iradial);

    dl0s = distcosmo2( lens[0].z, zradial );
    dos = distcosmo1( zradial );
    dlsds = dl0s / dos;
    OUT = fopen("radial.dat", "w");
    OUT2 = fopen("radial2.dat", "w");
    
    switch (iradial)
    {
        case(1):
            P.theta = theta;
            N = NPOINT;
            P.r = 0.0001;
            for (i = 0; i < NPOINT && P.r < F.rmax; i++)
            {
                P.r *= 1.05;

                if (theta > 0)
                {
                    Q.x = lens[0].C.x + P.r * cos(P.theta);
                    Q.y = lens[0].C.y + P.r * sin(P.theta);
                    ampli = e_unmag(&Q, dl0s, dos, zradial);
                    A = ampli.a * ampli.b;
                    fprintf(OUT, "%.4lf %.3lf %lf %lf %lf %lf %lf %lf\n",
                            P.r, P.theta, ampli.a, ampli.b, ampli.theta, A, 1. - (ampli.a + ampli.b) / 2.,
                            (ampli.a - ampli.b) / 2.);

                    MA = e_grad2(&Q, dl0s, zradial);
                    MA.a /= dos;
                    MA.b /= dos;
                    MA.c /= dos;
                    kappa = (MA.a + MA.c) / 2.;
                    ga1 = (MA.a - MA.c) / 2.;
                    ga2 = MA.b;
                    gam = sqrt(ga1 * ga1 + ga2 * ga2);
                    gp = gam / (1 - kappa);

                    if (i == 0)
                        kappa0 = kappa;

                    fprintf(OUT2, "%lf %.2lf %lf %g %lf %lf %lf %lf\n",
                            P.r, P.theta, kappa, kappa / kappa0, ga1, ga2, gam, gp);
                }
                else
                {
                    kappa = ga1 = ga2 = gam = gp = 0;

                    for (i = 0; i < 60; i++)
                    {
                        P.theta = 6 * i / 180 * PI;
                        Q.x = lens[0].C.x + P.r * cos(P.theta);
                        Q.y = lens[0].C.y + P.r * sin(P.theta);
                        MA = e_grad2(&Q, dl0s, zradial);
                        MA.a /= dos;
                        MA.b /= dos;
                        MA.c /= dos;
                        kappa += (MA.a + MA.c) / 2.;
                        ga1 += (MA.a - MA.c) / 2.;
                        ga2 += MA.b;
                    }
                    kappa /= 60;
                    ga1 /= 60;
                    ga2 /= 60;
                    gam = sqrt(ga1 * ga1 + ga2 * ga2);
                    gp = gam / (1 - kappa);

                    if (i == 0)
                        kappa0 = kappa;

                    fprintf(OUT2, "%lf %.2lf %lf %g %lf %lf %lf %lf\n",
                            P.r, 0.0, kappa, kappa / kappa0, ga1, ga2, gam, gp);
                }

            };
            break;
        case(2):
            P.theta = theta;
            N = NPOINT;
            for (i = 0; i < NPOINT; i++)
            {
                P.r = ((double)i) * F.rmax / N;
                Q.x = lens[0].C.x + P.r * cos(P.theta);
                Q.y = lens[0].C.y + P.r * sin(P.theta);
                ampli = e_unmag(&Q, dl0s, dos, zradial);
                A = ampli.a * ampli.b;
                qpot = fabs(ampli.b / ampli.a);
                epot = fabs(1. - qpot * qpot) / (1. + qpot * qpot);
                dpot = (1. + qpot * qpot) / 2. / qpot;
                tpot = fabs(1. - qpot * qpot) / 2. / qpot;
                fprintf(OUT, "%.3lf %.3lf %.3lf %.3lf %.3lf %.3lf \n", P.r,
                        epot, dpot, tpot, ampli.theta, A);
                fprintf(OUT2, "%lf %lf %lf\n", P.r, tpot*tpot / P.r,
                        tpot / P.r);
            };
            break;
        case(3):
            P.theta = theta;
            N = NPOINT;
            for (i = 0; i < NPOINT; i++)
            {
                P.r = ((double)(2 * i - NPOINT)) * F.rmax / N;
                Q.x = P.r * cos(P.theta);
                Q.y = P.r * sin(P.theta);
                e_dpl(&Q, dlsds, &QS);
                alpha = P.r / fabs(P.r) * dist(Q, QS);
                fprintf(OUT, "%lf %lf %lf %lf\n", P.r, alpha, P.r - lens[0].cr, P.r - lens[0].ct);
            };
            break;
        default:
            fprintf(stderr, "WARNING: radialprop mode not defined\n");
            break;
    };

    fclose(OUT);
}
