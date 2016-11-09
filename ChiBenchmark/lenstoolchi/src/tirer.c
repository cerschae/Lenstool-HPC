#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        tirer               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/


void    tirer(struct galaxie *source)
{
    const extern struct pot       lens[];
    const extern struct g_mode    M;
    const extern struct g_source  S;
    const extern struct g_grille  G;
    const extern struct g_frame   F;
    const extern struct point gsource_global[NGGMAX][NGGMAX];

    int  k = S.rand;
    double dlsds;
    double t, I, J, N;//z,
    double Txmin, Txmax, Tymin, Tymax, Trmax;
    struct point pi, ps1, ps2;
    double dx, dy, ex, ey, e2;
    long int i, j, l, ng;
    double I0max = 100.;

    ng = G.ngrid;

    /* determination de la fenetre de tirage dans le plan soure */

    dlsds = dratio(lens[0].z, S.zs);
    pi.x = F.xmin;
    pi.y = F.ymin;
    e_dpl(&pi, dlsds, &ps1);
    Txmin = ps1.x;
    Tymin = ps1.y;
    pi.x = F.xmax;
    pi.y = F.ymax;
    e_dpl(&pi, dlsds, &ps2);
    Txmax = ps2.x;
    Tymax = ps2.y;
    Trmax = dist(ps1, ps2) / 2.;

    if (S.grid == 0)
    {
        NPRINTF(stderr, "COMP: Sources catalog drawn randomly\n");
        if ( S.ns > NFMAX )
        {
            fprintf( stderr, "ERROR: Too many sources. NFMAX limit is %d.\n", NFMAX);
            exit(1);
        }
        dx = Txmax - Txmin;
        dy = Tymax - Tymin;
        for (i = 0; i < S.ns; i++)
        {
            fprintf(stderr, "Create source %ld/%ld\r", i + 1, S.ns);
            sprintf(source[i].n, "%ld", i + 1);
            source[i].C.x = d_random(&k) * dx + Txmin;
            source[i].C.y = d_random(&k) * dy + Tymin;
            if (S.distz == 0)
            {
                source[i].z = S.zs;
                dratio_gal(&source[i], lens[0].z);
                source[i].I0 = I0max * pow(d_random(&k), 2.);
                t = S.taille * (1. + d_random(&k) / 2.);
                source[i].mag = 25. - 2.5 * log10(d_random(&k) * t * 100.);
            }
            else if (S.distz == 1)
            {
                source[i].z = S.zsmin + (S.zsmax - S.zsmin) * d_random(&k);
                dratio_gal(&source[i], lens[0].z);
                source[i].I0 = I0max * pow(d_random(&k), 2.) / (1 + source[i].z);
                t = S.taille * (1. + d_random(&k) / 2.) / (1. + source[i].z);
                source[i].mag = 21. - 2.5 *
                                log10(d_random(&k) * t / pow(1. + source[i].z, 5.));
            }
            else
            {
                source[i].z = d_rndz(S.zsmax, &k);
                dratio_gal(&source[i], lens[0].z);
                source[i].magabs = d_rndschechter(&k);
                source[i].type = d_rndtype(&k);
                source[i].I0 = I0max * pow(d_random(&k), 2.) / (1 + source[i].z);
                t = S.taille * (1. + d_random(&k) / 2.) / (1. + source[i].z);
                s_compmag(&source[i], &k);
            };

            // TODO: Understand this distribution and possibly edit:
            // S.emax = 0.3 but is defined in set_default.c and by the user using the elip_max keyword
            do
            {
                ex = d_gauss(S.emax, &k);  // the e1 and e2 components are ~Gaussian at low ell
                ey = d_gauss(S.emax, &k);  // on the edges, the J11+J22 of small objects enforce
                // the tails, hence produce Cauchy-Lorentz distrib
                e2 = sqrt(ex * ex + ey * ey);
            }
            while ( e2 >= 0.999 );
            // t = galaxy size (taille)
            // e = (a*a-b*b)/(a*a+b*b)
            source[i].E.a = t * sqrt(1. + e2);
            source[i].E.b = t * sqrt(1. - e2);

            if (S.emax == 0.)
                source[i].E.theta = 0.;
            else
                source[i].E.theta = 0.5 * atan2(ey, ex);  // is this definition arbitrary?
        }
    }
    else if (S.grid > 0)
    {
        //TODO: Fix the regular grid part
        NPRINTF(stderr, "COMP: Sources catalog drawn on a regular grid\n");
        if (S.grid == 1)
        {
            if (ng > 99)
                ng = 99;
            else
                ng = S.ns;
        }
        else
        {
            for (i = 0; i < ng; i++)
            {
                fprintf(stderr, "Create sources on line %ld...\r", i);
                for (j = 0; j < ng; j++)
                {
                    l = j + i * ng;
                    I = i;
                    J = j;
                    N = (double)(ng - 1);
                    I = I / N;
                    J = J / N;
                    Txmin = Txmin + .0001;
                    Tymin = Tymin + .0001;
                    sprintf(source[l].n, "%ld", l + 1);
                    source[i].z = S.zs;
                    dratio_gal(&source[i], lens[0].z);
                    if (G.pol == 0)
                    {
                        source[l].C.x = Txmin + I * (Txmax - Txmin);
                        source[l].C.y = Tymin + J * (Tymax - Tymin);
                    }
                    else
                    {
                        I = I + 0.0001;
                        source[l].C.x = I * Trmax / 2.*cos(J * 2.*PI);
                        source[l].C.y = I * Trmax / 2.*sin(J * 2.*PI);
                    };
                    if (S.grid == 2)
                    {
                        source[l].C.x = gsource_global[i][j].x + .000;
                        source[l].C.y = gsource_global[i][j].y + .000;
                    };

                    source[l].E.a = S.taille;
                    source[l].E.b = S.taille;
                    source[l].E.theta = 0.;
                    source[l].I0 = 50.;
                };
            };
        }
    }
    fprintf(stderr, "\n");

}
