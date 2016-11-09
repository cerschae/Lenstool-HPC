#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/*
* local
* Jean-Paul Kneib
* jan 95
* Cambridge
*/

void    local(int type, char localfile[])
{
    const extern  struct  g_source    S;
    const extern  struct  pot lens[];

    FILE    *OUT;
    struct  galaxie alet[NAMAX];
    long int na;

    struct  ellipse ampli;
    struct  matrix  M;
    double  q, eps, tau, epsx, epsy, taux, tauy, ka, ga1, ga2, ga;

    register long int    i;


    na = 0;
    if (type > 0)
        f_shape(&na, alet, localfile, 1);
    else
        f_shape2(&na, alet, localfile);

    OUT = fopen("local.dat", "w");

    for (i = 0; i < na; i++)
    {
        if (alet[i].z <= lens[0].z)
            alet[i].dl0s = distcosmo2( lens[0].z, S.zs );
        else
            alet[i].dl0s = distcosmo2( lens[0].z, alet[i].z );

        alet[i].dos = distcosmo1( alet[i].z );
        alet[i].dr = alet[i].dl0s / alet[i].dos;

        ampli = e_unmag_gal(&alet[i]);
        q = alet[i].E.b / alet[i].E.a;
        eps = (1. - q * q) / (1. + q * q);
        tau = .5 * (1. / q - q);
        epsx = -eps * cos(2.*(alet[i].E.theta - ampli.theta));
        epsy = eps * sin(2.*(alet[i].E.theta - ampli.theta));
        taux = -tau * cos(2.*(alet[i].E.theta - ampli.theta));
        tauy = tau * sin(2.*(alet[i].E.theta - ampli.theta));

        M = e_grad2_gal(&alet[i], NULL);
        M.a /= alet[i].dos;
        M.b /= alet[i].dos;
        M.c /= alet[i]. dos;
        ka = (M.a + M.c) / 2.;
        ga1 = (M.a - M.c) / 2.;
        ga2 = M.b;
        ga = sqrt(ga1 * ga1 + ga2 * ga2);

        fprintf(OUT,
                "%s %.3lf %.3lf   %.3lf %.3lf %.3lf    %.3lf %.3lf %.3lf    %.3lf %.3lf %.3lf\n",
                alet[i].n, alet[i].C.x, alet[i].C.y,
                eps, epsx, epsy, tau, taux, tauy, ka, ga1, ga2
               );
    }

    fclose(OUT);
}
