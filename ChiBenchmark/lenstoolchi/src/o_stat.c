#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/********************************************************/
/*      fonction: o_stat            */
/*      auteur: jpk             */
/********************************************************/

void    o_stat(int na, struct galaxie arclet[NAMAX])
{
    register int    i;
    double  Dtauxs, Etauxs, Mtauxs, tauxs;
    double  Dtauys, Etauys, Mtauys, tauys;
    double  chi0;
    FILE    *OUT;

    o_mag(na, arclet);

    Mtauxs = 0.;
    Etauxs = 0.;
    Mtauys = 0.;
    Etauys = 0.;
    for (i = 0; i < na; i++)
    {
        tauxs = arclet[i].tau * cos(2.*(arclet[i].E.theta - arclet[i].thp));
        tauxs = tauxs * arclet[i].dp - arclet[i].dis * arclet[i].tp;
        Mtauxs += tauxs;
        Etauxs += tauxs * tauxs;
        tauys = arclet[i].tau * sin(2.*(arclet[i].E.theta - arclet[i].thp));
        Mtauys += tauys;
        Etauys += tauys * tauys;
    }
    Mtauxs = Mtauxs / na;
    Mtauys = Mtauys / na;
    Etauxs = Etauxs / na;
    Etauys = Etauys / na;
    Dtauxs = Etauxs - Mtauxs * Mtauxs;
    Dtauys = Etauys - Mtauys * Mtauys;

    chi0 = Etauxs + Etauys;

    OUT = fopen("stat.dat", "w");

    fprintf(OUT, "%lf %lf %lf\n", Mtauxs, Dtauxs, Etauxs);
    fprintf(OUT, "%lf %lf %lf\n", Mtauys, Dtauys, Etauys);
    fprintf(OUT, "%lf\n", chi0);
}
