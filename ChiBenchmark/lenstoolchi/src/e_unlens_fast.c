#include<stdio.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        e_unlens            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    e_unlens_fast(  int nima,
                        struct galaxie *ima,
                        struct galaxie *source )
{
    const extern  struct  g_source    S;
    //const extern    struct  g_mode  M;
    const extern  struct  pot lens[];
    struct  matrix  MM;
    struct  ellipse ampli;
    register int    i;
    double  A, B, C, gg, qq;

    for (i = 0; i < nima; i++)
    {
        if (ima[i].z == 0.)
            ima[i].z = S.zs;
        ima[i].dr = dratio(lens[0].z, ima[i].z);
        source[i].c = ima[i].c;

        if (ima[i].dr > 0)
        {
            MM = e_grad2_gal(&ima[i], NULL);
            MM.a *= ima[i].dr;
            MM.b *= ima[i].dr;
            MM.c *= ima[i].dr;
            MM.d *= ima[i].dr;
            A = 1. - MM.a;
            B = -MM.b;
            C = 1. - MM.c;
            ampli = formeli(A, B, C);

            ima[i].kappa = (MM.a + MM.c) / 2. / ima[i].dr;
            ima[i].gamma1 = (MM.a - MM.c) / 2. / ima[i].dr;
            ima[i].gamma2 = MM.b / ima[i].dr;
            ima[i].thp = 0.5 * atan2(ima[i].gamma2, ima[i].gamma1);
            gg = ima[i].dr * sqrt(ima[i].gamma1 * ima[i].gamma1
                                  + ima[i].gamma2 * ima[i].gamma2) / (1 - ima[i].dr * ima[i].kappa);
            ima[i].tp = 2 * gg / (1 - gg * gg);
            ima[i].dp = sqrt(1. + ima[i].tp * ima[i].tp);
            ima[i].ep = ima[i].tp / ima[i].dp;

            qq = ima[i].E.b / ima[i].E.a;
            ima[i].tau = (1. - qq * qq) / 2. / qq;
            ima[i].taux = ima[i].tau * cos(2.*(ima[i].E.theta - ima[i].thp));
            ima[i].tauy = ima[i].tau * sin(2.*(ima[i].E.theta - ima[i].thp));
            ima[i].dis = sqrt(1. + ima[i].tau * ima[i].tau);
            ima[i].eps = ima[i].tau / ima[i].dis;

            if (ima[i].E.b != 0.)
                isoima(&ima[i].E, &ampli, &source[i].E);
            else
            {
                source[i].E.a = source[i].E.b =
                                    ima[i].E.a * fabs(ampli.a * ampli.b);
                source[i].E.theta = 0.;
            }

            e_dpl(&ima[i].C, ima[i].dr, &source[i].C);
            strcpy(source[i].n, ima[i].n);
            source[i].z = ima[i].z;
            source[i].mag = ima[i].mag - 2.5 * log10(fabs(ampli.a * ampli.b));
            source[i].dr = ima[i].dr;
        }
        else
            source[i] = ima[i];
    };

}
