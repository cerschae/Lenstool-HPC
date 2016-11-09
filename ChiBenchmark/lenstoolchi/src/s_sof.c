#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

static void wr_sof(int nl, struct galaxie *ima, struct galaxie *source, char name[50]);

/****************************************************************/
/*      nom:        s_sof           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    s_sof()
{
    extern  struct  g_mode      M;
    long int nima;
    struct  galaxie     ima[NAMAX], sofim[NAMAX];

    nima = 0;
    if (M.sof == 1)
        f_shape(&nima, ima, M.imsfile,0);
    else if (M.sof == 2)
        f_shape2(&nima, ima, M.imsfile);

    e_unlens_fast(nima, ima, sofim);

    wr_sof(nima, ima, sofim, M.sfile);

}

static void wr_sof(int nl, struct galaxie *ima, struct galaxie *source, char name[50])
{
    FILE    *OUT;
    int i;
    double  jacobian;

    OUT = fopen(name, "w");

    for (i = 0; i < nl; i++)
    {
        if (source[i].E.theta > PI)
            source[i].E.theta -= PI;
        if (ima[i].thp > PI)
            ima[i].thp -= PI;

        if (ima[i].dis != 0.)
        {
            jacobian = ima[i].dp - ima[i].tp * ima[i].taux / ima[i].dis;
        }
        else
        {
            jacobian = 0.;
        }

        fprintf(OUT, "%s %.3lf %.3lf %.5lf %.5lf %.3lf %.3lf %.2lf %.4lf %.4lf %.4lf %.3lf %.4lf\n",
                source[i].n, source[i].C.x,
                source[i].C.y, source[i].E.a, source[i].E.b,
                RTD*source[i].E.theta, source[i].z, source[i].mag,
                ima[i].taux, ima[i].tauy, ima[i].tp, ima[i].thp*RTD, jacobian);
    };

    fclose(OUT);
}
