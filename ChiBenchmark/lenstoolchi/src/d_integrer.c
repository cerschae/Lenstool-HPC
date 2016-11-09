#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#define ITERMAX 1

/****************************************************************/
/*      programme   d_integrer.c            */
/*      auteur      Henri Bonnet            */
/*      place       OMP             */
/*      date        12.11.1991          */
/*      version     1               */
/****************************************************************/

/****************************************************************/
/*  d_integrer(x,y,k,t,f)                   */
/*      double x,y,t,f;                 */
/*      int   k;                    */
/*                              */
/*  evalue, sur un pixel centre en (x,y), de cote t     */
/*  l'integrale de la fonction:             */
/*      f(x,y)=Io/(1+alpha**2*(e**2*(x-X0)**2+(y-Y0)**2))   */
/*  dont la valeur en (x,y) vaut f.             */
/*  avec une precision absolue prec/nobj            */
/*                              */
/*  les parametres sont definis dans "para"         */
/*      alpha=1/(e*l)                   */
/****************************************************************/


double d_integrer(  struct galaxie A,
                    double x, double y,
                    double t, double f,
                    int n )
{
    const extern struct g_observ O;

    double X, Y;
    double res = 0;
    double tab[4];
    int q = 0;
    double inf;
    struct point I, S;



    for (X = x - t; X < x + 3.*t / 2.; X += 2.*t)
        for (Y = y - t; Y < y + 1.5*t; Y += 2.*t)
        {
            I.x = X;
            I.y = Y;
            e_dpl(&I, A.dr, &S);
            tab[q++] = d_profil(S.x, S.y, &A);
        };

    res = (tab[0] + tab[1] + tab[2] + tab[3]) / 4.;

    if (n < ITERMAX && ( res != 0. || f != 0.))
    {
        if (f == 0.)
        {
            res = 0.;
            q = 0;
            for (X = x - t; X < x + 3.*t / 2.; X += 2.*t)
                for (Y = y - t; Y < y + 1.5*t; Y += 2.*t)
                {
                    I.x = (X + x) / 2.;
                    I.y = (Y + y) / 2.;
                    e_dpl(&I, A.dr, &S);
                    tab[q] = d_profil(S.x, S.y, &A);
                    tab[q] = d_integrer(A, X, Y, t / 2., tab[q], n + 1);
                    q++;
                };
            res = (tab[0] + tab[1] + tab[2] + tab[3]) / 4.;
        }

        else
        {
            if (res == 0.)
                inf = fabs(f);
            else
            {
                inf = fabs(res);
                if (fabs(res) < fabs(f))
                    inf = fabs(res);
            };

            if (inf != 0.)
                if (fabs(res - f) / inf > O.prec)
                {
                    q = 0;
                    for (X = x - t; X < x + 3.*t / 2.; X += 2.*t)
                        for (Y = y - t; Y < y + 1.5*t; Y += 2.*t)
                        {
                            I.x = (X + x) / 2.;
                            I.y = (Y + y) / 2.;
                            e_dpl(&I, A.dr, &S);
                            tab[q] = d_profil(S.x, S.y, &A);
                            tab[q] = d_integrer(A, X, Y, t / 2., tab[q], n + 1);
                            q++;
                        };
                    res = (tab[0] + tab[1] + tab[2] + tab[3]) / 4.;
                };
        };
    };
    return((f + res) / 2.);

}
