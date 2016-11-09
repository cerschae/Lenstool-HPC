#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*  g_profil(x,y,k)                     */
/*      double x,y,k;                   */
/*                              */
/*  evalue, en (x,y), la fonction:              */
/*      f(x,y)=Io/(1+alpha**2*(e**2*(x-C.x)**2+(y-C.y)**2)) */
/*  k est l'indice de l'objet qu'on integre         */
/*                              */
/*  les parametres sont definis dans "para"         */
/*      alpha=1/(e*l)                   */
/****************************************************************/

double g_profil(double x, double y, struct galaxie gal)
{
    double xx, yy, xxx, yyy;
    double res;

    xx = ((x - gal.C.x) * sin(gal.E.theta) - (y - gal.C.y) * cos(gal.E.theta)) / gal.E.b;
    yy = ((x - gal.C.x) * cos(gal.E.theta) + (y - gal.C.y) * sin(gal.E.theta)) / gal.E.a;
    xxx = xx * xx;
    yyy = yy * yy;

    /*
    res=gal.I0*exp(-(xxx+yyy));
    */

    res = gal.I0 / (1. + (xxx + yyy));

    return(res);
}
