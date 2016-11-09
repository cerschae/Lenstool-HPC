#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        weight_baryc                */
/*      auteur:     Ghislain Golse          */
/*      date:       12/99               */
/*      place:      Toulouse            */
/****************************************************************/


/*
* barycenter of a list of points weighted by the amplification
*
* Global variables used :
* - amplifi
* - in e_amp() : G, lens, lens_table
* - in amplif() : multi, I, amplifi, G, lens, lens_table, C
*/
struct  point   weight_baryc(struct point *P, struct galaxie *multi, int n, int n_famille)
{
//    const extern  double   amplifi[NFMAX][NIMAX];


    struct  point   B;
    register int i;
    double    dlsds, Atot, A;

    dlsds = multi[0].dr;
    B.x = B.y = 0.;
    Atot = 0.;

    for (i = 0; i < n; i++)
    {
        /* A=1./fabs(e_amp(multi[i].C,dlsds)); */
//        A = fabs(amplifi[n_famille][i]);
        A = 1./fabs(e_amp_gal(&multi[i], NULL));
        Atot += sqrt(A);
        B.x += P[i].x * sqrt(A);
        B.y += P[i].y * sqrt(A);
    };

    B.x /= Atot;
    B.y /= Atot;

    return(B);
}

/* printf("B.x=%.3lf B.y=%.3lf\n",B.x,B.y); */
