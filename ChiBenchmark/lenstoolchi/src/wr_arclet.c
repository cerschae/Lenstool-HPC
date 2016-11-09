#include<stdio.h>
#include<math.h>
#include "fonction.h"
#include "constant.h"
#include"dimension.h"
#include "structure.h"

/****************************************************************/
/*      nom:        ecrire              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/
//NOT USED IN LENSTOOL

void wr_arclet(int nl, struct galaxie liste[NAMAX], char name[50])
{
    FILE    *OUT;
    int i;

    OUT = fopen(name, "w");

    for (i = 0; i < nl; i++)
    {
        if (liste[i].E.theta > PI)
            liste[i].E.theta -= PI;
        if (liste[i].c == 's')
        {
            liste[i].E.a = liste[i].E.b = 0.;
        };

        fprintf(OUT, "%4s %8.3lf %8.3lf %.5lf %.5lf %8.3lf %5.3lf %6.3lf\n",
                liste[i].n, liste[i].C.x,
                liste[i].C.y, liste[i].E.a, liste[i].E.b,
                RTD*liste[i].E.theta, liste[i].z, liste[i].mag);
    };

    fclose(OUT);
}
