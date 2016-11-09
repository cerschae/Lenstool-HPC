#include<stdio.h>
#include<math.h>
#include "fonction.h"
#include "constant.h"
#include"dimension.h"
#include "structure.h"

/****************************************************************/
/*      nom:        gs_ecrire           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void gs_ecrire( char name[20])
{
    const extern  struct  g_grille    G;
    const extern  struct  point   gsource[][NGGMAX];

    FILE    *OUT1;
    register    int i, j;

    OUT1 = fopen(name, "w");

    for (i = 0; i < G.ngrid; i++)
    {
        for (j = 0; j < G.ngrid; j++)
        {
            fprintf(OUT1, "%d\t%d\t%lf\t%lf\t\n", i, j, gsource[i][j].x, gsource[i][j].y);
        };
        if (++i < G.ngrid)
        {
            for (j = G.ngrid - 1; j >= 0; j--)
            {
                fprintf(OUT1, "%d\t%d\t%lf\t%lf\t\n", i, j, gsource[i][j].x, gsource[i][j].y);
            };
        };
    };

    fclose(OUT1);
}
