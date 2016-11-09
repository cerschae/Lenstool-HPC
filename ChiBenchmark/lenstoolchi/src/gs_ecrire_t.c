#include<stdio.h>
#include<math.h>
#include "fonction.h"
#include "constant.h"
#include"dimension.h"
#include "structure.h"

/****************************************************************/
/*      nom:        gs_ecrire_t         */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void gs_ecrire_t( char name[20])
{
    const extern  struct  g_grille    G;
    const extern  struct  point   gsource[][NGGMAX];

    FILE    *OUT1;
    register    int i, j;

    OUT1 = fopen(name, "w");

    for (j = 0; j < G.ngrid; j++)
    {
        for (i = 0; i < G.ngrid; i++)
        {
            fprintf(OUT1, "%d\t%d\t%lf\t%lf\t\n", i, j, gsource[i][j].x, gsource[i][j].y);
        };
        if (++j < G.ngrid)
        {
            for (i = G.ngrid - 1; i >= 0; i--)
            {
                fprintf(OUT1, "%d\t%d\t%lf\t%lf\t\n", i, j, gsource[i][j].x, gsource[i][j].y);
            };
        };
    };

    fclose(OUT1);
}
