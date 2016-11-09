#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        ecrirei             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    w_gianti(int ngs, char name[20])
{
    const extern  struct  pointgal    gianti[][NIMAX];
    register    int i, j;
    int l;
    FILE    *OUT;

    OUT = fopen(name, "w");

    l = 0;
    for (i = 0; i < ngs; i++)
        for (j = 0; (j < NIMAX) && (gianti[i][j].n != 0); l++, j++);

    fprintf(OUT, "%d\n", 1);
    fprintf(OUT, "%d\n", l);
    fprintf(OUT, "double txt complex\n");
    fprintf(OUT, "test\n");
    for (i = 0; i < ngs; i++)
        for (j = 0; (j < NIMAX) && (gianti[i][j].n != 0); j++)
            fprintf(OUT, "%lf %lf\n", gianti[i][j].C.x, gianti[i][j].C.y);

    fclose(OUT);
}
