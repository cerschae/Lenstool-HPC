#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        w_prop          */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    w_prop(int nprop, char propfile[])
{
    const extern  struct  g_frame     F;
    const extern struct   propertie   gprop[NTMAX][NTMAX];

    FILE    *OUT;
    register int    i, j;

    OUT = fopen(propfile, "w");

    fprintf(OUT, "%d %d %lf %lf %lf %lf \n",
            nprop, nprop, F.xmin, F.xmax, F.ymin, F.ymax);

    for (j = 0; j < nprop; j++)
        for (i = 0; i < nprop; i++)
        {
            fprintf(OUT, "%.3lf %.4lf %.4lf %.2lf %.4lf %.4lf %.4lf\n",
                    gprop[i][j].A, gprop[i][j].k, gprop[i][j].s, 180. / PI*gprop[i][j].theta,
                    gprop[i][j].t, gprop[i][j].q, gprop[i][j].g);
        };
    fclose(OUT);
}
