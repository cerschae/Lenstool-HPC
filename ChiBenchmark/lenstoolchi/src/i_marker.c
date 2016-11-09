#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*      nom:        i_marker            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    i_marker(char markfile[], double z)
{
    const extern  struct  pot lens[];
    int n = 0;
    char    line[128];
    double  dlsds;
    struct  point   A, B;
    FILE    *IN, *OUT;

    dlsds = dratio(lens[0].z, z);

    IN = fopen(markfile, "r");
    OUT = fopen("marker_s.dat", "w");

    if ( IN != NULL )
    {
        while ((fscanf(IN, "%d%lf%lf", &n, &A.x, &A.y) != -1))
        {
            flire(IN, line);
            e_dpl(&A, dlsds, &B);
            fprintf(OUT, "%d %lf\t%lf\t\n", n, B.x, B.y);
        };
    };

    fclose(IN);
    fclose(OUT);
}
