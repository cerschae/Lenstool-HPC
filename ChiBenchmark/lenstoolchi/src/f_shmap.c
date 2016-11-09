#include<stdio.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*      nom:        f_shmap             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       19/12/96            */
/*      place:      Toulouse            */
/* Read shear map files - as created by smg -           */
/****************************************************************/

void    f_shmap(int *istart, struct shear *liste, char *name)
{
    const extern  struct  g_mode          M;
    FILE    *IN;
    int i, j, ncom;
    double  x;
    char    word[10], line[256];

    i = 0;

    NPRINTF(stderr, "READ SHEAR MAP: %s ", name);
    IN = fopen(name, "r");

    if ( IN != NULL )
    {

        /* skipping the header */

        rewind(IN);
        for (j = 0, fscanf(IN, "%s", word); (strncmp(word, "#", 1) == 0); j++)
        {
            flire(IN, line);
            fscanf(IN, "%s", word);
        };
        ncom = j;
        rewind(IN);
        for (j = 0; j < ncom; j++, flire(IN, line));

        /* reading the data */

        while ((fscanf(IN, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                       &liste[i].n, &liste[i].C.x, &liste[i].C.y,
                       &x, &x, &x, &liste[i].mx, &liste[i].my, &liste[i].dx, &liste[i].dy,
                       &liste[i].err)
               ) != -1)
        {
            flire(IN, line);
            i++;
        };
    }

    fclose(IN);

    NPRINTF(stderr, ": %d shear points\n", i);

    (*istart) = i;
}
