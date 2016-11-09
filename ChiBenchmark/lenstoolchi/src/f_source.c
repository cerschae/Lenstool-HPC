#include<stdio.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*      nom:        f_source            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/* lecture de la premiere galaxy dans un fichier        */
/****************************************************************/

void    f_source(char *name, struct galaxie *galaxy, int *n)
{
    const extern  struct  g_mode          M;
    FILE    *IN;
    int     i, j, ncom;
    char    word[10], line[1280];

    i = j = 0;

    NPRINTF(stderr, "READ: %s ", name);
    IN = fopen(name, "r");
    if ( IN != NULL )
    {
        rewind(IN);
        for (j = 0, fscanf(IN, "%s", word); (strncmp(word, "#", 1) == 0); j++)
        {
            flire(IN, line);
            fscanf(IN, "%s", word);
        };

        ncom = j;
        rewind(IN);

        for (j = 0; j < ncom; j++, flire(IN, line));

        while (( fscanf(IN, "%s%lf%lf%lf%lf%lf%lf%lf",
                        galaxy[i].n, &galaxy[i].C.x, &galaxy[i].C.y,
                        &galaxy[i].E.a, &galaxy[i].E.b,
                        &galaxy[i].E.theta, &galaxy[i].z, &galaxy[i].I0)) != -1)
        {
            galaxy[i].E.theta *= DTR;
            flire(IN, line);
            i++;
        }

        NPRINTF(stderr, ": %d galaxy\n", i);
    }

    *n = i;
}
