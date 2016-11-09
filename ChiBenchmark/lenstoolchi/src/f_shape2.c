#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        f_shape2                */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/* lecture de fichiers ellipses, source ou image        */
/****************************************************************/

void    f_shape2(long int *istart, struct galaxie *liste, char *name)
{
    const extern  struct  g_mode          M;
    FILE    *IN;
    long int i, k = 0;
    char    line[128];

    i = (*istart);

    NPRINTF(stderr, "READ2: %s:", name);
    IN = fopen(name, "r");
    if ( IN != NULL )
        while ((fscanf(IN, "%s%lf%lf%lf%lf%lf%lf%lf",
                       liste[i].n, &liste[i].C.x, &liste[i].C.y,
                       &liste[i].E.a, &liste[i].E.b,
                       &liste[i].E.theta, &liste[i].z, &liste[i].mag)) != -1)
        {
            liste[i].E.theta *= DTR;
            if ((liste[i].E.a == 0.) || (liste[i].E.b == 0.))
                liste[i].c = 's';
            else
                liste[i].c = 'g';
            fgets(line, 128, IN);
            i++;
            k++;
        }

    fclose(IN);

    NPRINTF(stderr, "%ld\n", k);

    (*istart) = i;
}
