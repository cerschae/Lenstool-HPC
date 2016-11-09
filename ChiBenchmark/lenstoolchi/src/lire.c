#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        lire                */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

int lire(   struct galaxie *liste,
            char *name )
{
    FILE    *IN;
    int i = 0;

    printf("lecture de %s\n", name);
    IN = fopen(name, "r");
    if (  IN != NULL )
        while ((fscanf(IN, "%s%lf%lf%lf%lf%lf%lf\n",
                       liste[i].n, &liste[i].C.x, &liste[i].C.y,
                       &liste[i].E.a, &liste[i].E.b,
                       &liste[i].E.theta, &liste[i].z)) != -1)
            liste[i++].E.theta *= DTR;

    fclose(IN);
    return(i);
}
