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
/* Format of the catalog:
 * ID  RA DEC a b theta z mag deps dtheta
 */

void    f_shape3(long int *istart, struct galaxie *liste, char *name)
{
    const extern  struct  g_mode          M;
    FILE    *IN;
    long int i, k = 0;
    char    line[128];

    i = (*istart);

    NPRINTF(stderr, "READ3: %s:", name);
    IN = fopen(name, "r");
    if ( IN != NULL )
        while ((fscanf(IN, "%s%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                       liste[i].n, &liste[i].C.x, &liste[i].C.y,
                       &liste[i].E.a, &liste[i].E.b,
                       &liste[i].E.theta, &liste[i].z, &liste[i].mag,
                       &liste[i].var1, &liste[i].var2)) != -1)
        {
            liste[i].E.theta *= DTR;
            liste[i].var2 *= DTR;  // dtheta
            fgets(line, 128, IN);
            i++;
            k++;
        }

    fclose(IN);

    NPRINTF(stderr, "%ld\n", k);

    (*istart) = i;
}
