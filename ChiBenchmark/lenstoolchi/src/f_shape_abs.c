#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*      nom:        f_shape_abs             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       23/04/04            */
/*      place:      Toulouse            */
/* lecture de fichiers ellipses, source ou image
 *
 * Read a data file of elliptical regions defined in absolute coordinates
 * Input format :
 *   char* n, double Cx, double Cy, double a, double b, double theta, double z, flaot mag
 *
 * Cx,Cy,a,b and theta must be in degree.
 *
 * If the first character of the line is a # then the line is ignored.
 *
 * Parameters :
 * - istart : number of objects already counted and to increment
 * - liste : list of galaxie structure to fill with the data
 * - name : name of the file to read
 *
 */

void    f_shape_abs(long int *istart, struct galaxie *liste, char *name)
{
    const extern  struct  g_mode          M;
    FILE    *IN;
    long int     i, k = 0;
    char    line[128];

    i = (*istart);

    NPRINTF(stderr, "READ_ABS: %s:", name);
    IN = fopen(name, "r");

    while ( IN != NULL && !feof(IN) && !ferror(IN) )
    {
        flire(IN, line);
        if ( sscanf(line, "%s%lf%lf%lf%lf%lf%lf%lf",
                    liste[i].n, &liste[i].C.x, &liste[i].C.y,
                    &liste[i].E.a, &liste[i].E.b,
                    &liste[i].E.theta, &liste[i].z, &liste[i].mag) == 8 )
        {
            if ( liste[i].n[0] != '#' )
            {
                liste[i].C.x -= M.ref_ra;
                liste[i].C.x *= -3600 * cos(M.ref_dec * DTR);
                liste[i].C.y -= M.ref_dec;
                liste[i].C.y *= 3600;
                liste[i].E.theta *= DTR;
                if ((liste[i].E.a == 0.) || (liste[i].E.b == 0.))
                    liste[i].c = 's';
                else
                    liste[i].c = 'g';

                i++;
                k++;
            }
        }
    }

    if ( IN == NULL || ferror(IN) || k == 0 )
    {
        fprintf( stderr, "ERROR: Error reading the %s file\n", name);
        exit(-1);
    }

    fclose(IN);

    NPRINTF(stderr, "%ld\n", k);

    (*istart) = i;
}
