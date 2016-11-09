#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        e_tau               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse
 * Fill the tau parameter of the galaxies in the gal list.
 * Parameters :
 * - n : number of elements in gal
 * - gal : list of galaxie structure
 ****************************************************************/

void  e_tau(long int n, struct galaxie gal[NAMAX])
{
    long int    i;
    double  qi;

    for (i = 0; i < n; i++)
    {
        qi = gal[i].q = gal[i].E.b / gal[i].E.a;
        gal[i].tau = .5 * (1. / qi - qi);
    };
}
