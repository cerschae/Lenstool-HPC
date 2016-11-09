#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        pro_arclet          */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse
 *
 * Append physical parameters from input information contained
 * in the gal list of galaxies to the gal list of galaxies.
 *
 * Parameters :
 * - n : number of elements in gal
 * - gal : galaxie strcture containing the arclet input informations
 */
void    pro_arclet( long int n, struct galaxie *gal )
{
    long int    i;
    double  qi;

    for (i = 0; i < n; i++)
    {
        if ( n > 1000 )
            printf( "INFO: prepare arclet %ld/%ld\r", i + 1, n);

        qi = gal[i].q = gal[i].E.b / gal[i].E.a;

        if (qi > 0.)
        {
            gal[i].tau = .5 * (1. / qi - qi);
            gal[i].dis = .5 * (1. / qi + qi);
            gal[i].eps = (1. - qi) / (1. + qi);
        }
        else
        {
            gal[i].tau = gal[i].E.a;
            gal[i].dis = sqrt(1. + gal[i].tau * gal[i].tau);
            gal[i].eps = gal[i].tau / gal[i].dis;
        }
    }
    if ( n > 1000 )    printf( "\n" );
}
