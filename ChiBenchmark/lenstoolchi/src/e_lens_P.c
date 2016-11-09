#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        e_lens              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/
/* Return the number of images for a given source position.
 */

int e_lens_P(struct point ps, struct point pim[NIMAX], double dlsds)
{
    const extern struct point gsource_global[NGGMAX][NGGMAX];
   
    struct  bitriplet   Tsol[NIMAX];
   
    //struct    chaine *parent;

    int nimage;

    //Tsol=NULL;
    nimage = inverse(gsource_global, &ps, Tsol);

//  if ((nimage>0)&&(nimage<20))
    nimage = e_test_P(Tsol, nimage, &ps, pim, dlsds, 0.00001);
    /*  else
        {
            // free Tsol list until it is NULL value. (usefull in case nimage > NIMAX)
            // Go down to the last element of the Tsol list.
            while( Tsol != NULL )
            {
                parent = Tsol;
                Tsol=(Tsol->F);
                free((struct chaine*) parent);
            }
        }*/

    return(nimage);
}
