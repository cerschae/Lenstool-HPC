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
/****************************************************************
 * For source, use the 2 grids gsource and gimage that establish a
 * bijection between the source and the image plane to find a list
 * of arclets. The number of arclets must be lower than 20.
 *
 * Fill the image list of arclets and return the number of arclets.
 *
 * Parameters :
 * - source : the source
 * - image  : a list of arclets for the familly (modified)
 */
int e_lens( struct  galaxie source,
            struct  galaxie image[NIMAX] )
{
    const extern  struct  g_mode          M;
    const extern struct point gsource_global[NGGMAX][NGGMAX];
   
    struct  bitriplet   Tsol[NIMAX];

    int nimage = 0;
    //struct    chaine *parent;

    //Tsol=NULL;    // head of the Tsol list.
    if ( source.dr > PREC_DLSDS )
    {
        nimage = inverse(gsource_global, &source.C, Tsol);
        NPRINTF(stderr, "COMP: %s: multiplicity: %d at (%.3lf,%.3lf)",
                source.n, nimage, source.C.x, source.C.y);

//      if ((nimage>0)&&(nimage<NIMAX))
//      {
        /*Fill the image list with the arclets of the Tsol linked list*/
        nimage = e_test(Tsol, nimage, source, image);
        NPRINTF(stderr, ": found %d image(s)\n", nimage);
        /*      }
                else
                {
                    NPRINTF(stderr,"\n");
                    // free Tsol list until it is NULL value. (usefull in case nimage > NIMAX)
                    // Go down to the last element of the Tsol list.
                    while( Tsol != NULL )
                    {
                        parent = Tsol;
                        Tsol=(Tsol->F);
                        free((struct chaine*) parent);
                    }
                }*/

    }
    else
        NPRINTF(stderr, "COMP: %s: no multiplicity because its redshift is 0\n", source.n);

    return(nimage);
}
