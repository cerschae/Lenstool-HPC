#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        chsigne             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/*****************************************************************
 *  Return 0 if A and B are on the same side of a critical line
 * If A and B are on each side of a critical line then
 *      if A or B have 2 positive partial parities then Return 1 (radial CL)
 *      else Return 2  (tangential CL)
 *
 * Global variables used :
 * -  in e_unmag() : G, lens
 */
int chsigne( struct point A, struct point B, double dl0s, double dos, double zs )
{
    struct ellipse ampa, ampb;
    double sa, sb;  /*total parities of A and B*/

    ampa = e_unmag(&A, dl0s, dos, zs);
    ampb = e_unmag(&B, dl0s, dos, zs);
    sa = ampa.a * ampa.b;
    sb = ampb.a * ampb.b;

    if (sa*sb >= 0.)
        return(0);
    else
    {
        if ( ((sa > 0.) && (ampa.a > 0.)) || ((sb > 0.) && (ampb.a > 0.)) )
            return(1);
        else
            return(2);
    };
}
