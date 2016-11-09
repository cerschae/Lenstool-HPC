#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        polxy               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Global variables used :
 * - none
 */

struct  polar   polxy(struct point xy)
{
    struct  polar   pol;
    pol.r = sqrt(xy.x*xy.x + xy.y*xy.y);
    pol.theta = atan2(xy.y, xy.x);
    return(pol);
}
