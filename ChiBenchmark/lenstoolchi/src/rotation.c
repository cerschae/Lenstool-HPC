#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        rotation            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/
/* Global variables used :
 * - none
 * */
struct  point   rotation(struct point P, double theta)
{
    struct  point   Q;

    Q.x = P.x*cos(theta) + P.y*sin(theta);
    Q.y = P.y*cos(theta) - P.x*sin(theta);

    return(Q);
}
