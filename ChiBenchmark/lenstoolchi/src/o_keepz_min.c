#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        o_keepz_min         */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * In the zm_limit optimisation context :
 * ------------------------------------
 * - Modify the global variables x1min,x2min,y1min and y2min to keep
 * in memory the min and max Khi2 positions.
 * - In case of modification of any of those variables, assign to
 * the global variable izmin, the id of the zm_limit source.
 */
void  o_keepz_min(double x1, double x2, double y1, double y2, int iz)
{
    extern  double  x1min, x1max, y1min, y1max;
    extern  double  x2min, x2max, y2min, y2max;
    extern  int izmin, izmax;
    extern  int ipmin, ipmax;
    extern  int ilsmin, ilsmax;
//  extern  struct  z_lim   zlim[];

    if ((y1 < y1min) || (y2 < y1min))
    {
        ipmin = -1;
        ilsmin = -1;
        izmin = iz;
        if (y1 < y2)
        {
            x1min = x1;
            y1min = y1;
            x2min = x2;
            y2min = y2;
        }
        else
        {
            x1min = x2;
            y1min = y2;
            x2min = x1;
            y2min = y1;
        };
    }
    else if ((y1 > y1max) || (y2 > y1max))
    {
        ipmax = -1;
        ilsmax = -1;
        izmax = iz;
        if (y1 < y2)
        {
            x1max = x2;
            y1max = y2;
            x2max = x1;
            y2max = y1;
        }
        else
        {
            x1max = x1;
            y1max = y1;
            x2max = x2;
            y2max = y2;
        };
    }
    /*
    else
        { zlim[izmax].excu=.2; zlim[izmax].excd=.2; };
    */
}
