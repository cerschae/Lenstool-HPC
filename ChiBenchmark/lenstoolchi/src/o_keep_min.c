#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        o_keep_min          */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            *
 ***************************************************************
 * Keep in
 * -    ipmin, ipmax the lens potential index
 * -    ilsmin, ilsmax the parameter index
 * that has generated the min and max Khi2 values
 *
 * If y1 or y2 is smaller the min Khi2 (y1min) value then
 *  Keep in
 *      - x1min, y1min the minimum coordinates of this min Khi2
 *      - x2min, y2min the other coordinates
 * Else if
 *  Keep in
 *      - x1max, y1max the maximum coordinates of this max Khi2
 *      - x2min, y2max the other coordinates
 *
 */
void  o_keep_min(double x1, double x2, double y1, double y2, int ils, int ipx)
{
    extern  double  x1min, x1max, y1min, y1max;
    extern  double  x2min, x2max, y2min, y2max;
    extern  int izmin, izmax;
    extern  int ipmin, ipmax;
    extern  int ilsmin, ilsmax;
//  extern  double  excd[][NPAMAX];
//  extern  double  excu[][NPAMAX];

    if ((y1 < y1min) || (y2 <= y1min))
    {
        ipmin = ipx;
        ilsmin = ils;
        izmin = -1;
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
    else if ((y1 >= y1max) || (y2 > y1max))
    {
        ipmax = ipx;
        ilsmax = ils;
        izmax = -1;
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
        {
        excd[ils][ipx]=Min( excd[ils][ipx]*1.1,.5);
        excu[ils][ipx]=Min( excd[ils][ipx]*1.1,.5);
        };
    */
}
