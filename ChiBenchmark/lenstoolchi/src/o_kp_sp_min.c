#include<stdio.h>
#include<math.h>
#include "fonction.h"
#include "constant.h"
#include"dimension.h"
#include "structure.h"

/*
*       nom:        o_kp_sp_min
*       auteur:     Jean-Paul Kneib
*       date:       september 94
*       place:      Cambridge
*/

void  o_kp_sp_min(double x1, double x2, double y1, double y2,
                  int i, int j)
{
    extern  double  x1min, x1max, y1min, y1max;
    extern  double  x2min, x2max, y2min, y2max;
    extern  int izmin, izmax;
    extern  int imapmin, imapmax;
    extern  int jmapmin, jmapmax;

    if ((y1 < y1min) || (y2 <= y1min))
    {
        imapmin = i;
        jmapmin = j;
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
        imapmax = i;
        jmapmax = j;
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
}
