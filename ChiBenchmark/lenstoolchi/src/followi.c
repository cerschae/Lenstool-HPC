#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        followi             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * See follow()
 *
 * Global variables used :
 * - CL, ntline, flagt, tangent
 * - in next() : G, lens
 * - in e_dpl() : G, lens, lens_table
 */
void followi(struct point A,
             struct point B,
             struct point O,
             double dl0s, double dos, double zs)
{
    extern  struct  g_cline CL;
    extern  int ntline, flagt;
    extern  struct  biline  tangent[];
    struct  point   C;
    register int    i;

    for ( i = 0 ; i < NTLINEMAX && ntline < 5*NTLINEMAX ; i++ )
    {
        tangent[ntline].i = flagt;
        C = tangent[ntline].I = next(A, B, .5 * CL.cpas, dl0s, dos, zs);
        e_dpl(&C, dl0s / dos, &tangent[ntline++].S);

        if ( dist(O, C) > CL.cpas )
        {
            A = B;
            B = C;
        }
        else
            return;
    }
}
