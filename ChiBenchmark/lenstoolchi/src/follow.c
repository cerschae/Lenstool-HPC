#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        follow              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/
/* Follow the critical line using the trapeze algorithm. (see next function)
 * Algorithm stop when the next position is close enough to the origin O.
 *
 * Return the radial array with the sequence of image and source points
 * that are on the critical line.
 * Parameters :
 * - A and B are the initial direction of the critical line
 * - O is the origin and final position on the critical line
 *
 * Global variables used :
 * - CL, radial, nrline, flagr
 * - in next() : G, lens
 * - in e_dpl() : G, lens, lens_table
 * */
void    follow( struct point A,
                struct point B,
                struct point O,
                double dl0s, double dos, double zs )
{
    extern  struct  g_cline CL;
    extern  struct  biline  radial[NMAX];
    extern  int nrline, flagr;
    struct  point   C;
    register int    i;

    for ( i = 0; i < NRLINEMAX && nrline < 5*NRLINEMAX ; i++)
    {
        radial[nrline].i = flagr;
        C = radial[nrline].I = next(A, B, CL.cpas, dl0s, dos, zs);
        e_dpl(&C, dl0s / dos, &radial[nrline++].S);

        if ( dist(O, C) > CL.cpas )
        {
            A = B;
            B = C;
        }
        else
            return;
    };

}

/****************************************************************/
/* Return the next point of the critical line close to A and B.
 * Algorithm looks for the critical line to cross any of the 4 sides
 *  of a trapeze which symetric axis is the (AB) segment.
 * Parameters :
 * - A and B define the trapeze region of search
 * - dpl is the half size of the basis of the trapeze
 *
 * Global variables used :
 * - in e_zeroamp() : G, lens
 * - in ch_signe() : G, lens
 * */
struct  point   next(struct point A, struct point B, double dpl,
                     double dl0s, double dos, double zs)
{
    double  norme, dx, dy;
    struct  point   A1, A2, A3, A4;
    struct  point   col, perp;

    dx = B.x - A.x;
    dy = B.y - A.y;
    norme = sqrt(dx*dx + dy*dy);

    col.x = dx / norme;
    col.y = dy / norme;
    perp.x = -col.y;
    perp.y = col.x;


    A1.x = B.x + dpl*perp.x;
    A1.y = B.y + dpl*perp.y;
    A2.x = B.x + dpl*(col.x + perp.x / 2.);
    A2.y = B.y + dpl*(col.y + perp.y / 2.);
    A3.x = B.x + dpl*(col.x - perp.x / 2.);
    A3.y = B.y + dpl*(col.y - perp.y / 2.);
    A4.x = B.x - dpl*perp.x;
    A4.y = B.y - dpl*perp.y;

    if ( chsigne(A1, A2, dl0s, dos, zs) != 0 )
        return( e_zeroamp(A1, A2, dl0s, dos, zs) );
    else if ( chsigne(A2, A3, dl0s, dos, zs) != 0 )
        return( e_zeroamp(A2, A3, dl0s, dos, zs) );
    else if ( chsigne(A3, A4, dl0s, dos, zs) != 0 )
        return( e_zeroamp(A3, A4, dl0s, dos, zs) );
    else if ( dpl > .01 )
        return( next(A, B, dpl / 2., dl0s, dos, zs) );
    else
        return(A);

}
