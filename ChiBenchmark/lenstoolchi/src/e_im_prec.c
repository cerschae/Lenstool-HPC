#include<stdio.h>
#include<math.h>
#include "fonction.h"
#include "constant.h"
#include"dimension.h"
#include "structure.h"

//#define DDEBUG
/****************************************************************/
/*      nom:        e_im_prec               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse
 *
 * Return in *A the smallest couple of triangles (source, image plane) in which P is
 * located according to the minimum surface required for the source triangle.
 *
 * Stop when the triangle in the source plane contains the source P"and
 * has a surface lower than DMIN^2.
 *
 * it variable counts the number of loop needed to reach the smallest triangle.
 *
 * Warning : The function can return a triangle that doesnt comply
 * with those requirements if the it variable becomes larger
 * than NITMAX.
 *
 * Parameters :
 * - E : Triplet of simulated images and corresponding sources (6 points)
 * - P : Barycenter of the original sources
 * - dlsds : Lens efficiency
 * - it : loop counter
 * - A : returned bitriangle
 *
 * Global variables used :
 * - in e_inthere() : G, lens, lens_table
 * - in e_transform() : G, lens, lens_table
 */
void e_im_prec(const struct bitriplet *E, /* Triplet of simulated images and corresponding sources*/
               const struct point *P,     /* Barycenter of the original sources */
               double dlsds, int *it,
               struct bitriplet *A)

{
    struct bitriplet I;
    double det;

#ifdef DDEBUG
    // plot the E triangle in red in the image plane
    // and in blue in the source plane
    const extern struct g_mode M;
    FILE  *dbg;

    dbg = fopen( "imsearch.reg", "a");
    fprintf( dbg, "fk5;polygon(%lf,%lf,%lf,%lf,%lf,%lf) #color=red\n",
             M.ref_ra + E->i.a.x / -3600 / cos(M.ref_dec*DTR),
             M.ref_dec + E->i.a.y / 3600,
             M.ref_ra + E->i.b.x / -3600 / cos(M.ref_dec*DTR),
             M.ref_dec + E->i.b.y / 3600,
             M.ref_ra + E->i.c.x / -3600 / cos(M.ref_dec*DTR),
             M.ref_dec + E->i.c.y / 3600 );
    fprintf( dbg, "fk5;polygon(%lf,%lf,%lf,%lf,%lf,%lf) #color=blue\n",
             M.ref_ra + E->s.a.x / -3600 / cos(M.ref_dec*DTR),
             M.ref_dec + E->s.a.y / 3600,
             M.ref_ra + E->s.b.x / -3600 / cos(M.ref_dec*DTR),
             M.ref_dec + E->s.b.y / 3600,
             M.ref_ra + E->s.c.x / -3600 / cos(M.ref_dec*DTR),
             M.ref_dec + E->s.c.y / 3600 );
    fclose(dbg);
#endif

    interieur(&E->i, &I.i);       //return the triangle I.i inside the triangle E.i
    e_transform(&I.i, dlsds, &I.s); //return the triangle I.i in the source plane
    e_inthere(E, &I, P, dlsds, A);  //return the small couple of triangles A in which P is located

    det = determinant(&A->s.a, &A->s.b, &A->s.c); //det is the surface of the source triangle
    if ( fabs(det) < DMIN*DMIN || (*it)++ > NITMAX ) // cf also e_test.c::
        return;
    else
    {
        I.i = A->i;
        I.s = A->s;
        e_im_prec(&I, P, dlsds, it, A);
    }
}
