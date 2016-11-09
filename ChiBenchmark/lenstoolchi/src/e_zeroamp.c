#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

static struct   point   e_zerodich2(struct point *A, struct point *B, double dl0s, double dos, double zs);

/****************************************************************/
/*      nom:        e_zeroamp           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/
/* Return the critical line position between A and B
 * Global variables used :
 * in e_amp() : G, lens, lens_table
 * in e_zerodich2() : G, lens, lens_table
 */
struct  point   e_zeroamp(struct point A, struct point B, double dl0s, double dos, double zs)
{
    struct point C;

    /*order A and B so that amp(A)>amp(B)*/
    if ( e_amp(&A, dl0s, dos, zs) > e_amp(&B, dl0s, dos, zs))
    {
        C = A;
        A = B;
        B = C;
    };

    return( e_zerodich2(&A, &B, dl0s, dos, zs) );
}

/****************************************************************/
/* Return the critical line position between A and B
 * Global variables used :
 * - in e_amp() : G, lens, lens_table
 */
static struct   point   e_zerodich2(struct point *A, struct point *B, double dl0s, double dos, double zs)
{
    double  am;
    struct  point   M;

    M = milieu(A, B);
    am = e_amp(&M, dl0s, dos, zs);
    if ( dist(*A, *B) < PREC_ZERO || fabs(am) < PREC_ZERO )
        return(M);
    else if ( am * e_amp(A, dl0s, dos, zs) > 0 )
        /*M and A are on the same side of the critical line (same parity)*/
        return(e_zerodich2(&M, B, dl0s, dos, zs));
    else
        return(e_zerodich2(A, &M, dl0s, dos, zs));
}
