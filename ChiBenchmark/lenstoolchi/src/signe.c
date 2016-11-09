#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        signe               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

/* Return >0 if P and T->c are on the same side of the [T->a,T->b]
 * segment, otherwise <0.
 *
 * Return 0, if P or c are on the [T->a,T->b] segment.
 *
 * Global variables used :
 * -  none
 */
double signe(const struct triplet *T, const struct point *P)
{
    return(determinant(&T->a, &T->b, P)*determinant(&T->a, &T->b, &T->c));
}

/*********************************************************************/
/* Return the scalar triple product (a*b).c of the 3 vectors A[x,y,1],
 * B[x,y,1], C[x,y,1].
 * If 2 of the 3 vectors are equal, colinear or form an orthogonal basis,
 * the triple product is 0.
 * This is also the determinant of the matrix
 *   | Ax  Bx  Cx |
 *   | Ay  By  Cy |
 *   |  1   1   1 |
 *
 * Global variables used :
 * - none
 */
double determinant(const struct point *A,
                   const struct point *B,
                   const struct point *C)
{
    return( B->x * C->y - B->y * C->x +
            A->x * B->y - A->y * B->x +
            A->y * C->x - A->x * C->y );
}
