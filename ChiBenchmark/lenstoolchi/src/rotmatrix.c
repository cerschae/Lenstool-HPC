#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

static void prdm(struct matrix *A, struct matrix *B, struct matrix *C);
static void Mrot(double theta, struct matrix *R);

/****************************************************************/
/*      nom:        rotmatrix           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Global variables used :
 * - none
 */
struct matrix rotmatrix(struct matrix *P, double theta)
/*  Matrix form | a  b |
 *              | d  c |
 *  Returns the matrix P in a new reference frame rotated by theta radians
 *  relatively to the original reference frame.
 *
 *  Result = Mrot(theta)*P*Mrot(-theta)   (see also Eq 2.45 of JPK Phd Thesis)
 *
 *  R(theta)*P*R(-theta) = [ R(theta)*P ] * R(-theta)
 */
{
    struct matrix Q, M, R;

    Mrot(-theta, &M); // Mrot(-theta)
    prdm(P, &M, &Q);
    M.b *= -1;
    M.d *= -1;  // Mrot(theta)
    prdm(&M, &Q, &R);

    return(R);
}

/****************************************************************/
/*  Matrix form | a  b |
 *              | d  c |
 * Perform the standard matrix product C = A*B
 *
 * Global variables used :
 * - none
 */
static void prdm(struct matrix *A, struct matrix *B, struct matrix *C)
{
    C->a = A->a * B->a + A->b * B->d;
    C->b = A->a * B->b + A->b * B->c;
    C->c = A->d * B->b + A->c * B->c;
    C->d = A->d * B->a + A->c * B->d;
}
/****************************************************************/
/* Return a the rotation matrix | cos(t) -sin(t) |
 *                              | sin(t)  cos(t) |
 *
 * Global variables used :
 * - none
 */
static void Mrot(double theta, struct matrix *R)
{
    R->a = cos(theta);
    R->b = -sin(theta);
    R->c = cos(theta);
    R->d = sin(theta);
}
