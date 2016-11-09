#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        inside              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Return 1 if P is inside the triangle T, 0 otherwise.
 * Parameters :
 * - P : a point
 * - T : a triplet of points.
 */
int inside(struct point *P, struct triplet *T)
{
    double  s, s1, s2, d;

    d = determinant(&T->a, &T->b, &T->c);
    s = determinant(&T->a, &T->b, P) * d;
    s1 = determinant(&T->b, &T->c, P) * d;
    s2 = determinant(&T->c, &T->a, P) * d;

    return((s > 0.) && (s1 > 0.) && (s2 > 0.));
}

/****************************************************************
 * Return 1 if P is inside the triangle T or on its border, 0 otherwise.
 * Parameters :
 * - P : a point
 * - T : a triplet of points.
 */
int insidebord(struct point *P, struct triplet *T)
{
    double  s, s1, s2, d;

    d = determinant(&T->a, &T->b, &T->c);
    s = determinant(&T->a, &T->b, P) * d;
    s1 = determinant(&T->b, &T->c, P) * d;
    s2 = determinant(&T->c, &T->a, P) * d;

    return((s >= 0.) && (s1 >= 0.) && (s2 >= 0.));
}
