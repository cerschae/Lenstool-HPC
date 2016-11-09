#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        dist                */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

/* Euclidean distance between 2 points
 * Global variables used :
 * - none
 */
double dist(struct point A, struct point B)
{
    double  x, y;
    x = A.x - B.x;
    y = A.y - B.y;
    return(sqrt(x*x + y*y));
}

/**************************************************************
 * Global variables used :
 * - none
 */
double  dist2(struct point A, struct point B)
{
    double  x, y;
    x = A.x - B.x;
    y = A.y - B.y;
    return(x*x + y*y);
}
/*
*   Barycentre of a triplet/triangle -----------------------------
* A is a structure triplet that contains 3 structures point a,b and c
* Return value B is a point
*
* Global variables used :
* - none
*/
struct  point   barycentre(struct triplet *A)
{
    struct  point   B;

    B.x = (A->a.x + A->b.x + A->c.x) / 3.;
    B.y = (A->a.y + A->b.y + A->c.y) / 3.;
    return(B);
}

/*
 *  center of two points ----------------------------------------*/
struct  point   milieu(struct point *A, struct point *B)
{
    struct  point   C;

    C.x = (A->x + B->x) / 2.;
    C.y = (A->y + B->y) / 2.;
    return(C);
}

/*
* weighted center
* Global variables used :
* - none
*/
struct  point   wcenter(struct point A, double wa, struct point B, double wb)
{
    struct  point   C;
    double  w;
    w = wa + wb;
    C.x = (wa*A.x + wb*B.x) / w;
    C.y = (wa*A.y + wb*B.y) / w;
    return(C);
}

/*
* barycenter of a list of points P
* n is the number of values in P
*
* Global variables used :
* - none
*/
struct  point   bcentlist(struct point *P, int n)
{
    struct  point   B;
    register int i;

    B.x = B.y = 0.;
    for (i = 0; i < n; i++)
    {
        B.x += P[i].x;
        B.y += P[i].y;
    }
    B.x /= n;
    B.y /= n;
    return(B);
}

