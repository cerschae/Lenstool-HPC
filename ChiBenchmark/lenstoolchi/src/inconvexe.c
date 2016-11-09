#include <grille.h>

/*
 * Parameters :
 * - P point that is inside or outside the contour
 * - np number of points of the contour
 * - I[NPOINTs] polygonal contour
 *
 * Return 1 if P is inside the contour, otherwise 0
 */

int inconvexe(struct point P, int np, struct point I[NPOINT])
{
    register int    k;
    double  xi, xj, yi, yj, theta = 0.;

    I[np].x = I[0].x;
    I[np].y = I[0].y;
    for (k = 0; k < np; k++)
    {
        xi = I[k].x - P.x;
        xj = I[k+1].x - P.x;
        yi = I[k].y - P.y;
        yj = I[k+1].y - P.y;
        /* in a triangle OIJ, (OI,OJ)=atan(Det(OI,OJ)/Pscal(OI,OJ) */
        theta += atan2((xj * yi - xi * yj), (xi * xj + yi * yj));
    };

    if (fabs(theta) > PI) /* if P is on the edge of the polygon theta=PI */
        return(1);
    else    /*theta is close to 0*/
        return(0);
}
