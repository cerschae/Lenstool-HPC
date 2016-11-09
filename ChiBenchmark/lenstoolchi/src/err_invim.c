#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

static double   incaustic(struct point P);

/****************************************************************/
/*                 Program  : err_invim             */
/*                 Version  : 1 mars 1994           */
/*                 Location : ESO La Silla          */
/*                 Auteur   : jean-paul             */
/*                 But   : calcul de l'erreur d'inversion   */
/****************************************************************
 * Global variables used :
 * - ps
 * */
double  err_invim(double **errim, int **imult)
{
    const extern struct   g_pixel ps;

    register    int ii, jj;
    double  norm = 1., err = 0., inc;

    struct  point   P;

    for (ii = 0; ii < ps.ny; ii++)
    {
        P.y = ps.ymin + ps.pixely * ii;

        for (jj = 0; jj < ps.nx; jj++)
        {
            P.x = ps.xmin + ps.pixelx * jj;
            if ((errim[ii][jj] > 0.) || (imult[ii][jj] > 1))
            {
                inc = incaustic(P);
                err += inc * errim[ii][jj];
                norm += inc;
            };
        };
    };


    return(err / norm);
}

/*************************************************************/
/* Return
 * - 1 if P is inside a radial caustic line,
 * - 2 if P is also inside a tangeantial caustic line
 * - 0 otherwise
 * Modified : EJ 27/12/05, replace nrline by ntline in tangent
 *
 * Global variables used :
 * - biline, radial, tangent, nrline, ntline
 * */
static double   incaustic(struct point P)
{
    extern struct biline radial[], tangent[];
    extern int nrline, ntline;
    register int    k;
    double  xi, xj, yi, yj, theta = 0.;
    double  in = 0.;

    radial[nrline].S.x = radial[0].S.x;
    radial[nrline].S.y = radial[0].S.y;
    for (k = 0; k < nrline; k++)
    {
        xi = radial[k].S.x - P.x;
        xj = radial[k+1].S.x - P.x;
        yi = radial[k].S.y - P.y;
        yj = radial[k+1].S.y - P.y;
        theta += atan2((xj * yi - xi * yj), (xi * xj + yi * yj));
    };

    if (fabs(theta) > PI)
        in += 1.;

    tangent[ntline].S.x = tangent[0].S.x;
    tangent[ntline].S.y = tangent[0].S.y;
    for (k = 0; k < ntline; k++)
    {
        xi = tangent[k].S.x - P.x;
        xj = tangent[k+1].S.x - P.x;
        yi = tangent[k].S.y - P.y;
        yj = tangent[k+1].S.y - P.y;
        theta += atan2((xj * yi - xi * yj), (xi * xj + yi * yj));
    };

    if (fabs(theta) > PI)
        in += 1.;

    return(in);
}
