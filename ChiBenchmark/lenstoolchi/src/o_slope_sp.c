#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/*
*       nom:        o_slope_sp
*       auteur:     Jean-Paul Kneib
*       date:       september 94
*       place:      Cambridge
*/
const extern  struct  g_grille    G;
extern  struct  pot lens[];

int  o_slope_sp(double *y0)
{
//  extern  struct  g_mode          M;
    extern  double  **map_p;
    extern  double  **tmp_p;
//  extern  int imapmin,imapmax;
//  extern  int jmapmin,jmapmax;

    register int    i, j;

    double  x0, x1, x2, y1, y2;
    int flag = 0;
    double  dx, dy, xc, yc;
    int dirx, diry, ndx, ndy;

    dx = G.echant * G.dx;
    dy = G.echant * G.dy;
    xc = lens[0].C.x;
    yc = lens[0].C.y;
    dirx = 1;
    diry = 1;
    ndx = 1;
    ndy = 1;

    do
    {
        for (i = 0; i < ndx; i++)
        {
            if ((xc < G.xmax) && (xc > G.xmin))
            {
                lens[G.nlens-1].C.x = xc;
                tmp_p = map_p;
                x1 = G.exc;
                add_pm(map_p, G.nx, G.ny, xc, yc, x1);
                y1 = o_chi();
                x2 = -G.exc;
                add_pm(map_p, G.nx, G.ny, xc, yc, 2.*x2);
                y2 = o_chi();
                x0 = 0.;
                map_p = tmp_p;
                if ((y1 < (*y0)) && (y1 < y2))
                {
                    add_pm(map_p, G.nx, G.ny, xc, yc, x1);
                    (*y0) = y1;
                    flag++;
                }
                else if ((y2 < (*y0)) && (y2 < y1))
                {
                    add_pm(map_p, G.nx, G.ny, xc, yc, x2);
                    (*y0) = y2;
                    flag++;
                }
            }
            xc += dirx * dx;
        }
        dirx *= -1;
        ndx++;

        for (j = 0; j < ndy; j++)
        {
            if ((yc < G.ymax) && (yc > G.ymin))
            {
                lens[G.nlens-1].C.y = yc;
                tmp_p = map_p;
                x1 = G.exc;
                add_pm(map_p, G.nx, G.ny, xc, yc, x1);
                y1 = o_chi();
                x2 = -G.exc;
                add_pm(map_p, G.nx, G.ny, xc, yc, 2.*x2);
                y2 = o_chi();
                x0 = 0.;
                map_p = tmp_p;
                if ((y1 < (*y0)) && (y1 < y2))
                {
                    add_pm(map_p, G.nx, G.ny, xc, yc, x1);
                    (*y0) = y1;
                    flag++;
                }
                else if ((y2 < (*y0)) && (y2 < y1))
                {
                    add_pm(map_p, G.nx, G.ny, xc, yc, x2);
                    (*y0) = y2;
                    flag++;
                }
            }
            yc += diry * dy;
        }
        diry *= -1;
        ndy++;

    }
    while ( (xc < G.xmax) && (xc > G.xmin) && (yc < G.ymax) && (yc > G.ymin) );

    return(flag);
}

void    add_pm(double **map, int nx, int ny, double x0, double y0, double m)
{
    register int ii, jj;
    struct point    P;
    double  z;

    z = lens[G.nlens-1].rc * lens[G.nlens-1].rc;
    for (ii = 0; ii < nx; ii++)
    {
        P.x = G.xmin + ii * G.dx - x0;
        for (jj = 0; jj < ny; jj++)
        {
            P.y = G.ymin + jj * G.dy - y0;
            map[ii][jj] += m * log(1. + (P.x * P.x + P.y * P.y) / z);
        }
    }

}

/* point mass
    if ((ii==i0)&&(jj==j0))
        {
        map[ii][jj]+=m/2.*log(G.dx*G.dx/4.+G.dy*G.dy/4.);
        }
    else
        {
        y=(jj-j0)*G.dy;
        map[ii][jj]+=m/2.*log(x*x+y*y);
        }
*/
