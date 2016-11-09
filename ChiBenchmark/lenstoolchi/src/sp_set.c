#include <stdlib.h>
#include "lt.h"

static void sp_ixx(const double *x, double **map, int nx, int ny, double p1, double pn, double **map_xx);
static void sp_iyy(const double *y, double **map, int nx, int ny, double p1, double pn, double **map_yy);


/*Global variables used :
 * - v_xx, v_yy
 * */
void sp_set(double **map, int nx, int ny, double **map_xx, double **map_yy)
{
    const extern double *v_xx;
    const extern double *v_yy;

    sp_ixx(v_yy, map, nx, ny, 1.e30, 1.e30, map_xx);
    sp_iyy(v_xx, map, nx, ny, 1.e30, 1.e30, map_yy);

}

//NOT USED IN LENSTOOL
void sp_ix( double *x, double **map, int nx, int ny, double p1, double pn,
            double **map_x, double **map_xx)
{
    register int i, j, k;
    double p, qn, sig, un, *u;
    double deltax, deltam;
    double dx, d2x;

    u = (double *) alloc_vector_double(nx);

    for (j = 0; j < ny; j++)
    {
        if (p1 > 0.99e30)
            map_xx[0][j] = u[0] = 0.0;
        else
        {
            map_xx[0][j] = -0.5;
            u[0] = (3.0 / (x[1] - x[0])) * ((map[1][j] - map[0][j]) / (x[1] - x[0]) - p1);
        }
        for (i = 1; i < nx - 1; i++)
        {
            dx = x[i] - x[i-1];
            d2x = x[i+1] - x[i-1];
            sig = dx / d2x;
            p = sig * map_xx[i-1][j] + 2.0;
            map_xx[i][j] = (sig - 1.0) / p;
            u[i] = (map[i+1][j] - map[i][j]) / (x[i+1] - x[i])
                   - (map[i][j] - map[i-1][j]) / dx;
            u[i] = (6.0 * u[i] / d2x - sig * u[i-1]) / p;
        }
        if (pn > 0.99e30)
            qn = un = 0.0;
        else
        {
            qn = 0.5;
            un = (3.0 / (x[nx-1] - x[nx-2])) * (pn - (map[nx-1][j] - map[nx-2][j]) / (x[nx-1] - x[nx-2]));
        }
        map_xx[nx-1][j] = (un - qn * u[nx-2]) / (qn * map_xx[nx-2][j] + 1.0);
        for (k = nx - 2; k >= 0; k--)
            map_xx[k][j] = map_xx[k][j] * map_xx[k+1][j] + u[k];

        /* first derivatives */

        deltax = deltam = 0;
        for (k = 0; k < nx - 1; k++)
        {
            deltax = x[k+1] - x[k];
            deltam = map[k+1][j] - map[k][j];
            map_x[k][j] = deltam / deltax - deltax / 3.*(map_xx[k][j] + map_xx[k+1][j] / 2.);
        }

        map_x[nx-1][j] = deltam / deltax + deltax / 3.*(map_xx[nx-2][j] / 2. + map_xx[nx-1][j]);


    }

    free((double *) u);
}

//NOT USED IN LENSTOOL
void sp_iy( double *y, double **map, int nx, int ny, double p1, double pn,
            double **map_y, double **map_yy)
{

    register int i, j, k;
    double p, qn, sig, un, *u;
    double deltay, deltam;
    double dy, d2y;

    u = (double *) alloc_vector_double(ny);

    for (i = 0; i < nx; i++)
    {
        if (p1 > 0.99e30)
            map_yy[i][0] = u[0] = 0.0;
        else
        {
            map_yy[i][0] = -0.5;
            u[0] = (3.0 / (y[1] - y[0])) * ((map[i][1] - map[i][0]) / (y[1] - y[0]) - p1);
        }
        for (j = 1; j < ny - 1; j++)
        {
            dy = y[j] - y[j-1];
            d2y = y[j+1] - y[j-1];
            sig = dy / d2y;
            p = sig * map_yy[i][j-1] + 2.0;
            map_yy[i][j] = (sig - 1.0) / p;
            u[j] = (map[i][j+1] - map[i][j]) / (y[j+1] - y[j])
                   - (map[i][j] - map[i][j-1]) / dy;
            u[j] = (6.0 * u[j] / d2y - sig * u[j-1]) / p;
        }
        if (pn > 0.99e30)
            qn = un = 0.0;
        else
        {
            qn = 0.5;
            un = (3.0 / (y[ny-1] - y[ny-2])) * (pn - (map[i][ny-1] - map[i][ny-2]) / (y[ny-1] - y[ny-2]));
        }
        map_yy[i][ny-1] = (un - qn * u[ny-2]) / (qn * map_yy[i][ny-2] + 1.0);

        for (k = ny - 2; k >= 0; k--)
            map_yy[i][k] = map_yy[i][k] * map_yy[i][k+1] + u[k];

        /* first derivatives */

        deltam = deltay = 0;
        for (k = 0; k < ny - 1; k++)
        {
            deltay = y[k+1] - y[k];
            deltam = map[i][k+1] - map[i][k];
            map_y[i][k] = deltam / deltay - deltay / 3.*(map_yy[i][k] + map_yy[i][k+1] / 2.);
        }

        map_y[i][ny-1] = deltam / deltay + deltay / 3.*(map_yy[i][ny-2] / 2. + map_yy[i][ny-1]);

    }

    free((double *) u);
}

//NOT USED IN LENSTOOL
void sp_ixx(const double *x, double **map, int nx, int ny, double p1, double pn, double **map_xx)
{

    register int i, j, k;
    double p, qn, sig, un, *u;

    u = (double *) alloc_vector_double(nx);

    for (j = 0; j < ny; j++)
    {
        if (p1 > 0.99e30)
            map_xx[0][j] = u[0] = 0.0;
        else
        {
            map_xx[0][j] = -0.5;
            u[0] = (3.0 / (x[1] - x[0])) * ((map[1][j] - map[0][j]) / (x[1] - x[0]) - p1);
        }
        for (i = 1; i < nx - 1; i++)
        {
            sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
            p = sig * map_xx[i-1][j] + 2.0;
            map_xx[i][j] = (sig - 1.0) / p;
            u[i] = (map[i+1][j] - map[i][j]) / (x[i+1] - x[i]) - (map[i][j] - map[i-1][j]) / (x[i] - x[i-1]);
            u[i] = (6.0 * u[i] / (x[i+1] - x[i-1]) - sig * u[i-1]) / p;
        }
        if (pn > 0.99e30)
            qn = un = 0.0;
        else
        {
            qn = 0.5;
            un = (3.0 / (x[nx-1] - x[nx-2])) * (pn - (map[nx-1][j] - map[nx-2][j]) / (x[nx-1] - x[nx-2]));
        }
        map_xx[nx-1][j] = (un - qn * u[nx-2]) / (qn * map_xx[nx-2][j] + 1.0);
        for (k = nx - 2; k >= 0; k--)
            map_xx[k][j] = map_xx[k][j] * map_xx[k+1][j] + u[k];
    }

    free((double *)u);
}

//NOT USED IN LENSTOOL
void sp_iyy(const double *y, double **map, int nx, int ny, double p1, double pn, double **map_yy)
{

    register int i, j, k;
    double p, qn, sig, un, *u;

    u = (double *) alloc_vector_double(ny);

    for (i = 0; i < nx; i++)
    {
        if (p1 > 0.99e30)
            map_yy[i][0] = u[0] = 0.0;
        else
        {
            map_yy[i][0] = -0.5;
            u[0] = (3.0 / (y[1] - y[0])) * ((map[i][1] - map[i][0]) / (y[1] - y[0]) - p1);
        }
        for (j = 1; j < ny - 1; j++)
        {
            sig = (y[j] - y[j-1]) / (y[j+1] - y[j-1]);
            p = sig * map_yy[i][j-1] + 2.0;
            map_yy[i][j] = (sig - 1.0) / p;
            u[j] = (map[i][j+1] - map[i][j]) / (y[j+1] - y[j]) - (map[i][j] - map[i][j-1]) / (y[j] - y[j-1]);
            u[j] = (6.0 * u[j] / (y[j+1] - y[j-1]) - sig * u[j-1]) / p;
        }
        if (pn > 0.99e30)
            qn = un = 0.0;
        else
        {
            qn = 0.5;
            un = (3.0 / (y[ny-1] - y[ny-2])) * (pn - (map[i][ny-1] - map[i][ny-2]) / (y[ny-1] - y[ny-2]));
        }
        map_yy[i][ny-1] = (un - qn * u[ny-2]) / (qn * map_yy[i][ny-2] + 1.0);
        for (k = ny - 2; k >= 0; k--)
            map_yy[i][k] = map_yy[i][k] * map_yy[i][k+1] + u[k];
    }

    free((double *) u);
}
