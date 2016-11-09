#include <math.h>
#include <fonction.h>

/*
* spline
* jean-paul kneib
* september 94
* IoA cambridge
* lens tool
*/

static void sp_loc(struct point P, int *i1, int *i2, int *j1, int *j2);
static void sp_lin(int i1, int i2, int j1, int j2, struct point P, const double **map, double *ret);
static void sp_cxy(const double *xx, const double *yy, const double **map, const double **map_yy,
                   int nx, int ny, int i1, int i2, int j1, int j2,
                   double x, double y, double *dzdx);
static void spint(const double *vx, const double *za, const double *z2a, int k1, int k2, double x, double *z);
static void spint_x(const double *vx, const double *za, const double *z2a, int k1, int k2, double x, double *dzdx);
static void sp_comp(const double *xx, const double *yy, const double **map, const double **map_yy,
                    int nx, int ny, int i1, int i2, int j1, int j2,
                    double x, double y, double *zz);
static void sp_cx(const double *xx, const double *yy, const double **map, const double **map_yy,
                  int nx, int ny, int i1, int i2, int j1, int j2,
                  double x, double y, double *dzdx);
static void sp_cy(const double *xx, const double *yy, const double **map, const double **map_yy,
                  int nx, int ny, int i1, int i2, int j1, int j2,
                  double x, double y, double *dzdy);


const extern  struct  g_grille    G;
const extern  double  **map_p;
const extern  double  **map_axx;
const extern  double  **map_axy;
const extern  double  **map_ayy;
const extern  double  **map_x4;
const extern  double  **map_y4;
const extern  double  *v_xx;
const extern  double  *v_yy;

/*
* spline interpolation of the potential
*/

double  sp_pot(struct point P)
{
    double  pot = 0.;
    int i1, i2, j1, j2;

    /* localisation of the current point in the map */

    sp_loc(P, &i1, &i2, &j1, &j2);

    if ((i1 > 0) && (i2 > 1) && (i1 < G.nx - 2) && (i2 < G.nx - 1)
        && (j1 > 0) && (j2 > 1) && (j1 < G.ny - 2) && (j2 < G.ny - 1) )
    {
        sp_comp(v_xx, v_yy, map_p, map_ayy, G.nx, G.ny, i1, i2, j1, j2, P.x, P.y, &pot);
    }
    return(pot);
}

/*
* spline interpolation of the potential first derivates
*
* Global variables used :
* - none
*/

struct point    sp_grad(struct point P)
{
    struct point    dpl;
    int i1, i2, j1, j2;

    dpl.x = dpl.y = 0.;

    sp_loc(P, &i1, &i2, &j1, &j2);

    if ((i1 > 0) && (i2 > 1) && (i1 < G.nx - 2) && (i2 < G.nx - 1)
        && (j1 > 0) && (j2 > 1) && (j1 < G.ny - 2) && (j2 < G.ny - 1) )
    {
        sp_cx(v_xx, v_yy, map_p, map_ayy, G.nx, G.ny, i1, i2, j1, j2, P.x, P.y, &dpl.x);
        sp_cy(v_xx, v_yy, map_p, map_ayy, G.nx, G.ny, i1, i2, j1, j2, P.x, P.y, &dpl.y);
    }
    return(dpl);
}

/*
* spline interpolation of the potential
* Global variables used :
* - none
*/

struct matrix sp_grad2(struct point P)
{
    struct matrix   mag;
    int i1, i2, j1, j2;

    mag.a = mag.b = mag.c = mag.d = 0.;

    sp_loc(P, &i1, &i2, &j1, &j2);

    if ((i1 > 0) && (i2 > 1) && (i1 < G.nx - 2) && (i2 < G.nx - 1)
        && (j1 > 0) && (j2 > 1) && (j1 < G.ny - 2) && (j2 < G.ny - 1) )
    {
        sp_lin(i1, i2, j1, j2, P, map_axx, &mag.a);
        sp_lin(i1, i2, j1, j2, P, map_ayy, &mag.c);
        /*  sp_comp(v_xx,v_yy,map_axx,map_x4,G.nx,G.ny,i1,i2,j1,j2,P.x,P.y,&mag.a);
            sp_comp(v_xx,v_yy,map_ayy,map_y4,G.nx,G.ny,i1,i2,j1,j2,P.x,P.y,&mag.a);*/
        sp_cxy(v_xx, v_yy, map_p, map_ayy, G.nx, G.ny, i1, i2, j1, j2, P.x, P.y, &mag.b);
        mag.d = mag.b;
    }

    return(mag);
}

/* Global variables used :
 * - none
 */
static void sp_loc(struct point P, int *i1, int *i2, int *j1, int *j2)
{
    double q;

    q = (P.x - G.xmin) / G.dx;
    *i1 = (int)(floor(q));
    /* *i2=(int) (ceil(q));*/
    *i2 = (*i1) + 1;
    q = (P.y - G.ymin) / G.dy;
    *j1 = (int) floor(q);
    *j2 = (*j1) + 1;
}

/* Global variables used :
 * - none
 */
static void sp_lin(int i1, int i2, int j1, int j2, struct point P, const double **map, double *ret)
{
    double  a, b;
    double  x1, y1;

    x1 = G.xmin + i1 * G.dx;
    y1 = G.ymin + j1 * G.dy;

    a = map[i1][j1] + (P.x - x1) * (map[i2][j1] - map[i1][j1]) / G.dx;
    b = map[i1][j2] + (P.x - x1) * (map[i2][j2] - map[i1][j2]) / G.dx;

    *ret = a + (P.y - y1) * (b - a) / G.dy;
}

static void sp_cxy(const double *xx, const double *yy, const double **map, const double **map_yy,
                   int nx, int ny, int i1, int i2, int j1, int j2,
                   double x, double y, double *dzdx)
{
    register int i;
    double mp2[100], mp[100];
    /*  double *mp2,*mp;

        mp2=(double *) alloc_vector_double(ny);
        mp=(double *) alloc_vector_double(ny);*/

    for (i = 0; i < nx; i++)
        spint_x(yy, map[i], map_yy[i], j1, j2, y, &mp[i]);

    spline(xx, mp, nx, 1.0e30, 1.0e30, mp2);
    spint_x(xx, mp, mp2, i1, i2, x, dzdx);

    /*
        free((double *) mp);
        free((double *) mp2);
    */

}

static void spint(const double *vx, const double *za, const double *z2a, int k1, int k2, double x, double *z)
{
    double h, b, a;

    h = vx[k2] - vx[k1];
    if (h == 0.0) fprintf(stderr, "Bad vx input to routine splint\n");
    a = (vx[k2] - x) / h;
    b = (x - vx[k1]) / h;
    *z = a * za[k1] + b * za[k2]
         + (h * h) / 6.*((a * a - 1.) * a * z2a[k1] + (b * b - 1.) * b * z2a[k2]);

}

static void spint_x(const double *vx, const double *za, const double *z2a, int k1, int k2, double x, double *dzdx)
{
    double h, b, a;

    h = vx[k2] - vx[k1];
    if (h == 0.0) fprintf(stderr, "Bad vx input to routine splint\n");
    a = (vx[k2] - x) / h;
    b = (x - vx[k1]) / h;
    *dzdx = (za[k2] - za[k1]) / h
            - h / 6.*((3.*a * a - 1.) * z2a[k1] - (3.*b * b - 1.) * z2a[k2]);

}

static void sp_comp(const double *xx, const double *yy, const double **map, const double **map_yy,
                    int nx, int ny, int i1, int i2, int j1, int j2,
                    double x, double y, double *zz)
{
    register int i;
    double mp2[200], mp[200];

    /*mp2=(double *) alloc_vector_double(ny);
    mp=(double *) alloc_vector_double(ny);*/

    for (i = 0; i < nx; i++)
        spint(yy, map[i], map_yy[i], j1, j2, y, &mp[i]);
    spline(xx, mp, nx, 1.0e30, 1.0e30, mp2);
    spint(xx, mp, mp2, i1, i2, x, zz);

    /*
        free((double *) mp);
        free((double *) ytmp);
    */

}

static void sp_cx(const double *xx, const double *yy, const double **map, const double **map_yy,
                  int nx, int ny, int i1, int i2, int j1, int j2,
                  double x, double y, double *dzdx)
{
    register int i;
    double mp2[200], mp[200];
    /*  double  *mp2,*mp;

        mp2=(double *) alloc_vector_double(ny);
        mp=(double *) alloc_vector_double(ny);*/

    for (i = 0; i < nx; i++)
        spint(yy, map[i], map_yy[i], j1, j2, y, &mp[i]);
    spline(xx, mp, nx, 1.0e30, 1.0e30, mp2);
    spint_x(xx, mp, mp2, i1, i2, x, dzdx);


    /*
        free((double *) mp);
        free((double *) mp2);
    */
}

static void sp_cy(const double *xx, const  double *yy, const double **map, const double **map_yy,
                  int nx, int ny, int i1, int i2, int j1, int j2,
                  double x, double y, double *dzdy)
{
    register int i;
    double mp2[200], mp[200];
    /*  double *mp2,*mp;

        mp2=(double *) alloc_vector_double(ny);
        mp=(double *) alloc_vector_double(ny);*/

    for (i = 0; i < nx; i++)
        spint_x(yy, map[i], map_yy[i], j1, j2, y, &mp[i]);
    spline(xx, mp, nx, 1.0e30, 1.0e30, mp2);
    spint(xx, mp, mp2, i1, i2, x, dzdy);


    /*
        free((double *) mp);
        free((double *) mp2);
    */

}
