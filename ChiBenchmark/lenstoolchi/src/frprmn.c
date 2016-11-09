#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fonction.h>

#define ITMAX 200
#define EPS 1.0e-10
#define TOL 2.0e-4
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);


static double f1dim(double x);
static void linmin(double *p, double *xi, int n, double *fret, double (*func)(double*));
static void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
                   double *fc, double (*func)(double));
static double (*nrfunc)(double*);

int ncom;   /* defined in LINMIN */
double *pcom, *xicom;

void frprmn(double *p, int n, double ftol, int *iter, double *fret,
            double (*func)(double*), void (*dfunc)(double*, double*) )
{
    int j, its;
    double gg, gam, fp, dgg;
    double *g, *h, *xi;

    g = (double *) malloc(n * sizeof(double));
    h = (double *) malloc(n * sizeof(double));
    xi = (double *) malloc(n * sizeof(double));
    fp = (*func)(p);
    (*dfunc)(p, xi);
    for (j = 0; j < n; j++)
    {
        g[j] = -xi[j];
        xi[j] = h[j] = g[j];
    }
    for (its = 0; its < ITMAX; its++)
    {
        *iter = its;
        fprintf(stderr, "%d %.2lf %.2lf\n", its, fp, *fret);
        linmin(p, xi, n, fret, func);
        if (2.0*fabs(*fret - fp) <= ftol*(fabs(*fret) + fabs(fp) + EPS))
        {
            free((double *) xi);
            free((double *) h);
            free((double *) g);
            return;
        }
        fp = (*func)(p);
        (*dfunc)(p, xi);
        dgg = gg = 0.0;
        for (j = 0; j < n; j++)
        {
            gg += g[j] * g[j];
            /*        dgg += xi[j]*xi[j];   */
            dgg += (xi[j] + g[j]) * xi[j];
        }
        if (gg == 0.0)
        {
            free((double *) xi);
            free((double *) h);
            free((double *) g);
            return;
        }
        gam = dgg / gg;
        for (j = 0; j < n; j++)
        {
            g[j] = -xi[j];
            xi[j] = h[j] = g[j] + gam * h[j];
        }
    }
    fprintf(stderr, "Too many iterations in FRPRMN\n");
}

static void linmin(double *p, double *xi, int n, double *fret, double (*func)(double*))
{
    ncom = 0; /* defining declarations */
    pcom = NULL;
    xicom = NULL;

    int j;
    double xx, xmin, fx, fb, fa, bx, ax;

    ncom = n;
    pcom = (double *) malloc(n * sizeof(double));
    xicom = (double *) malloc(n * sizeof(double));
    nrfunc = func;
    for (j = 0; j < n; j++)
    {
        pcom[j] = p[j];
        xicom[j] = xi[j];
    }
    ax = 0.0;
    xx = 1.0;
    bx = 2.0;
    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, f1dim);
    *fret = brent(ax, xx, bx, f1dim, TOL, &xmin);
    for (j = 0; j < n; j++)
    {
        xi[j] *= xmin;
        p[j] += xi[j];
    }
    free((double *) xicom);
    free((double *) pcom);
}

static double f1dim(double x)
{
    extern int ncom;    /* defined in LINMIN */
    extern double *pcom, *xicom;

    int j;
    double f, *xt;

    xt = (double *) malloc(ncom * sizeof(double));

    for (j = 0; j < ncom; j++) xt[j] = pcom[j] + x * xicom[j];
    f = (*nrfunc)(xt);
    free((double *) xt);
    return f;
}

static void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
                   double *fc, double (*func)(double))
{
    double ulim, u, r, q, fu, dum;

    *fa = (*func)(*ax);
    *fb = (*func)(*bx);
    if (*fb > *fa)
    {
        SHFT(dum, *ax, *bx, dum)
        SHFT(dum, *fb, *fa, dum)
    }
    *cx = (*bx) + GOLD * (*bx - *ax);
    *fc = (*func)(*cx);
    while (*fb > *fc)
    {
        r = (*bx - *ax) * (*fb - *fc);
        q = (*bx - *cx) * (*fb - *fa);
        u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) /
            (2.0 * SIGN(MAX(fabs(q - r), TINY), q - r));
        ulim = (*bx) + GLIMIT * (*cx - *bx);
        if ((*bx - u)*(u - *cx) > 0.0)
        {
            fu = (*func)(u);
            if (fu < *fc)
            {
                *ax = (*bx);
                *bx = u;
                *fa = (*fb);
                *fb = fu;
                return;
            }
            else if (fu > *fb)
            {
                *cx = u;
                *fc = fu;
                return;
            }
            u = (*cx) + GOLD * (*cx - *bx);
            fu = (*func)(u);
        }
        else if ((*cx - u)*(u - ulim) > 0.0)
        {
            fu = (*func)(u);
            if (fu < *fc)
            {
                SHFT(*bx, *cx, u, *cx + GOLD*(*cx - *bx))
                SHFT(*fb, *fc, fu, (*func)(u))
            }
        }
        else if ((u - ulim)*(ulim - *cx) >= 0.0)
        {
            u = ulim;
            fu = (*func)(u);
        }
        else
        {
            u = (*cx) + GOLD * (*cx - *bx);
            fu = (*func)(u);
        }
        SHFT(*ax, *bx, *cx, u)
        SHFT(*fa, *fb, *fc, fu)
    }
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
#undef SIGN
#undef SHFT
#undef ITMAX
#undef EPS
#undef TOL
