#include <stdio.h>
#include "fonction.h"
#include "lt.h"

void bicubic_int_z(double y[], double y1[], double y2[], double y12[],
                   double x1l, double x1u, double x2l, double x2u,
                   double x1, double x2,
                   double *ansy)
{
    int i;
    double t, u, d1, d2, **c;

    c = alloc_square_double(4, 4);
    d1 = x1u - x1l;
    d2 = x2u - x2l;
    bicubic_coef(y, y1, y2, y12, d1, d2, c);
    if (x1u == x1l || x2u == x2l)
        fprintf(stderr, "Bad input in routine bicubic_int");
    t = (x1 - x1l) / d1;
    u = (x2 - x2l) / d2;
    *ansy = 0.0;
    for (i = 3; i >= 0; i--)
        *ansy = t * (*ansy) + ((c[i][3] * u + c[i][2]) * u + c[i][1]) * u + c[i][0];
    free_square_double(c, 4);
}


/* gradient */

void bicubic_int_zz(double y[], double y1[], double y2[], double y12[],
                    double x1l, double x1u, double x2l, double x2u,
                    double x1, double x2,
                    double *ansy1, double *ansy2)
{
    int i;
    double t, u, d1, d2, **c;

    c = alloc_square_double(4, 4);
    d1 = x1u - x1l;
    d2 = x2u - x2l;
    bicubic_coef(y, y1, y2, y12, d1, d2, c);
    if (x1u == x1l || x2u == x2l)
        fprintf(stderr, "Bad input in routine bicubic_int");
    t = (x1 - x1l) / d1;
    u = (x2 - x2l) / d2;
    (*ansy2) = (*ansy1) = 0.0;
    for (i = 3; i >= 0; i--)
    {
        *ansy2 = t * (*ansy2) + (3.0 * c[i][3] * u + 2.0 * c[i][2]) * u + c[i][1];
        *ansy1 = u * (*ansy1) + (3.0 * c[3][i] * t + 2.0 * c[2][i]) * t + c[1][i];
    }

    *ansy1 /= d1;
    *ansy2 /= d2;
    free_square_double(c, 4);
}

/* second derivatives */

void bicubic_int_zxy(double y[], double y1[], double y2[], double y12[],
                     double x1l, double x1u, double x2l, double x2u,
                     double x1, double x2,
                     double *ansy12)
{
    int i;
    double t, u, d1, d2, **c;

    c = alloc_square_double(4, 4);
    d1 = x1u - x1l;
    d2 = x2u - x2l;
    bicubic_coef(y, y1, y2, y12, d1, d2, c);
    if (x1u == x1l || x2u == x2l)
        fprintf(stderr, "Bad input in routine bicubic_int");
    t = (x1 - x1l) / d1;
    u = (x2 - x2l) / d2;

    (*ansy12) = 0.;
    for (i = 3; i >= 1; i--)
        *ansy12 = t * (*ansy12) + (3.*i * c[i][3] * u + 2.*i * c[i][2]) * u + i * c[i][1];

    *ansy12 /= (d1 * d2);

    free_square_double(c, 4);
}

/* everything */

void bicubic_int(double y[], double y1[], double y2[], double y12[],
                 double x1l, double x1u, double x2l, double x2u,
                 double x1, double x2,
                 double *ansy, double *ansy1, double *ansy2,
                 double *ansy11, double *ansy12, double *ansy22)
{
    int i;
    double t, u, d1, d2, **c;

    c = alloc_square_double(4, 4);
    d1 = x1u - x1l;
    d2 = x2u - x2l;
    bicubic_coef(y, y1, y2, y12, d1, d2, c);
    if (x1u == x1l || x2u == x2l)
        fprintf(stderr, "Bad input in routine bicubic_int");
    t = (x1 - x1l) / d1;
    u = (x2 - x2l) / d2;
    *ansy = (*ansy2) = (*ansy1) = 0.0;
    (*ansy22) = (*ansy11) = 0.0;
    for (i = 3; i >= 0; i--)
    {
        *ansy = t * (*ansy) + ((c[i][3] * u + c[i][2]) * u + c[i][1]) * u + c[i][0];
        *ansy2 = t * (*ansy2) + (3.0 * c[i][3] * u + 2.0 * c[i][2]) * u + c[i][1];
        *ansy1 = u * (*ansy1) + (3.0 * c[3][i] * t + 2.0 * c[2][i]) * t + c[1][i];
        *ansy22 = t * (*ansy22) + (3.*c[i][3] * u + c[i][2]) * 2.;
        *ansy11 = u * (*ansy11) + (3.*c[3][i] * t + c[2][i]) * 2.;
    }
    (*ansy12) = 0.;
    for (i = 3; i >= 1; i--)
        *ansy12 = t * (*ansy12) + (3.*i * c[i][3] * u + 2.*i * c[i][2]) * u + i * c[i][1];

    *ansy1 /= d1;
    *ansy2 /= d2;
    *ansy11 /= (d1 * d1);
    *ansy12 /= (d1 * d2);
    *ansy22 /= (d2 * d2);
    free_square_double(c, 4);
}
