#include<stdio.h>
#include<math.h>

double  interpol(double xx, const double *fxx, const double *fyy, const double *fy2, int imax)
{
    int klo, khi, k;
    double yy, h, b, a;

    klo = 0;
    khi = imax - 1;
    while (khi - klo > 1)
    {
        k = (khi + klo) >> 1;
        if (fxx[k] > xx)
            khi = k;
        else
            klo = k;
    }
    h = fxx[khi] - fxx[klo];
    if (h == 0.0)
        fprintf(stderr, "Bad xx table's input to routine interpol\n");


    a = (fxx[khi] - xx) / h;
    b = (xx - fxx[klo]) / h;
    if (fy2[klo]<1e30 && fy2[klo]> -1e30 && fy2[khi] < 1e30 && fy2[khi] - 1e30)
    {
        yy = a * fyy[klo] + b * fyy[khi] +
             ((a * a * a - a) * fy2[klo] + (b * b * b - b) * fy2[khi]) * (h * h) / 6.0;
    }
    else
    {
        yy = a * fyy[klo] + b * fyy[khi];
    }

    fprintf(stderr, "%d %.3lf %d %.3lf %.3lf -> %.3lf %.3lf %.3lf\n",
            khi, fxx[khi], klo, fxx[klo], xx, fy2[khi], fy2[klo], yy);

    return(yy);
}

