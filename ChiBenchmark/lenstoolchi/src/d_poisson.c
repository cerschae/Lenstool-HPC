#include <stdio.h>
#include <math.h>
#include <fonction.h>

#define PI 3.141592653589793238462643

static double lngam(double a);

double d_poisson(double xm, int *idum)
{
    static double sq, alxm, g, oldm = (-1.0);
    double em, t , y;

    if (xm < 12.0)
    {
        if (xm != oldm)
        {
            oldm = xm;
            g = exp(-1.*xm);
        };
        em = (-1);
        t = 1.0;
        do
        {
            em += 1.0;
            t *= d_random(idum);
        }
        while (t > g);
    }
    else
    {
        if (xm != oldm)
        {
            oldm = xm;
            sq = sqrt(2.0 * xm);
            alxm = log(xm);
            g = xm * alxm - lngam(xm + 1.0);
        };
        do
        {
            do
            {
                y = tan(PI * d_random(idum));
                em = sq * y + xm;
            }
            while (em < 0.0);
            em = floor(em);
            t = .9 * (1.0 + y * y) * exp(em * alxm - lngam(em + 1.0) - g);
        }
        while (d_random(idum) > t);
    };
    return (em);
}

/******************************************************/

static double lngam(double a)
{
    double cof[6], stp, half = .5, one = 1., fpf = 5.5, x, tmp, ser;
    double pi;
    double atemp;
    double gama, gama1;
    int i;

    pi = 3.141592654;
    atemp = 10.;
    stp = 2.50662827465;
    cof[0] = 76.18009173;
    cof[1] = (-86.50532033);
    cof[2] = 24.01409822;
    cof[3] = (-1.231739516);
    cof[4] = 0.120858003e-2;
    cof[5] = (-0.536382e-5);

    if (a == 1)
        gama = 1.;
    else
    {
        if (a < 1.)
        {
            atemp = 1.;
            a = 2. - a;
        };
        x = a - one;
        tmp = x + fpf;
        tmp = (x + half) * log(tmp) - tmp;
        ser = one;
        for (i = 0; i < 6; i++)
        {
            x += one;
            ser += cof[i] / x;
        };
        gama1 = tmp + log(stp * ser);
        gama = gama1;
        if (atemp < 1.)
        {
            pi = pi * (a - 1.);
            gama = log(pi / sin(pi)) - gama;
            a = atemp;
        };
    };
    return(gama);
}

