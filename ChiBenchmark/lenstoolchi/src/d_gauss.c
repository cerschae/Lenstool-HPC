#include<math.h>
#include<fonction.h>

/* -------------------------------------------------------------------
double  gauss()

purpose: draw numbers under an gaussian law
-------------------------------------------------------------------*/

double d_gauss(double sig, int *idum)
{
    double x, z, r;

    while ( (z = pow(x = d_random(idum) - 0.5, 2.0) + pow(d_random(idum) - 0.5, 2.0) ) > 0.25 );
    while ((r = d_random(idum)) <= 0.0);

    return  sig*sqrt(-2.0*log(r) / z)*x;

}
