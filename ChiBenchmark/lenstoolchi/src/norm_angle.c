#include<stdio.h>
#include<math.h>
#include "fonction.h"
#include "constant.h"
#include"dimension.h"
#include "structure.h"

double  norm_angle(double theta, double npi)
{
    double N;
    double  t, d;
    int n;

    n = (int) (theta / npi);
    N = n;
    d = theta - N * npi;
    if (d > .5*npi)
        N += 1.;

    t = theta - N * npi;
    return(t);
}
