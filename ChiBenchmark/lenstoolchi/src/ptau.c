#include<math.h>

double   ptau(double x, double y)
{
    return( exp(-(x*x + y*y) / 0.0625));
}
