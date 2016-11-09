#include<math.h>
#include<fonction.h>

/* -------------------------------------------------------------------
generate a redshift for a source in a homogeneous, non-evolutive universe.
normalized so that z < zmax.
-------------------------------------------------------------------*/

double d_rndz(double zmax, int *idum)
{
    const extern  struct  g_cosmo C;
    double zrnd, dl, prnd;

    do
    {
        zrnd = pow(d_random(idum), 1.0 / 3.0) * zmax;      /* majorating func.: */
        prnd = zrnd * zrnd * d_random(idum);   /* Euclidean static approx.  */
        dl = dlumcosmo1(zrnd);
    }
    while ( dl*dl / (pow(1.0 + zrnd, 3.0)*sqrt(1.0 + C.omegaM*zrnd))
            < prnd);


    return  zrnd;
}
