#include<math.h>
#include<fonction.h>

/* -------------------------------------------------------------------
return a random B absolute magnitude according to Schechter (1976)
distribution.
-------------------------------------------------------------------*/

double d_rndschechter(int *idum)
{
    const extern  struct  g_cosmo C;
    const extern  struct  g_source S;
    double mrnd, prnd, dm;

    dm = 5.0 * log10(C.h); /* C.h is in units of H0=50 */
    do
    {
        mrnd = S.lfm_min + (S.lfm_max - S.lfm_min) * d_random(idum);
        prnd = 2.0 * d_random(idum);
    }
    while (exp(0.921*(S.lfalpha + 1)*(S.lfm_star - mrnd)
               - exp(0.921*(S.lfm_star - mrnd))) < prnd);

    return  mrnd + dm;
}
