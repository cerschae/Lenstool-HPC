#include<stdlib.h>
#include<math.h>
#include<fonction.h>

/* -------------------------------------------------------------------
return a random Hubble type between -5.0 and 10.0 following Spiekermann (1992).
-------------------------------------------------------------------*/

int d_rndtype(int *idum)
{
    double  rnd1, rnd2;

    rnd1 = d_random(idum);
    rnd2 = d_random(idum);

    if ((rnd1 -= 0.15) < 0.0)                     /* "pure" Ellipticals */
        return (int)(rnd2 * 2.0 - 6.0);
    else if ((rnd1 -= 0.21) < 0.0)                /* S0 galaxies */
        return (int)(rnd2 * 3.0 - 4.0);
    else if ((rnd1 -= 0.12) < 0.0)                /* S0/a, Sa */
        return (int)(rnd2 * 3.0);
    else if ((rnd1 -= 0.30) < 0.0)                /* Sab, Sb, Sbc */
        return (int)(rnd2 * 3.0 + 3.0);
    else                                  /* Sc, Scd, Sd, Sdm, Sm */
        return (int)(rnd2*5.0 + 6.0);

    fprintf(stderr, "RND: problem with rndtype()\n");
    exit (-1);

}
