#include<math.h>
#include<constant.h>
#include<fonction.h>

#undef PC
#define PC  3.085678e16             /* one parsec in MKS */

/* estimate magnitude knowing type,z and absolute magnitude,
as well as compute other important parameters
*/

void    s_compmag(struct galaxie *gal, int *idum)
{
    const extern  struct  g_cosmo C;
    double  dm, dl, dltpc, da;
    double  id, idc, ib, mabsb, mabsd;
    double  ar;
    double  color = 1;

    /*distance modulus difference with respect to h=0.5 */
    dm = 5.0 * log10(C.h);

    /*define luminosity distance dl as in Weinberg (1972) */
    dl = D0pc / C.h * dlumcosmo1(gal->z);
    dltpc = dl / 10.0; /* to use absolute magnitudes. */

    /*idem for angular distance da */
    da = dl / ((1 + gal->z) * (1 + gal->z));

    if (gal->type <= -4)      /*- pure elliptical profile */
    {
        mabsd = 0.0;
        mabsb = gal->magabs;
    }
    else               /* - B/D ratio from Simien (1988) */
    {
        if (gal->type == -3)
            mabsb = gal->magabs + 0.713;
        else
            mabsb = gal->magabs + 0.713 + 0.160 * gal->type + 0.042 * gal->type * gal->type;
        mabsd = -2.5 * log10(pow(10.0, -0.4 * gal->magabs)
                             - pow(10.0, -0.4 * mabsb));
    }

    /*compute axis ratio from inclination */
    /*  ar = fabs(cos(gal->inclination*DEG));
      if (ar<0.15)*/
    ar = 0.15;

    /*compute the total observed intensities */
    ib = pow(10.0, -0.4 * (mabsb - 1.5 * color)) / (dltpc * dltpc);
    if (gal->type > 0)
        id = pow(10.0, -0.4 * (mabsd - 0.7 * color)) / (dltpc * dltpc);
    else if (gal->type > -4)
        id = pow(10.0, -0.4 * (mabsd - 1.5 * color)) / (dltpc * dltpc);
    else
        id = 0.0;

    /*now the K+e-correction from Metcalfe et al. (1991) */
    ib *= pow(10.0, -0.4 * (      /* = E galaxy */
                  -0.109 + 0.061 * color + (5.700 - 5.463 * color) * gal->z
                  + (-4.919 + 7.518 * color) * gal->z * gal->z
                  + (-0.286 - 3.560 * color) * gal->z * gal->z * gal->z
                  + (0.558 + 0.603 * color) * gal->z * gal->z * gal->z * gal->z));

    if (gal->type > 7 && gal->type < 20)
    {
        id *= pow(10.0, -0.4 * (    /* = Sdm galaxy */
                      -0.003 + 0.003 * color
                      + (2.411 - 2.187 * color) * gal->z
                      + (-4.288 + 6.167 * color) * gal->z * gal->z
                      + (2.272 - 4.332 * color) * gal->z * gal->z * gal->z
                      + (-0.4158 + 0.940 * color) * gal->z * gal->z * gal->z * gal->z));

        /*disk extinction following De Vaucouleurs et al. 1991 et Cardelli et al. 1989 */
        idc = id * pow(ar,
                       -0.4 * (color < 0.5 ? 1.0 : 0.564) * (1.5 - 0.03 * (gal->type - 5.0) * (gal->type - 5.0)));
    }

    else if (gal->type > 0 && gal->type <= 7)
    {
        id *= pow(10.0, -0.4 * (    /* = Scd galaxy */
                      -0.037 + 0.062 * color
                      + (4.044 - 4.124 * color) * gal->z
                      + (-4.660 + 8.629 * color) * gal->z * gal->z
                      + (1.863 - 5.273 * color) * gal->z * gal->z * gal->z
                      + (-0.247 + 1.039 * color) * gal->z * gal->z * gal->z * gal->z));


        /*disk extinction following De Vaucouleurs et al. 1991 et Cardelli et al. 1989 */
        idc = id * pow(ar,
                       -0.4 * (color < 0.5 ? 1.0 : 0.564) * (1.5 - 0.03 * (gal->type - 5.0) * (gal->type - 5.0)));

    }
    else
    {
        id *= pow(10.0, -0.4 * (    /* = S0 galaxy */
                      -0.109 + 0.061 * color
                      + (5.700 - 5.463 * color) * gal->z
                      + (-4.919 + 7.518 * color) * gal->z * gal->z
                      + (-0.286 - 3.560 * color) * gal->z * gal->z * gal->z
                      + (0.558 + 0.603 * color) * gal->z * gal->z * gal->z * gal->z));
        idc = id;
    }

    gal->mag = -2.5 * log10(ib + idc);

}
