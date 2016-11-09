#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

const extern lensdata *lens_table;

/****************************************************************/
/*              nom:            e_nfwgt                         */
/*              auteur:         Jean-Paul Kneib  + David Sand   */
/*              date:              12/05                        */
/*              place:          Toulouse                        */
/****************************************************************/
/* Global variables used :
 * in tnfwg_dpl() : lens_table
 */
double   nfwg_dpl(double r, double rs, double kappas, double alpha)
{
    // I divide by 2 because in lenstool_tab, it is multiplied by 2
    double dpl = kappas * rs / r * tnfwg_dpl(r / rs, alpha) / 2.;
    return( dpl );
}

/* --------------------------------------------------------------*
 * Global variable used :
 * - in tnfwg_kappa() : lens_table
 */
double   nfwg_kappa(double r, double rs, double kappas, double alpha)
{
    double kappa = kappas / 2. / rs * tnfwg_kappa(r / rs, alpha);
    return( kappa );
}

/* --------------------------------------------------------------*
 * Global variables used :
 * - in nfwg_kappa() : lens_table
 * - in nfwg_dpl() : lens_table
 */
double   nfwg_gamma(double r, double rs, double kappas, double alpha)
{
    return( nfwg_dpl(r, rs, kappas, alpha) / r - nfwg_kappa(r, rs, kappas, alpha) );
}

/* --------------------------------------------------------------*
 * Global variables used :
 * in nfwg_kappa() : lens_table
 * in nfwg_gamma() : lens_table
 */
double   nfwg_kappa_eps(double r, double rs, double theta, double kappas, double eps, double alpha)
{
    double kappa_eps;

    kappa_eps = nfwg_kappa(r, rs, kappas, alpha) + eps * cos(2 * theta)
                * nfwg_gamma(r, rs, kappas, alpha);


    return(kappa_eps);
}

/* --------------------------------------------------------------*
 * Return gamma1 and gamma2 for a pseudo-elliptical potential
 * Global variables used :
 * - in nfwg_dpl() : lens_table
 * - in nfwg_kappa() : lens_table
 * - in nfwg_gamma() : lens_table
 * - in nfwg_kappa_eps() : lens_table
 */

struct point   nfwg_gamma_eps(double r, double rs, double theta, double kappas, double eps, double alpha)
{
    struct point gamma_eps;
    double kappa, gamma;

    kappa = nfwg_kappa(r, rs, kappas, alpha);
    gamma = nfwg_gamma(r, rs, kappas, alpha);

//    gamma_eps = sqrt(gamma * gamma + 2 * eps * cos(2 * theta) * kappa * gamma + eps * eps * (kappa * kappa - pow(sin(2 * theta) * gamma, 2.)));
//
    gamma_eps.x = gamma * cos(2. * theta) + eps * kappa;  // gamma1 
    gamma_eps.y = -sqrt(1. - eps * eps) * gamma * sin(2. * theta); // gamma2

    return(gamma_eps);
}

/* --------------------------------------------------------------*
 * Return the NFW f function of x and alpha
 * Global variables used :
 * - lens_table
 */
double tnfwg_kappa(double xwant, double alphawant)
{
    long int xgoint, alphagoint, index_found, index_found2;
    long int n_xstepsint, intnalphasteps;
    double alphago, xgo, kappa_1, kappa_2, slope; //,intercept,dpl_1,dpl_2,kappatest,dpl_1test;
    double alphalow, alphastep, a_step, x_not, x_max, alphahigh;
    double nfwg_kappa_want;

    /*read in the first two lines of the array to find the parameters for the rang
    e of x and alpha*/
    alphalow = lens_table[0].alpha_now;
    intnalphasteps = (long int)lens_table[0].x_now;
    alphastep = lens_table[0].kappa;
    alphahigh = alphalow + alphastep * (intnalphasteps - 1);

    x_not = lens_table[1].alpha_now;
    n_xstepsint = (long int)lens_table[1].x_now;
    a_step = lens_table[1].kappa;
    x_max = x_not * pow(a_step, n_xstepsint - 1);

    /*first find the two neighboring index positions in lens_table
            from xwant and al phawant*/


    if (alphawant > alphahigh)
    {
//      fprintf(stderr,"alpha is out of bounds!!!! %lf\n",alphawant);
        alphawant = alphahigh - alphastep / 2;
    }
    else if (alphawant < alphalow)
    {
//      fprintf(stderr,"alpha is out of bounds!!!! %lf\n",alphawant);
        alphawant = alphalow + alphastep / 2;
    }

    alphago = (alphawant - alphalow) / alphastep;
    alphagoint = (long int)alphago;

    if (xwant > x_max)
    {
        xwant = x_max / sqrt(a_step);
    }
    else if (xwant < x_not)
    {
//      fprintf(stderr,"xwant is out of bounds!!!! %lf\n",xwant);
        xwant = x_not * sqrt(a_step);
    }
    xgo = log10(xwant / x_not) / log10(a_step);
    xgoint = (long int)xgo;


    /*the '+2' accounts for the first two lines being used to determine the size of the file*/

    index_found = alphagoint * n_xstepsint + xgoint + 2;
    index_found2 = index_found + n_xstepsint;  //+2;

    /*first find value of nfwg_kappa_want at first interpolation step*/

    slope = (lens_table[index_found].x_now - lens_table[index_found+1].x_now) / (lens_table[index_found].kappa - lens_table[index_found+1].kappa);

    kappa_1 = lens_table[index_found].kappa - (lens_table[index_found].x_now - xwant) / slope;


    /*now find value of nfwg_kappa_want at second interpolation step*/

    slope = (lens_table[index_found2].x_now - lens_table[index_found2+1].x_now)
            / (lens_table[index_found2].kappa - lens_table[index_found2+1].kappa);

    kappa_2 = lens_table[index_found2].kappa - (lens_table[index_found2].x_now - xwant)
              / slope;

    /*now interpolate between the two values of kappa (kappa_1 and kappa_2) that i
     just found*/

    slope = alphastep / (kappa_2 - kappa_1);

    nfwg_kappa_want = kappa_1 - (lens_table[index_found].alpha_now - alphawant) / slope;

    return nfwg_kappa_want;
}

/* -------------------------------------------------------------------- *
 * Return the NFW g function of x and alpha.
 * Global variables used :
 * - lens_table
 * */
double tnfwg_dpl(double xwant, double alphawant)
{
    long int xgoint, alphagoint, index_found, index_found2;
    long int n_xstepsint, intnalphasteps;

    double alphago, xgo, dpl_1, dpl_2, slope;
    double alphalow, alphastep, a_step, x_not, x_max, alphahigh;
    double nfwg_dpl_want;

    /*read in the first two lines of the array to find the parameters for
     * the range of x and alpha*/
    alphalow = lens_table[0].alpha_now;
    alphastep = lens_table[0].kappa;
    intnalphasteps = (long int)lens_table[0].x_now;
    alphahigh = alphalow + alphastep * (intnalphasteps - 1);

    x_not = lens_table[1].alpha_now;
    a_step = lens_table[1].kappa;
    n_xstepsint = (long int)lens_table[1].x_now;
    x_max = x_not * pow(a_step, n_xstepsint - 1);

    /*first find the two neighboring index positions in lens_table from xwant and alphawant*/

    if (alphawant > alphahigh)
    {
//      fprintf(stderr,"alpha is out of bounds!!!! %lf\n",alphawant);
        alphawant = alphahigh - alphastep / 2;
    }
    else if (alphawant < alphalow)
    {
//      fprintf(stderr,"alpha is out of bounds!!!! %lf\n",alphawant);
        alphawant = alphalow + alphastep / 2;
    }

    alphago = (alphawant - alphalow) / alphastep;
    alphagoint = (long int)alphago;


    if (xwant > x_max)
    {
//      fprintf(stderr,"xwant is out of bounds!!!! %lf\n",xwant);
        xwant = x_max / sqrt(a_step);
    }

    if (xwant < x_not)
    {
//      fprintf(stderr,"xwant is out of bounds!!!! %lf\n",xwant);
        xwant = x_not * sqrt(a_step);
    }

    xgo = log10(xwant / x_not) / log10(a_step);
    xgoint = (long int)xgo;

    /*the '+2' accounts for the first two lines being used to determine the size of the file*/
    index_found = alphagoint * n_xstepsint + xgoint + 2;
    index_found2 = index_found + n_xstepsint; //+2;

    /*first find value of nfwg_dpl_want at first interpolation step*/
    slope = (lens_table[index_found].x_now - lens_table[index_found+1].x_now) / (lens_table[index_found].dpl - lens_table[index_found+1].dpl);

    dpl_1 = lens_table[index_found].dpl - (lens_table[index_found].x_now - xwant) / slope;

    /*now find value of nfwg_dpl_want at second interpolation step*/
    slope = (lens_table[index_found2].x_now - lens_table[index_found2+1].x_now) / (lens_table[index_found2].dpl - lens_table[index_found2+1].dpl);

    dpl_2 = lens_table[index_found2].dpl - (lens_table[index_found2].x_now - xwant) / slope;


    /*now interpolate between the two values of dpl (dpl_1 and dpl_2) that i just found*/

    slope = alphastep / (dpl_2 - dpl_1);

    nfwg_dpl_want = dpl_1 - (lens_table[index_found].alpha_now - alphawant) / slope;
    return nfwg_dpl_want;
}
