#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

static double   nfw_f(double r);
static double   nfw_g(double r);

/****************************************************************/
/*              nom:            e_nfw                           */
/*              auteur:         Jean-Paul Kneib                 */
/*              date:              06/99                        */
/*              place:          Toulouse                        */
/****************************************************************
 * r and rs are in radians.
 *
 * theta(in ") = r * DOL
 */

/* Return the displacement in radians, ie the gradient of the projected
 * lens potential for a circular NFW potential.
 *
 * Warning : You have to multiply by DLS/DS to get the right value
 *
 * in GK2002 Eq 6: dpl(x) = 4*kappas*r(")/x^2 * g(x)
 *               kappa(x) = 2*kappas
 *       with      kappas = rho_s*rs*Sigma_crit^-1
 *       and      vdisp^2 = 8/3*G*rs^2*rho_s (the 8/3 factor scales for the LT b0 definition)
 *       and       x = r/rs
 *
 * ==> 4*kappas*r(")/x^2 = (6*PI/c^2 * vdisp^2) * DLS/DS * rs/r
 *     2*kappas = (6*PI/c^2 * vdisp^2) * DLS/DS * 1 / 2rs
 *
 * here : kappas_LT = (6*PI/c^2 * vdisp^2) (see: set_lens.c, e_grad.c)
 * ==> dpl(x) = kappas_LT * rs/r * nfw_g(r/rs) * DLS/DS
 * ==> kappa(x) = kappas_LT / 2rs * nfw_f(r/rs) * DLS/DS
 *
 * Global variables used :
 * - none
 * */
double   nfw_dpl(double r, double rs, double kappas)
{
    double dpl = kappas * rs / r * nfw_g(r / rs);
    return( dpl );
}

/* --------------------------------------------------------------*/
/* Global variables used :
 * - none
 * */
double   nfw_kappa(double r, double rs, double kappas)
{
    double kappa = kappas / 2. / rs * nfw_f(r / rs);
    return( kappa );
}

/* --------------------------------------------------------------*
 * Return kappa for elliptical NFW (cf. Golse & Kneib 2002, eq 17)
 * Global variables used :
 * - none
 * */
double   nfw_kappa_eps(double r, double rs, double theta, double kappas, double eps)
{
    double kappa_eps;

    kappa_eps = nfw_kappa(r, rs, kappas) + eps * cos(2 * theta) * nfw_gamma(r, rs, kappas);

    return(kappa_eps);
}

/* --------------------------------------------------------------*
 * Return gamma1 and gamma2 for elliptical NFW 
 * (cf. Drumet-Montoyat et al. 2012 Eq A10-A12, Golse & Kneib 2002, eq 18)
 * Global variables used :
 * - none
 */
struct point   nfw_gamma_eps(double r, double rs, double theta, double kappas, double eps)
{
    struct point gamma_eps;
    double kappa, gamma;

    kappa = nfw_kappa(r, rs, kappas);
    gamma = nfw_gamma(r, rs, kappas);

//    gamma_eps = sqrt( pow(gamma, 2.) + 2 * eps * cos(2 * theta) * kappa * gamma + eps * eps * ( pow(kappa, 2.) - pow(sin(2 * theta) * gamma, 2.) ) );  //NOT USED ANYMORE

    gamma_eps.x = gamma * cos(2. * theta) + eps * kappa;  // gamma1 
    gamma_eps.y = -sqrt(1. - eps * eps) * gamma * sin(2. * theta); // gamma2

    return(gamma_eps);
}

/* --------------------------------------------------------------*/
/* Global variables used :
 * - none
 * */
double   nfw_gamma(double r, double rs, double kappas)
{
    /* double nfw_g();
       double nfw_f();

    return( kappas/rs*( nfw_g(r/rs)/r/r*rs*rs -nfw_f(r/rs)/2./(r*r/rs/rs-1) )
            );
    */
    return( nfw_dpl(r, rs, kappas) / r - nfw_kappa(r, rs, kappas) );

}

/* --------------------------------------------------------------*/
/* Global variables used :
 * - none
 * */
static double   nfw_f(double r)
{
    double f = 0;

    if ( r > 1 )
        f = ( 1. - 2. / sqrt(r * r - 1) * atan(sqrt((r - 1.) / (r + 1.))) ) / (r * r - 1);    // same as acos(1/r)
    else if ( r < 1 )
        f = ( 1. - 2. / sqrt(1 - r * r) * atanh(sqrt((1. - r) / (r + 1.))) ) / (r * r - 1);   // same as acosh(1/r)
    else
        f = 1. / 3.;

    return(f);
}

/* --------------------------------------------------------------*/
/* Global variables used :
 * - none
 * */
static double   nfw_g(double r)
{
    double g = 1;

    if ( r > 1 )
        g = log(r / 2.) + 2. / sqrt(r * r - 1) * atan(sqrt((r - 1.) / (r + 1.)));   // same as arcos(1/r)
    else if ( r < 1 )
        g = log(r / 2) + 2. / sqrt(1 - r * r) * atanh(sqrt((1. - r) / (r + 1.)));   // same as arcosh(1/r)
    else
        g = 1. + log(0.5);

    return(g);
}

/* --------------------------------------------------------------
 * Conversion rs,sigmas <--> concentration
 * Tomas Verdugo Gonzalez :  Wed, 24 Jan 2007 22:07:18 -0600
 */

/* --------------------------------------------------------------*
 * Return the critical density of the Universe in Msol / Mpc3
 */

double rho_cri(double z)
{
    const extern  struct  g_cosmo C;
    const extern  struct  pot lens[];
    double Hcuad, rho_cri;

    Hcuad = C.H0 / chiz(z);
    Hcuad = Hcuad * Hcuad;  //  H^2 (z)

    rho_cri = 3.0 / 8.0 / PI * INVG * Hcuad;   //Critical density

    return (rho_cri);
}

/* --------------------------------------------------------------*
 * Return the Overdensity (rho-rho_cri)/rho_cri (Lacey & Cole 1993; Eke et al. 1996)
 */
double DDelta_c()
{
    const extern  struct  g_cosmo C;
    const extern  struct  pot lens[];
    double  DDelta_c = 0;
    double  E, Omega;

    E = 1. / chiz(lens[0].z);
    Omega = C.omegaM * pow( 1.0 + lens[0].z, 3) / E / E; //Omega at time z

    if (  C.omegaX == 0.0  )
        DDelta_c =  178.0 *  pow(Omega, 0.30);                  //Overdensity for lambda=0
    else if ( C.omegaX + C.omegaM == 1 )
        DDelta_c =  178.0 *  pow(Omega, 0.45);                  //Overdensity for Omega+lambda=1
    // TODO else

    return (DDelta_c);
}

/* --------------------------------------------------------------*
 * Bisection method to find the solution of Fc
 */
#define JMAX 40     //Maximum allowed number of bisections.
double rtbis2(double (*func)(double), double x1, double x2, double xacc)
//Using bisection, ﬁnd the root of a function func known to lie between x1 and x2. The root,
//returned as rtbis, will be reﬁned until its accuracy is ±xacc.
{
    int j;
    double dx, f, fmid, xmid, rtb;

    f = (*func)(x1);
    fmid = (*func)(x2);
    if (f*fmid >= 0.0)  fprintf(stderr, "ERROR: (NFW:rtbis2) Root must be bracketed for bisection\n");
    rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
    for (j = 1; j <= JMAX; j++)
    {
        fmid = (*func)(xmid = rtb + (dx *= 0.5));
        if (fmid <= 0.0) rtb = xmid;
        if (fabs(dx) < xacc || fmid == 0.0) return rtb;
    }
    fprintf(stderr, "ERROR: (NFW:rtbis2) Too many bisections\n");
    return 0.0;
}
#undef JMAX

static double delta_c_ext; //Characteristic density contrast

/* --------------------------------------------------------------*
* The next function return de F(c) value
* Global variables used :
*
*/
static double efe_c(double cc)
{
    double A, Fc;

    A = cc * cc * cc / (log(1. + cc) - cc / (1. + cc));
    Fc = delta_c_ext -  200. / 3.0 * A; // TODO : DDelta()

    return (Fc);
}

/* --------------------------------------------------------------*
 * Return one formulation of NFW given sigma_s and r_s
 * :
 *
 */
void  e_nfw_rs2c(double sigma_s, double r_s, double *rhos, double *c, double *M_200, double z)
{

    extern double delta_c_ext;
    double rho_c;
    double r_200;
    double r_s_Mpc;
    double a = 0.0001;
    double b = 30.0;
    double err = 0.0001;

    rho_c = rho_cri(z);
    r_s_Mpc = r_s / 1000.;    //rs in Megaparsecs
    *rhos =  (3. / 8. * INVG) * sigma_s * sigma_s / r_s_Mpc / r_s_Mpc;  
    delta_c_ext = *rhos / rho_c;

    *c = rtbis2(efe_c, a, b, err);   //This is the concentration value

    r_200 = (*c) * r_s_Mpc ;    //The Radius 200xrhocrit

    //Virial mass in solar masses
    *M_200 = (4.0 / 3.0) * PI * r_200 * r_200 * r_200 * rho_c * 200.; //DDelta_c(); 
}

/* --------------------------------------------------------------*
 * Return one formulation of NFW given c and rs (in kpc)
 * Global variables used :
 *
 */
void e_nfw_crs2sig(double c, double rs, double *sigma_s, double z)
{
    double sigma_s0, rhos;

    rhos = 200. / 3. * c * c * c / ( log( 1 + c ) - c / ( 1 + c ) ) * rho_cri(z); // in Msol/Mpc3

    rs /= 1000.;   // in Mpc
    sigma_s0 = 8. / 3. / INVG * rs * rs * rhos;
    *sigma_s = sqrt(sigma_s0);  // in km/s
}

/* --------------------------------------------------------------*
 * Return one formulation of NFW given c and m200 (in Msol) to
 * sigma_s in km/s and rs in kpc. (See NFW1996, Eq 4)
 * Global variables used :
 *
 */
void e_nfw_cm200_sigrs( double c, double m200, double *sigma_s, double *r_s, double z )
{
    double rhos, r200, sigma_s0, rho_c;

    rho_c = rho_cri(z);
    rhos = 200. / 3. * c * c * c / ( log( 1 + c ) - c / ( 1 + c ) ) * rho_c; // in Msol/Mpc3
    r200 = 3. * m200 / 4. / PI / rho_c / 200.;  // in Mpc^3
    r200 = pow( r200, 0.333333);

    *r_s = r200 / c;  // in Mpc
    sigma_s0 = 8.0 / INVG * (*r_s) * (*r_s) * rhos / 3.0;
    *sigma_s = sqrt(sigma_s0);              //in km/s

    *r_s *= 1000.;  // in kpc
}

/* --------------------------------------------------------------*
 * Return one formulation of NFW given c and r200 (in kpc) to
 * sigma_s in km/s and rs in kpc.
 * Global variables used :
 *
 */
void e_nfw_cr200_sigrs( double c, double r200, double *sigma_s, double *r_s, double z )
{
    double rhos, sigma_s0;

    rhos = 200. / 3. * c * c * c / ( log( 1 + c ) - c / ( 1 + c ) ) * rho_cri(z); // in Msol/Mpc3
    r200 /= 1000.; // in Mpc

    *r_s = r200 / c;  // in Mpc

    sigma_s0 = 8. / 3. / INVG * (*r_s) * (*r_s) * rhos;
    *sigma_s = sqrt(sigma_s0);              //in km/s

    *r_s *= 1000.;  // in kpc
}

/* --------------------------------------------------------------*
 * Return one formulation of NFW given rs (in kpc) and m200 (in Msol) to
 * sigma_s in km/s.
 * Global variables used :
 *
 */
void e_nfw_rsm200_sigrs( double rs, double m200, double *sigma_s, double z )
{
    double r200, c, rhos, sigma_s0, rho_c;

    rho_c = rho_cri(z);
    r200 = 3. * m200 / 4. / PI / rho_c / 200.;  // in Mpc
    r200 = pow( r200, 0.333333);

    c = r200 / rs * 1000.;

    rhos = 200. / 3. * c * c * c / ( log( 1 + c ) - c / ( 1 + c ) ) * rho_c; // in Msol/Mpc3

    sigma_s0 = 8. / 3. / INVG * rs * rs * rhos;
    *sigma_s = sqrt(sigma_s0);              // in km/s
}

/* --------------------------------------------------------------*
 * Return one formulation of NFW given rs and r200 (in kpc) to
 * sigma_s in km/s.
 * Global variables used :
 *
 */
void e_nfw_rsr200_sigrs( double rs, double r200, double *sigma_s, double z )
{
    double c, rhos, sigma_s0;

    c = r200 / rs;

    rhos = 200. / 3. * c * c * c / ( log( 1 + c ) - c / ( 1 + c ) ) * rho_cri(z); // in Msol/Mpc3

    sigma_s0 = 8. / 3. / INVG * rs * rs * rhos;
    *sigma_s = sqrt(sigma_s0);              // in km/s
}

