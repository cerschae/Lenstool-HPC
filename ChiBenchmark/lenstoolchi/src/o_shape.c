#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

static double   est_sx(double tau, double t, double dis, double tp, double tm);
static double   est_sy(double tau, double t, double dis, double tp, double tm);


/********************************************************/
/*      fonction:   o_shape             */
/*      auteur:         jpk         */
/********************************************************
 * Return the differences (dx, dy, da) between the distortion-deformation
 * coefficients and equivalent surfaces (2*a*b) of arclet i and
 * arclet i+1 in the source plane.
 *
 * sigx2 and sigy2 are set to 0.6.
 *
 * Global variables used :
 * - none
 */
void    o_shape(int i, struct galaxie *gali,
                double *dx, double *sigx2, double *dy, double *sigy2, double *da)
{
    int j;

    j = i + 1;

    if (gali[i].E.a*gali[j].E.a != 0)
    {
        *dx = est_sx(gali[i].tau, gali[i].E.theta, gali[i].dp, gali[i].tp, gali[i].thp)
              - est_sx(gali[j].tau, gali[j].E.theta, gali[j].dp, gali[j].tp, gali[j].thp);
        *sigx2 = 0.6;

        *dy = est_sy(gali[i].tau, gali[i].E.theta, gali[i].dp, gali[i].tp, gali[i].thp)
              - est_sy(gali[j].tau, gali[j].E.theta, gali[j].dp, gali[j].tp, gali[j].thp);
        *sigy2 = 0.6;

        *da = gali[i].E.a * gali[i].E.b * gali[i].A
              - gali[j].E.a * gali[j].E.b * gali[j].A;
    }
    else
    {
        *da = *dx = *dy = 0;
        *sigx2 = 1.;
        *sigy2 = 1.;
    };
}
/****************************************************************/

void    o_dmag(int i, struct galaxie *gali, double *da)
{
    if (gali[i].E.a != 0)
        *da = gali[i].E.a * gali[i].A - gali[i+1].E.a * gali[i+1].A;
    else
        *da = 0;
}

void    o_shape_rond(int i, struct galaxie *gali,
                     double *dx, double *sigx2, double *dy, double *sigy2, double *da, double *siga2)
{
    if (gali[i].E.a*gali[i].E.b != 0)
    {
        *dx = est_sx(gali[i].tau, gali[i].E.theta, gali[i].dp, gali[i].tp, gali[i].thp);
        *sigx2 = .1;

        *dx = est_sy(gali[i].tau, gali[i].E.theta, gali[i].dp, gali[i].tp, gali[i].thp);
        *sigy2 = .1;

    }
    else
    {
        *da = *dx = *dy = 0;
        *sigx2 = 1.;
        *sigy2 = 1.;
        *siga2 = 1.;
    };
}


/****************************************************************
 * Compute the distortion-deformation matrix x coefficient source
 * plane.
 *
 * Parameters :
 * - tau
 * - t : arclet orientation (theta)
 * - dis : arclet distortion (K^2 + G^2)/(K^2 - G^2)
 * - tp : tauPot 2KG/(K^2 - G^2)
 * - tm : shear orientation (theta_pot)
 *
 * Global variables used :
 * - none
 */
static double   est_sx(double tau, double t, double dis, double tp, double tm)
{
    return( tau*sin(2.*(t - tm))*cos(2.*tm) +
            (dis*tau*cos(2.*(t - tm)) - tp*sqrt(tau*tau + 1.))*sin(2.*tm) );
}

/****************************************************************
 * Compute the distortion-deformation matrix y coefficient source
 * plane.
 *
 * Parameters :
 * - tau
 * - t : arclet orientation (theta)
 * - dis : arclet distortion (K^2 + G^2)/(K^2 - G^2)
 * - tp : tauPot 2KG/(K^2 - G^2)
 * - tm : shear orientation (theta_pot)
 *
 * Global variables used :
 * - none
 */
static double   est_sy(double tau, double t, double dis, double tp, double tm)
{
    return( -tau*sin(2.*(t - tm))*sin(2.*tm) +
            (dis*tau*cos(2.*(t - tm)) - tp*sqrt(tau*tau + 1.))*cos(2.*tm) );
}
