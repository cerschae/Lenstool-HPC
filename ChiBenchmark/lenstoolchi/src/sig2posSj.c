#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        sig2posSj           */
/*      auteur:     Ghislain Golse          */
/*      date:       10/99               */
/*      place:      Toulouse            */
/****************************************************************/
/* Calcul de sig2pos dans le plan source
 *
 * Global variables used :
 * - amplifi
 */
double   sig2posSj(double sig2pos, struct galaxie *multi, int n, int j, int n_famille)
{
//  const extern struct   g_mode M;
//  const extern  struct  pot lens[NLMAX];
//  const extern  struct  g_cosmo     C;
    const extern  double   amplifi[NIMAX][NIMAX];

    double  sigma2, ainv;

    /* sigma=sqrt(fabs(e_amp(multi[j].C,dlsds))); */

//    ainv = 1. / fabs(amplifi[n_famille][j]);
    ainv = fabs(e_amp_gal(&multi[j], NULL)); //1. / fabs(amplifi[n_famille][i]);

    sigma2 = ainv * sig2pos;

    /* printf("s=%.3lf Sigma2=%.9lf\n",s,sigma2); */

    return(sigma2);

}
