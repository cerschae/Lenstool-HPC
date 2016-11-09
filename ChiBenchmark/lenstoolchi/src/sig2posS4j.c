#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        sig2posS4j          */
/*      auteur:     Ghislain Golse          */
/*      date:       12/99               */
/*      place:      Toulouse            */
/****************************************************************
 * Calcul de sig2pos dans le plan source en faisant intervenir
 * 4 points distants de sigma_I de l'image I.
 *
 * Global variables used :
 * - in o_dpl() : G, lens, lens_table
 */
double   sig2posS4j(double sig2pos, struct galaxie *multi, struct point ps, int j)
{
//  const extern  struct   g_mode M;
//  const extern  struct  pot lens[NLMAX];
//  const extern  struct  g_cosmo     C;

    struct  galaxie Image[4];
    struct  point   Source[4];
    double  dlsds, sigma, d1, d;
    double  sigmaArcsec;
    int     i, n;

    n = 4;
    dlsds = multi[0].dr;
    Image[0].dr = dlsds;
    sigmaArcsec = sqrt(sig2pos);

    Image[0].C.x = multi[j].C.x;
    Image[0].C.y = multi[j].C.y - sigmaArcsec;
    Image[1].C.x = multi[j].C.x - sigmaArcsec;
    Image[1].C.y = multi[j].C.y;
    Image[2].C.x = multi[j].C.x;
    Image[2].C.y = multi[j].C.y + sigmaArcsec;
    Image[3].C.x = multi[j].C.x + sigmaArcsec;
    Image[3].C.y = multi[j].C.y;

    for ( i = 0; i < 4; i++ )
        Image[i].grad.x = Image[i].grad.y = 0;

    /* printf("  1/A0[%d]=%.3lf \n",j,e_amp(multi[j].C,dlsds)); */

    for (i = 0; i < n; i++)
    {
        /* sigma+=fabs(e_amp(multi[i].C,dlsds)); */
        /* NPRINTF(stderr,"  1/A[%d]=%.3lf \n",i,e_amp(Image[i].C,dlsds)); */
    };

    /* for(i=0;i<n;i++)
      NPRINTF(stderr,"Image[%d].x=%.3lf Image[%d].y=%.3lf \n",i,Image[i].C.x,i,Image[i].C.y); */

    o_dpl(n, Image, Source, NULL);

    /* for(i=0;i<n;i++)
      NPRINTF(stderr,"Source[%d].x=%.3lf Source[%d].y=%.3lf \n",i,Source[i].x,i,Source[i].y);  */


    d1 = dist(ps, Source[0]);
    /* NPRINTF(stderr,"d=%.3lf\n",d1); */

    for (i = 1; i < n; i++)
    {
        d = dist(ps, Source[i]);
        /* if(d>d1)
          d1=d; */
        d1 += d;
        /* NPRINTF(stderr,"d=%.3lf\n",d);*/
    };

    d1 /= 4;

    sigma = d1 * d1;

    /* NPRINTF(stderr,"Sigma2=%.3lf\n",sigma); */

    return(sigma);

}
