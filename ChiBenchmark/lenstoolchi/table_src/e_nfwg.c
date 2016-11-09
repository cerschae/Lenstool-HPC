/*modified slightly by djs....10/18/04*/

#include<stdio.h>
#include<math.h>
/*#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>*/


/****************************************************************/
/*              nom:            e_nfwg                          */
/*              auteur:         Jean-Paul Kneib                 */
/*              date:              09/04                        */
/*              place:          Toulouse                        */
/****************************************************************/

static double alpha_e;
/*modified declaration to work with my C syntax

double   nfwg_dpl(r, rs,kappas,alpha)

double   r,rs,kappas,alpha;*/
double   nfwg_dpl(double r,double rs,double kappas,double alpha)
{
	
	double scmass(double x);
	double qromo(double (*func)(double), double a, double b);
	/*JP's declaration
	double scmass();
	double qromo();*/
	double x;
	
	x=r/rs;
	alpha_e=alpha;
	
	if ( r !=0.0)
	{
		return ( kappas*2/r*qromo(scmass,0.0,x));
		/* return ( 2*kappas/r*qromo(scmass,0.0,x));
		return ( 2*kappas/r*qromb(scmass,0.0,x)); */
	}
	else
	{
		return ( 0.);
	}

}

/* --------------------------------------------------------------*/

/*JP's declaration
double   nfwg_kappa(r,rs,kappas,alpha)

double   r,rs,kappas,alpha;*/
double nfwg_kappa(double r,double rs, double kappas, double alpha)
{
	double surfdens2(double x, double alpha);
	double x;
	
	x=r/rs;
	
	return( kappas*pow(x,1-alpha)*surfdens2(x,alpha) );
}

/* --------------------------------------------------------------*/
/*****comment this out...not needed
double   nfwg_kappa_eps(r,rs,theta,kappas,eps,alpha)

double   r,rs,theta,kappas,eps,alpha;
{
double kappa_eps;
double nfwg_kappa(),nfwg_gamma();

kappa_eps=nfwg_kappa(r,rs,kappas,alpha)+eps*cos(2*theta)
	  *nfwg_gamma(r,rs,kappas,alpha);

return(kappa_eps);
}
*****/
/* --------------------------------------------------------------*/
/*****comment this out...not needed
double   nfwg_gamma_eps(r,rs,theta,kappas,eps,alpha)

double   r,rs,theta,kappas,eps,alpha;
{
double gamma_eps;
double nfwg_dpl(),nfwg_kappa(),nfwg_gamma(),nfwg_kappa_eps();

gamma_eps=sqrt(pow(nfwg_gamma(r,rs,kappas,alpha),2.)+
2*eps*cos(2*theta)*nfwg_kappa(r,rs,kappas,alpha)*nfwg_gamma(r,rs,kappas,alpha)
+eps*eps*(pow(nfwg_kappa(r,rs,kappas,alpha),2.)-
pow(sin(2*theta)*nfwg_gamma(r,rs,kappas,alpha),2.)));

return(gamma_eps);
}
*/
/* --------------------------------------------------------------*/
/*****comment this out...not needed
double   nfwg_gamma(r,rs,kappas,alpha)

double   r,rs,kappas,alpha;
{
double nfwg_kappa();
double nfwg_dpl();


return( nfwg_dpl(r,rs,kappas,alpha)/r - nfwg_kappa(r,rs,kappas,alpha) );

}
*/
/* --------------------------------------------------------------*/


/* surfdens2 ----------------------------------*/

/*JP's declaration
double surfdens2(x,alpha)
double x,alpha;
*/
double surfdens2(double x, double alpha)
{
	int iter=1001;/*number of Simpson's rule steps*/
	int i=1;/*loop counter*/
	double thetastep;/*step size*/
	double th;/*current value of theta*/
	double sfdens=0.;/*quantity proportional to the projected surface density*/
	
	thetastep = 3.14156/2.0/(iter-1.0);
	th = 0.;
	
	
	sfdens = thetastep / 3.* pow( 1. + x, alpha - 3. ); /*first and last terms
								 +of sum*/
	for( i = 1 ; i <= (iter-1)/2 - 1 ; i++ ) 
	{
		th += thetastep;
		sfdens += thetastep * 4./3. * sin(th) * pow( sin(th)+x, alpha-3. );
		th += thetastep;
		sfdens += thetastep * 2./3. * sin(th) * pow( sin(th)+x, alpha-3. );
	}
		
	th += thetastep;
	sfdens += thetastep * 4./3. * sin(th) * pow( sin(th)+x, alpha-3. );
			    
	/* sfdens = sfdens*2; */
				  
	return sfdens;

}

/*scmass */

/*JP's declaration
double scmass(x)

double x;
*/
double scmass(double x)
{
	double surfdens2(double x, double alpha);
	double scmass_dm;
	
	
	scmass_dm = pow(x,2.0-alpha_e)*surfdens2(x,alpha_e);
	
	
	
	return scmass_dm;

}
