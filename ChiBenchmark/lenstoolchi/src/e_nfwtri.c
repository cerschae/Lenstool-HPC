#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*              nom:            e_nfwtri                          */
/*              auteur:              TV & EJ            */
/*              date:              2012                        */
/*              place:                                       */
/****************************************************************
 *  
 * Adding triaxiality in the NFW profile
 * following  Oguri, Lee & Suto 2003 (OLS03)
 * 
 */

/* --------------------------------------------------------------*
 * Return Eq. 28, A value (OLS03)
 */
double AA_tri(const struct pot *ilens)
{   
	
	double AA_tri1,AA_tri2,AA_tri3,AA_tri4, AA_tri;	
	
	
	AA_tri1 = cos(ilens->theta)*cos(ilens->theta);
	AA_tri2 = (1./ilens->epot)*(1./ilens->epot)*sin(ilens->phi)*sin(ilens->phi);
	AA_tri3 = (1./ilens->emass)*(1./ilens->emass)*cos(ilens->phi)*cos(ilens->phi);
	AA_tri4 = (1./ilens->epot)*(1./ilens->epot)*(1./ilens->emass)*(1./ilens->emass)*sin(ilens->theta)*sin(ilens->theta);
	
	AA_tri =  AA_tri1*(AA_tri2+AA_tri3)  +  AA_tri4;
	
  	
	return (AA_tri);    
}

/* --------------------------------------------------------------*
 * Return Eq. 29, B value (OLS03)
 */
double BB_tri(const struct pot *ilens)
{   
	double BB_tri1,BB_tri2,BB_tri3, BB_tri;	
	
	
	BB_tri1 = cos(ilens->theta)*sin(2.*(ilens->phi));
	BB_tri2 = (1./ilens->epot)*(1./ilens->epot);
	BB_tri3 = (1./ilens->emass)*(1./ilens->emass);
	
	BB_tri = BB_tri1*(BB_tri2 - BB_tri3);
  	
	
	return (BB_tri);    
}

/* --------------------------------------------------------------*
 * Return Eq. 30, C value (OLS03)
 */
double CC_tri(const struct pot *ilens)
{   
	double CC_tri1,CC_tri2, CC_tri;	
	
	
	CC_tri1 = (1./ilens->emass)*(1./ilens->emass)*sin(ilens->phi)*sin(ilens->phi);
	CC_tri2 = (1./ilens->epot)*(1./ilens->epot)*cos(ilens->phi)*cos(ilens->phi);
	
	
	CC_tri = CC_tri1 + CC_tri2;
  	
	
	return (CC_tri);    
}

/* --------------------------------------------------------------*
 * Return Eq. 31, Psi value (OLS03)
 */
double Psi_tri(const struct pot *ilens)
{   
	
    double Psi_tri;

	Psi_tri = (1./2.)*atan(  BB_tri(ilens)/( AA_tri(ilens)-CC_tri(ilens) ) );
  	
	
	return (Psi_tri);  
}
/* --------------------------------------------------------------*
 * Return the ellipticity of the projected ellipsoids, i.e. the ellipticities of the ellipses
 */
double elli_tri(const struct pot *ilens)
{   
	double elli_tri1,elli_tri2, elli_tri;
	
	
	elli_tri1 = 1./( AA_tri(ilens)+CC_tri(ilens)-sqrt((AA_tri(ilens)-CC_tri(ilens))*(AA_tri(ilens)-CC_tri(ilens))+BB_tri(ilens)*BB_tri(ilens)) );
  	elli_tri2 = 1./( AA_tri(ilens)+CC_tri(ilens)+sqrt((AA_tri(ilens)-CC_tri(ilens))*(AA_tri(ilens)-CC_tri(ilens))+BB_tri(ilens)*BB_tri(ilens)) );
	
	elli_tri = (elli_tri1-elli_tri2)/(elli_tri1+elli_tri2);
	
	return (elli_tri);  
}
/* --------------------------------------------------------------*
 * Return the 2D concentration, acording to Jing & Suto 2002
 */

void e_nfw_c3D2c2D(double c3D,double *c2D)
{
	
	*c2D = c3D/0.45;
		
}
