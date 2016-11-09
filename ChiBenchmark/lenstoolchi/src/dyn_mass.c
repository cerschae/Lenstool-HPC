#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<structure.h>
#include<constant.h>
//#include <cmath>
//#include <iostream.h>
#include "lt.h"
#include <gsl/gsl_sf_gamma.h>
#include "gsl/gsl_integration.h"
//#include <gsl/gsl_sf_result.h>

/****************************************************************/
/*		nom:		dyn_mass			*/
/*		auteur:		TV         			*/
/*		date:		2011-2013 			        */
/*		place:					        */
/****************************************************************/
/****************************************************************/
static double FFF(double rr1, double rr2);
//static double MMM(double rr1, double rr2);
static double MMM(double rr1, double rr2, double rr3);
static double MNA(double x, double y,double yy); 
/**************************************************************/

static double Int2(double z1, double z2, double z3);
static double Intz_gsl(double z, void* param1);
static double integral_Intz_ab(double a, double b, double rrr);



static double Int22(double z1, double z2, double z3, double z4);
static double Intz2_gsl(double z, void* param1);
static double integral_Intz2_ab(double a, double b, double rrr, double rrrbetta);


static double Int23(double z1, double z2, double z3, double z4);
static double Intz3_gsl(double z, void* param1);
static double integral_Intz3_ab(double a, double b, double rrr, double rrrbetta);


static double Int24(double z1, double z2, double z3, double z4);
static double Intz4_gsl(double z, void* param1);
static double integral_Intz4_ab(double a, double b, double rrr, double rrrbetta);




struct f_params { double Ref; double Bet; };



double mass2d_NFW(double velocity_disp,double reference_rad, double scale_rad)
/*
 *This function calculate de 2D mass for a NFW profile
*/
{
    double m2d_nfw;
	double velocity_cuad;
	
	reference_rad  /= 1000.;      //in Mpc
	scale_rad   /= 1000.;      //in Mpc
	
	velocity_cuad = velocity_disp*velocity_disp;
	
	
    m2d_nfw = (3.0*PI*INVG)*scale_rad*velocity_cuad*FFF(reference_rad,scale_rad);

    return(m2d_nfw);

}
/**************************************************************
 *See below
 *  
 */
double GG1(double x)
{
	
	double yy;
	double yyy;
	
	yy = 1./( x*x - 1. );
	
	yyy =  yy*(1.  -  sqrt(-1.*yy)*acosh(1./x) );
	
	return(yyy);
	
	
}
/**************************************************************
 * See below
 *  
 */
double GG2(double x)
{
	
	double yy;
	double yyy;
	
	yy = 1./( x*x - 1. );
	
	yyy =  yy*(1.  -  sqrt(yy)*acos(1./x) );
	
	return(yyy);
		
	
}
/**************************************************************
 * Return the function F, Eq. A4, Verdugo et al. 2011
 * There is an error in the paper. In the first factor of both integrands is necesary to eliminate  the square root.
 * The integration using three_eighths is precise enough for these calculations. If we compare with other
 *methods like the Gauss-Kronrod rule our result is precise up to the third decimal. The error in the measured
 *mass is considerable greater than the error propagated by the use of this integration method. In addition,
 *is more easy to implement the three_eights than the Gauus-Kronrond method in the code.
*/
static double FFF(double rr1,double rr2)
{
	double rrr;
	double FFFs;
    rrr = rr1/rr2;
	
	FFFs = three_eighths(0.0,0.9999,100000,GG1) + three_eighths(1.0001,rrr,100000,GG2);
	
	return(FFFs);
}

//
/////
//
double mass3d_NFW(double velocity_disp,double reference_rad, double scale_rad)
/*
 *This function calculate de 3D mass for a NFW profile
 */
{
    double m3d_nfw;
	double velocity_cuad;
	double zz;
	double hhh;
	
	
	reference_rad  /= 1000.;      //in Mpc
	scale_rad   /= 1000.;      //in Mpc
	
	zz = reference_rad/scale_rad;
	
	velocity_cuad = velocity_disp*velocity_disp;
	
	hhh = log(1. + zz) - (zz/(1.+zz));
	
	
	
	
    m3d_nfw = (3.0*PI*INVG)*scale_rad*velocity_cuad*hhh/2.0;
	
    return(m3d_nfw);
	
}
//
///
//////////
//////////////
//////////
/////
///