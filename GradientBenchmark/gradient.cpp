#include <iostream>
#include <string.h>
//#include <cuda_runtime.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <immintrin.h>

#include "simd_math.h"
#include "gradient.hpp"
#include "structure.h"
//#include "iacaMarks.h"

//#define __INV RCP
//#define __INV RCP_1NR
#define __INV  RCP_2NR
//#define __SQRT _mm256_sqrt_pd
//#define __SQRT SQRT
//#define __SQRT SQRT_1NR
#define __SQRT SQRT_2NR



//
//
//

/**** usefull functions for PIEMD profile : see old lenstool ****/

/**  I*w,v=0.5 Kassiola & Kovner, 1993 PIEMD, paragraph 4.1
 *
 * Global variables used :
 * - none
 */

/** @brief This module function calculates profile depended information like the impactparameter b0 and the potential ellipticity epot
 * 
 * @param lens: mass distribution for which to calculate parameters
 */
void module_readParameters_calculatePotentialparameter(Potential *lens)
{
        //std::cout << "module_readParameters_calculatePotentialparameter..." << std::endl;
        switch (lens->type)
        {

                case(5): /*Elliptical Isothermal Sphere*/
                        //impact parameter b0
                        lens->b0 = 4* pi_c2 * lens->vdisp * lens->vdisp ;
                        //ellipticity_potential 
                        lens->ellipticity_potential = lens->ellipticity/3 ;
                        break;

                case(8): /* PIEMD */
                        //impact parameter b0
                        lens->b0 = 6.*pi_c2 * lens->vdisp * lens->vdisp;
                        //ellipticity_parameter
                        if ( lens->ellipticity == 0. && lens->ellipticity_potential != 0. ){
                                // emass is (a2-b2)/(a2+b2)
                                lens->ellipticity = 2.*lens->ellipticity_potential / (1. + lens->ellipticity_potential * lens->ellipticity_potential);
                                //printf("1 : %f %f \n",lens->ellipticity,lens->ellipticity_potential);
                        }
                        else if ( lens->ellipticity == 0. && lens->ellipticity_potential == 0. ){
                                lens->ellipticity_potential = 0.00001;
                                //printf("2 : %f %f \n",lens->ellipticity,lens->ellipticity_potential);
                        }
                        else{
                                // epot is (a-b)/(a+b)
                                lens->ellipticity_potential = (1. - sqrt(1 - lens->ellipticity * lens->ellipticity)) / lens->ellipticity;
                                //printf("3 : %f %f \n",lens->ellipticity,lens->ellipticity_potential);
                        }
                        break;

                default:
                        std::cout << "ERROR: LENSPARA profil type of clump "<< lens->name << " unknown : "<< lens->type << std::endl;
                        //printf( "ERROR: LENSPARA profil type of clump %s unknown : %d\n",lens->name, lens->type);
                        break;
        };
}
//
//
//
void module_readParameters_calculatePotentialparameter_SOA(Potential_SOA *lens, int ii)
{
        //std::cout << "module_readParameters_calculatePotentialparameter..." << std::endl;
        switch (lens->type[ii])
        {

                case(5): /*Elliptical Isothermal Sphere*/
                        //impact parameter b0
                        lens->b0[ii] = 4* pi_c2 * lens->vdisp[ii] * lens->vdisp[ii] ;
                        //ellipticity_potential 
                        lens->ellipticity_potential[ii] = lens->ellipticity[ii]/3 ;
                        break;

                case(8): /* PIEMD */
                        //impact parameter b0
                        lens->b0[ii] = 6.*pi_c2 * lens->vdisp[ii] * lens->vdisp[ii];
                        //ellipticity_parameter
                        if ( lens->ellipticity[ii] == 0. && lens->ellipticity_potential[ii] != 0. ){
                                // emass is (a2-b2)/(a2+b2)
                                lens->ellipticity[ii] = 2.*lens->ellipticity_potential[ii] / (1. + lens->ellipticity_potential[ii] * lens->ellipticity_potential[ii]);
                                //printf("1 : %f %f \n",lens->ellipticity[ii],lens->ellipticity_potential[ii]);
                        }
                        else if ( lens->ellipticity[ii] == 0. && lens->ellipticity_potential[ii] == 0. ){
                                lens->ellipticity_potential[ii] = 0.00001;
                                //printf("2 : %f %f \n",lens->ellipticity[ii],lens->ellipticity_potential[ii]);
                        }
                        else{
                                // epot is (a-b)/(a+b)
                                lens->ellipticity_potential[ii] = (1. - sqrt(1 - lens->ellipticity[ii] * lens->ellipticity[ii])) / lens->ellipticity[ii];
                                //printf("3 : %f %f \n",lens->ellipticity[ii],lens->ellipticity_potential[ii]);
                        }
                        break;

                default:
                        std::cout << "ERROR: LENSPARA profil type of clump "<< lens->name << " unknown : "<< lens->type << std::endl;
                        //printf( "ERROR: LENSPARA profil type of clump %s unknown : %d\n",lens->name, lens->type);
                        break;
        };
}


//
/**@brief Return the gradient of the projected lens potential for one clump
 * !!! You have to multiply by dlsds to obtain the true gradient
 * for the expressions, see the papers : JP Kneib & P Natarajan, Cluster Lenses, The Astronomy and Astrophysics Review (2011) for 1 and 2
 * and JP Kneib PhD (1993) for 3
 * 
 * @param pImage 	point where the result is computed in the lens plane
 * @param lens		mass distribution
 */

//
/// Useful functions
//
//static
complex
piemd_1derivatives_ci05(double x, double y, double eps, double rc)
{
        double  sqe, cx1, cxro, cyro, rem2;
        complex zci, znum, zden, zis, zres;
        double norm;
        //
        //std::cout << "piemd_lderivatives" << std::endl;
        sqe  = sqrt(eps);
        cx1  = (1. - eps) / (1. + eps);
        cxro = (1. + eps) * (1. + eps);
        cyro = (1. - eps) * (1. - eps);
        //
        rem2 = x * x / cxro + y * y / cyro;
        /*zci=cpx(0.,-0.5*(1.-eps*eps)/sqe);
          znum=cpx(cx1*x,(2.*sqe*sqrt(rc*rc+rem2)-y/cx1));
          zden=cpx(x,(2.*rc*sqe-y));
          zis=pcpx(zci,lncpx(dcpx(znum,zden)));
          zres=pcpxflt(zis,b0);*/

        // --> optimized code
        zci.re = 0;
        zci.im = -0.5 * (1. - eps * eps) / sqe;
        znum.re = cx1 * x;
        znum.im = 2.*sqe * sqrt(rc * rc + rem2) - y / cx1;
        zden.re = x;
        zden.im = 2.*rc * sqe - y;
        norm = zden.re * zden.re + zden.im * zden.im;     // zis = znum/zden
        zis.re = (znum.re * zden.re + znum.im * zden.im) / norm;
        zis.im = (znum.im * zden.re - znum.re * zden.im) / norm;
        norm = zis.re;
        zis.re = log(sqrt(norm * norm + zis.im * zis.im));  // ln(zis) = ln(|zis|)+i.Arg(zis)
        zis.im = atan2(zis.im, norm);
        //  norm = zis.re;
        zres.re = zci.re * zis.re - zci.im * zis.im;   // Re( zci*ln(zis) )
        zres.im = zci.im * zis.re + zis.im * zci.re;   // Im( zci*ln(zis) )
        //
        //zres.re = zis.re*b0;
        //zres.im = zis.im*b0;
        //
        return(zres);
}
//
//// changes the coordinates of point P into a new basis (rotation of angle theta)
////      y'      y    x'
////       *      |   /
////         *    |  /  theta
////           *  | /
////             *|--------->x
//
struct point rotateCoordinateSystem(struct point P, double theta)
{
	struct  point   Q;

	Q.x = P.x*cos(theta) + P.y*sin(theta);
	Q.y = P.y*cos(theta) - P.x*sin(theta);

	return(Q);
}

__attribute__((noinline)) 
static struct point 
grad_halo(const struct point *pImage, const struct Potential *lens)
{
	struct point true_coord, true_coord_rotation, result;
	double R, angular_deviation;
	complex zis;
	//
	//std::cout << "grad_halo..." << lens->type << std::endl;
	result.x = result.y = 0.;
	//
	/*positionning at the potential center*/
	// Change the origin of the coordinate system to the center of the clump
	true_coord.x = pImage->x - lens->position.x;  
	true_coord.y = pImage->y - lens->position.y;
	//
	switch (lens->type)
	{
		case(5): /*Elliptical Isothermal Sphere*/
			/*rotation of the coordiante axes to match the potential axes*/
			true_coord_rotation = rotateCoordinateSystem(true_coord, lens->ellipticity_angle);
			//
			R=sqrt(true_coord_rotation.x*true_coord_rotation.x*(1 - lens->ellipticity/3.) + true_coord_rotation.y*true_coord_rotation.y*(1 + lens->ellipticity/3.));	//ellippot = ellipmass/3
			result.x = (1 - lens->ellipticity/3.)*lens->b0*true_coord_rotation.x/(R);
			result.y = (1 + lens->ellipticity/3.)*lens->b0*true_coord_rotation.y/(R);
			break;
		case(8): /* PIEMD */
			/*rotation of the coordiante axes to match the potential axes*/
			true_coord_rotation = rotateCoordinateSystem(true_coord, lens->ellipticity_angle);
			/*Doing something....*/
			zis = piemd_1derivatives_ci05(true_coord_rotation.x, true_coord_rotation.y, lens->ellipticity_potential, lens->rcore);
			//
			result.x = lens->b0*zis.re;
			result.y = lens->b0*zis.im;
			break;
		default:
			std::cout << "ERROR: Grad 1 profil type of clump "<< lens->name << " unknown : "<< lens->type << std::endl;
			break;
	};
	return result;
}
//
//
//
struct point module_potentialDerivatives_totalGradient(const runmode_param *runmode, const struct point *pImage, const struct Potential *lens)
{
        struct point grad, clumpgrad;
        grad.x = 0;
        grad.y = 0;
	std::cout << "nhalos = " << runmode->nhalos << std::endl;
        for(int i = 0; i < runmode->nhalos; i++)
        {
                clumpgrad = grad_halo(pImage, &lens[i]);  //compute gradient for each clump separately
		//nan check
		//std::cout << clumpgrad.x << " " << clumpgrad.y << std::endl;
                if(clumpgrad.x == clumpgrad.x or clumpgrad.y == clumpgrad.y)
                { 
			// add the gradients
                        grad.x += clumpgrad.x;
                        grad.y += clumpgrad.y;
                }  
        }
        //
        return(grad);
}

struct point grad_halo_piemd_SOA(const struct point *pImage, int iterator, const struct Potential_SOA *lens)
{
        struct point clumpgrad;
        clumpgrad.x = 0;
        clumpgrad.y = 0;


		struct point true_coord, true_coord_rotation; //, result;
		//double       R, angular_deviation;
		complex      zis;
		//
		//result.x = result.y = 0.;
		//
		true_coord.x = pImage->x - lens->position_x[iterator];
		true_coord.y = pImage->y - lens->position_y[iterator];
		//
		// std::cout << "grad_halo..." << lens->type << std::endl;
		//
		/*positionning at the potential center*/
		// Change the origin of the coordinate system to the center of the clump
		true_coord_rotation = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[iterator]);
		// std::cout << i << ": " << true_coord.x << " " << true_coord.y << " -> "
		//	<< true_coord_rotation.x << " " << true_coord_rotation.y << std::endl;
		//
		//double  sqe, cx1, cxro, cyro, rem2;
		//
		double x   = true_coord_rotation.x;
		double y   = true_coord_rotation.y;
		double eps = lens->ellipticity_potential[iterator];
		double rc  = lens->rcore[iterator];
		//
		//std::cout << "piemd_lderivatives" << std::endl;
		//
		double sqe  = sqrt(eps);
		//
		double cx1  = (1. - eps)/(1. + eps);
		double cxro = (1. + eps)*(1. + eps);
		double cyro = (1. - eps)*(1. - eps);
		//
		double rem2 = x*x/cxro + y*y/cyro;
		//
		/*zci=cpx(0.,-0.5*(1.-eps*eps)/sqe);
		  znum=cpx(cx1*x,(2.*sqe*sqrt(rc*rc+rem2)-y/cx1));
		  zden=cpx(x,(2.*rc*sqe-y));
		  zis=pcpx(zci,lncpx(dcpx(znum,zden)));
		  zres=pcpxflt(zis,b0);*/
		// --> optimized code
		complex zci, znum, zden, zres;
		double norm;
		//	
		zci.re  = 0;
		zci.im  = -0.5*(1. - eps*eps)/sqe;
		//
		znum.re = cx1*x;
		znum.im = 2.*sqe*sqrt(rc*rc + rem2) - y/cx1;
		//
		zden.re = x;
		zden.im = 2.*rc*sqe - y;
		norm    = (zden.re*zden.re + zden.im*zden.im);     // zis = znum/zden
		//
		zis.re  = (znum.re*zden.re + znum.im*zden.im)/norm;
		zis.im  = (znum.im*zden.re - znum.re*zden.im)/norm;
		norm    = zis.re;
		zis.re  = log(sqrt(norm*norm + zis.im*zis.im));  // ln(zis) = ln(|zis|)+i.Arg(zis)
		zis.im  = atan2(zis.im, norm);
		//  norm = zis.re;
		zres.re = zci.re*zis.re - zci.im*zis.im;   // Re( zci*ln(zis) )
		zres.im = zci.im*zis.re + zis.im*zci.re;   // Im( zci*ln(zis) )
		//
		zis.re  = zres.re;
		zis.im  = zres.im;
		//
		//zres.re = zis.re*b0;
		//zres.im = zis.im*b0;
		//
		//
		clumpgrad.x = lens->b0[iterator]*zis.re;
		clumpgrad.y = lens->b0[iterator]*zis.im;
		//nan check
	//
	return(clumpgrad);
}

struct point grad_halo_sis_SOA(const struct point *pImage, int iterator, const struct Potential_SOA *lens)
{
    struct point true_coord, true_coord_rotation, result;
    double R, angular_deviation;
    complex zis;

    result.x = result.y = 0.;

    /*positionning at the potential center*/
    true_coord.x = pImage->x - lens->position_x[iterator];  // Change the origin of the coordinate system to the center of the clump
    true_coord.y = pImage->y - lens->position_y[iterator];

    /*Elliptical Isothermal Sphere*/
	/*rotation of the coordiante axes to match the potential axes*/
	true_coord_rotation = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[iterator]);

	R=sqrt(true_coord_rotation.x*true_coord_rotation.x*(1-lens->ellipticity[iterator]/3.)+true_coord_rotation.y*true_coord_rotation.y*(1+lens->ellipticity[iterator]/3.));	//ellippot = ellipmass/3
	result.x=(1-lens->ellipticity[iterator]/3.)*lens->b0[iterator]*true_coord_rotation.x/(R);
	result.y=(1+lens->ellipticity[iterator]/3.)*lens->b0[iterator]*true_coord_rotation.y/(R);

	return result;
}

struct point module_potentialDerivatives_totalGradient_SOA(const int *Nlens, const struct point *pImage, const struct Potential_SOA *lens)
{
    struct point grad, clumpgrad;
	grad.x=0;
	grad.y=0;

	//This here could be done with function pointer to better acomodate future ass distributions functions
	// However I'm unsure of the time of function pointers -> ask gilles
	//for the moment lens and Nlens is organised the following way :  1. SIS, 2. PIEMD

	//SIS is the first
	for(int i=0; i<Nlens[0]; i++){
		clumpgrad=grad_halo_sis_SOA(pImage,i,&lens[0]);  //compute gradient for each clump separately

		if(clumpgrad.x == clumpgrad.x or clumpgrad.y == clumpgrad.y){ //nan check
		grad.x+=clumpgrad.x;
		grad.y+=clumpgrad.y;
		}  // add the gradients
	}

	//PIEMD is the second
	for(int i=0; i<Nlens[1]; i++){
		clumpgrad=grad_halo_piemd_SOA(pImage,i,&lens[1]);  //compute gradient for each clump separately

		if(clumpgrad.x == clumpgrad.x or clumpgrad.y == clumpgrad.y){ //nan check
		grad.x+=clumpgrad.x;
		grad.y+=clumpgrad.y;
		}  // add the gradients
	}

    return(grad);
}


struct point grad_halo_piemd_SOA_AVX(const int Nlens, const struct point *pImage, const struct Potential_SOA *lens)
{
        struct point grad, clumpgrad;
        grad.x = 0;
        grad.y = 0;
        //std::cout << "optimized gradient: nhalos = " << runmode->nhalos << std::endl;
        //
        // smearing the image coordinates on registers
        __m256d image_x = _mm256_set1_pd(pImage->x);
        __m256d image_y = _mm256_set1_pd(pImage->y);
        //
//#pragma unroll
        for(int i = 0; i < Nlens; i = i + 4)
        {
		//IACA_START;
                //
                //std::cout << "starting loop " << std::endl;
		__m256d two   = _mm256_set1_pd(2.);
		__m256d one   = _mm256_set1_pd(1.);
		__m256d zero  = _mm256_set1_pd(0.);
		__m256d half  = _mm256_set1_pd(0.5);
		__m256d mhalf = _mm256_set1_pd(-0.5);
                //
		// 2 loads
                __m256d true_coord_x      = _mm256_sub_pd(image_x, _mm256_loadu_pd(&lens->position_x[i]));
                __m256d true_coord_y      = _mm256_sub_pd(image_y, _mm256_loadu_pd(&lens->position_y[i]));
		// 2 loafs
                __m256d rc                = _mm256_loadu_pd(&lens->rcore[i]);
		__m256d b0                = _mm256_loadu_pd(&lens->b0[i]);
		// 2 loads
                __m256d eps               = _mm256_loadu_pd(&lens->ellipticity_potential[i]);
		//
		__m256d one_minus_eps     = _mm256_sub_pd(one, eps);
		__m256d one_plus_eps      = _mm256_add_pd(one, eps);
		__m256d one_plus_eps_rcp  = __INV(one_plus_eps);
                // 1 load
                __m256d theta = _mm256_loadu_pd(&lens->ellipticity_angle[i]);
                /*positionning at the potential center*/
		__m256d cos_theta = _mm256_cos_pd(theta);		
		__m256d sin_theta = _mm256_sin_pd(theta);		
		// rotation: 6 ops
		__m256d x = _mm256_add_pd(_mm256_mul_pd(true_coord_x, cos_theta), _mm256_mul_pd(true_coord_y, sin_theta));
		__m256d y = _mm256_sub_pd(_mm256_mul_pd(true_coord_y, cos_theta), _mm256_mul_pd(true_coord_x, sin_theta));
		//
                __m256d sqe   = __SQRT(eps);
                //
		// (1. - eps)/(1. + eps); 3 ops
                __m256d cx1  = one_minus_eps*one_plus_eps_rcp; 
		// (1. + eps)*(1. + eps); 3 ops
                __m256d cxro = one_plus_eps*one_plus_eps; 
		// (1. - eps)*(1. - eps); 3 ops
                __m256d cyro = one_minus_eps*one_minus_eps; 
                //__m256d rem2 = x*x/(cxro) + y*y/(cyro);
                // ~5 ops
                __m256d rem2 = x*x*__INV(cxro) + y*y*__INV(cyro);
                //      
                __m256d zci_re  = zero;
                __m256d zci_im  = mhalf*(one - eps*eps)*__INV(sqe); // ~4 ops
             	//  2.*sqe*sqrt(rc*rc + rem2) - y/cx1 
             	//  7 ops
                __m256d znum_re = cx1*x;
                __m256d znum_im = two*sqe*__SQRT(rc*rc + rem2) - y*__INV(cx1); // ~4 ops
                //
                __m256d zden_re = x;
                __m256d zden_im = _mm256_mul_pd(_mm256_set1_pd(2.), _mm256_mul_pd(rc, sqe)); 
		zden_im = _mm256_sub_pd(zden_im, y);
                //norm    = (zden.re*zden.re + zden.im*zden.im); 3 ops
                __m256d norm    = (zden_re*zden_re + zden_im*zden_im);     
                __m256d zis_re  = (znum_re*zden_re + znum_im*zden_im)*__INV(norm); // 3 ops
                __m256d zis_im  = (znum_im*zden_re - znum_re*zden_im)*__INV(norm); // 3 ops
                norm    	= zis_re;
		//
                zis_re  = _mm256_log_pd(__SQRT(norm*norm + zis_im*zis_im));  // 3 ops// ln(zis) = ln(|zis|)+i.Arg(zis)
                zis_re  = _mm256_log_pd(__SQRT(norm*norm + zis_im*zis_im));  // 3 ops// ln(zis) = ln(|zis|)+i.Arg(zis)
                zis_im  = _mm256_atan2_pd(zis_im, norm); //
		//
                __m256d zres_re = zci_re*zis_re - zci_im*zis_im;   // Re( zci*ln(zis) ) 3 ops
                __m256d zres_im = zci_im*zis_re + zis_im*zci_re;   // Im( zci*ln(zis) ) 3 ops
                //
                zis_re  = zres_re;
                zis_im  = zres_im;
                //
		__m256d b0_x, b0_y;
		//
		b0_x = _mm256_mul_pd(b0, zis_re);    // 1 ops
                clumpgrad.x  = ((double*) &b0_x)[0];
                clumpgrad.x += ((double*) &b0_x)[1];
                clumpgrad.x += ((double*) &b0_x)[2];
                clumpgrad.x += ((double*) &b0_x)[3];
		//
		b0_y = _mm256_mul_pd(b0, zis_im);    // 1 ops
                clumpgrad.y  = ((double*) &b0_y)[0];
                clumpgrad.y += ((double*) &b0_y)[1];
                clumpgrad.y += ((double*) &b0_y)[2];
                clumpgrad.y += ((double*) &b0_y)[3];
		//
                //nan check
                if(clumpgrad.x == clumpgrad.x or clumpgrad.y == clumpgrad.y)
                {
                        // add the gradients
                        grad.x += clumpgrad.x;
                        grad.y += clumpgrad.y;
                }
        }
        for(int i = Nlens - Nlens%4; i < Nlens; ++i){
        	//std::cout << i << std::endl;
        	clumpgrad=grad_halo_piemd_SOA(pImage,i,&lens[1]);  //compute gradient for each clump separately
        	if(clumpgrad.x == clumpgrad.x or clumpgrad.y == clumpgrad.y)
        	{ //nan check
        		grad.x+=clumpgrad.x;
        		grad.y+=clumpgrad.y;
        	}



        }
	//IACA_END;
        //
        return(grad);
}

struct point module_potentialDerivatives_totalGradient_SOA_AVX(const int *Nlens, const struct point *pImage, const struct Potential_SOA *lens)
{
    struct point grad, clumpgrad;
    grad.x = 0;
    grad.y = 0;

	//SIS is the first
	for(int i=0; i<Nlens[0]; i++){
		clumpgrad=grad_halo_sis_SOA(pImage,i,&lens[0]);  //compute gradient for each clump separately

		if(clumpgrad.x == clumpgrad.x or clumpgrad.y == clumpgrad.y){ //nan check
		grad.x+=clumpgrad.x;
		grad.y+=clumpgrad.y;
		}  // add the gradients
	}

	clumpgrad = grad_halo_piemd_SOA_AVX(Nlens[1],pImage,&lens[1]);
	//std::cout << clumpgrad.x << "   " << clumpgrad.y << "   " << Nlens[1] <<std::endl;
	grad.x+=clumpgrad.x;
	grad.y+=clumpgrad.y;

    return grad;

}
