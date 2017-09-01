#include <iostream>
#include <string.h>
//#include <cuda_runtime.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <immintrin.h>
#include <map>
/*
#ifdef __AVX__
#include "simd_math_avx.h"
#endif
#ifdef __AVX512F__
#include "simd_math_avx512f.h"
#endif
*/
#include "structure_hpc.h"
#include "gradient.hpp"
#include "utils.hpp"
//#include "iacaMarks.h"
//
//
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
        zci.re  = 0;
        zci.im  = -0.5*(1. - eps*eps)/sqe;
        znum.re = cx1 * x;
        znum.im = 2.*sqe*sqrt(rc*rc + rem2) - y/cx1;
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
/*
inline
struct point rotateCoordinateSystem(struct point P, double theta)
{
	struct  point   Q;

	Q.x = P.x*cos(theta) + P.y*sin(theta);
	Q.y = P.y*cos(theta) - P.x*sin(theta);

	return(Q);
}
*/
//
//
struct point
grad_halo(const struct point *pImage, const struct Potential *lens)
{
	struct point true_coord, true_coord_rot, result;
	double X, Y, R, angular_deviation, u;
	complex zis;
	//
	result.x = result.y = 0.;
	//
	/*positionning at the potential center*/
	// Change the origin of the coordinate system to the center of the clump
	true_coord.x   = pImage->x - lens->position.x;  
	true_coord.y   = pImage->y - lens->position.y;
	//printf("x, y = %f, %f\n", lens->position.x, lens->position.y);
	//
	true_coord_rot = rotateCoordinateSystem(true_coord, lens->ellipticity_angle);
	//
	switch (lens->type)
	{
		case(5): /*Elliptical Isothermal Sphere*/
			/*rotation of the coordiante axes to match the potential axes*/
			//true_coord_rotation = rotateCoordinateSystem(true_coord, lens->ellipticity_angle);
			//
			R=sqrt(true_coord_rot.x*true_coord_rot.x*(1 - lens->ellipticity_potential) + true_coord_rot.y*true_coord_rot.y*(1 + lens->ellipticity_potential));	//ellippot = ellipmass/3
			result.x = (1 - lens->ellipticity_potential)*lens->b0*true_coord_rot.x/(R);
			result.y = (1 + lens->ellipticity_potential)*lens->b0*true_coord_rot.y/(R);
			break;
		case(8): /* PIEMD */
			/*rotation of the coordiante axes to match the potential axes*/
			/*Doing something....*/
			complex zis = piemd_1derivatives_ci05(true_coord_rot.x, true_coord_rot.y, lens->ellipticity_potential, lens->rcore);
			//
			result.x = lens->b0*zis.re;
			result.y = lens->b0*zis.im;
			break;
		case(81): //PIEMD Kassiola & Kovner,1993 I0.5c-I0.5cut
			double t05;
			if ( lens->ellipticity_potential > 2E-4 )
			{
				//printf("1 ");
				t05 = lens->b0*lens->rcut/(lens->rcut - lens->rcore);
				complex zis     = piemd_1derivatives_ci05(true_coord_rot.x, true_coord_rot.y, lens->ellipticity_potential, lens->rcore);
				complex zis_cut = piemd_1derivatives_ci05(true_coord_rot.x, true_coord_rot.y, lens->ellipticity_potential, lens->rcut);
				result.x = t05 * (zis.re - zis_cut.re);
				result.y = t05 * (zis.im - zis_cut.im);
				//printf("	g   = %f %f\n", result.x, result.y);
			}
			else if (( u = true_coord_rot.x*true_coord_rot.x + true_coord_rot.y*true_coord_rot.y) > 0. )
			{
				//printf("2 ");
				// Circular dPIE Elliasdottir 2007 Eq A23 slighly modified for t05
				X = lens->rcore;
				Y = lens->rcut;
				t05  = sqrt(u + X * X) - X - sqrt(u + Y * Y) + Y;  // Faster and equiv to Elliasdottir (see Golse PhD)
				t05 *= lens->b0 * Y / (Y - X) / u; // 1/u because t05/sqrt(u) and normalised Q/sqrt(u)
				result.x = t05*true_coord_rot.x;
				result.y = t05*true_coord_rot.y;
			}
			else
			{
				//printf("3 ");
				result.x = 0.;
				result.y = 0.;
			}
			break;


		default:
			std::cout << "ERROR: Grad 1 profil type of clump "<< lens->name << " unknown : "<< lens->type << std::endl;
			break;
	};
	result = rotateCoordinateSystem(result, -lens->ellipticity_angle); 
	//printf("	rot grad = %f %f\n", result.x, result.y); 
	return result;
}
//
//
//
struct point module_potentialDerivatives_totalGradient(const int nhalos, const struct point* pImage, const struct Potential* lens)
{
	struct point grad, clumpgrad;
	//
	grad.x = 0;
	grad.y = 0;
	//
	for(int i = 0; i < nhalos; i++)
	{
		clumpgrad = grad_halo(pImage, &lens[i]);  //compute gradient for each clump separately
		//std::cout << clumpgrad.x << " " << clumpgrad.y << std::endl;
		//nan check
		//if(clumpgrad.x == clumpgrad.x or clumpgrad.y == clumpgrad.y)
		{ 
			// add the gradients
			grad.x += clumpgrad.x;
			grad.y += clumpgrad.y;
		}  
	}
	//	
	return(grad);
}
//
// SOA versions, vectorizable
//
struct point module_potentialDerivatives_totalGradient_5_SOA(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos)
{
        asm volatile("# module_potentialDerivatives_totalGradient_SIS_SOA begins");
        //printf("# module_potentialDerivatives_totalGradient_SIS_SOA begins\n");
	//
	struct point grad, result;
        grad.x = 0;
        grad.y = 0;
	for(int i = shalos; i < shalos + nhalos; i++)
	{
		//
		struct point true_coord, true_coord_rotation;
		//
		true_coord.x = pImage->x - lens->position_x[i];
                true_coord.y = pImage->y - lens->position_y[i];
		//
		true_coord_rotation = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[i]);
		//
		double ell_pot = lens->ellipticity_potential[i];
                //
		double R = sqrt(true_coord_rotation.x*true_coord_rotation.x*(1 - ell_pot) + true_coord_rotation.y*true_coord_rotation.y*(1 + ell_pot));
		//
		result.x = (1 - ell_pot)*lens->b0[i]*true_coord_rotation.x/R;
		result.y = (1 + ell_pot)*lens->b0[i]*true_coord_rotation.y/R;
		//
		result = rotateCoordinateSystem(result, -lens->ellipticity_angle[i]);
		//
		grad.x += result.x;
		grad.y += result.y;
		
	}
	return grad;
}
//
//
//
struct point module_potentialDerivatives_totalGradient_5_SOA_v2(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos)
{
        asm volatile("# module_potentialDerivatives_totalGradient_SIS_SOA_v2 begins");
        //printf("# module_potentialDerivatives_totalGradient_SIS_SOA_v2 begins\n");
        //
        struct point grad, result;
        grad.x = 0;
        grad.y = 0;
        for(int i = shalos; i < shalos + nhalos; i++)
        {
                //
                struct point true_coord, true_coord_rotation;
                //
                true_coord.x = pImage->x - lens->position_x[i];
                true_coord.y = pImage->y - lens->position_y[i];
                //
                //true_coord_rotation = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[i]);
		double cose = lens->anglecos[i];
                double sine = lens->anglesin[i];
		//
		double x = true_coord.x*cose + true_coord.y*sine;
                double y = true_coord.y*cose - true_coord.x*sine;
		//
		double ell_pot = lens->ellipticity_potential[i];
		//
		double R = sqrt(x*x*(1 - ell_pot) + y*y*(1 + ell_pot));
		//
		result.x = (1 - ell_pot)*lens->b0[i]*x/R;
                result.y = (1 + ell_pot)*lens->b0[i]*y/R;
                //
                grad.x += result.x*cose - result.y*sine; 
                grad.y += result.y*cose + result.x*sine; 
        }
        return grad;
}
//
//
//
struct point module_potentialDerivatives_totalGradient_5_SOA_v2_novec(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos)
{
        asm volatile("# module_potentialDerivatives_totalGradient_SIS_SOA begins");
        //
        struct point grad, result;
        grad.x = 0;
        grad.y = 0;
#pragma novec
	for(int i = shalos; i < shalos + nhalos; i++)
        {
                //
                struct point true_coord, true_coord_rotation;
                //
                true_coord.x = pImage->x - lens->position_x[i];
                true_coord.y = pImage->y - lens->position_y[i];
                //
                //true_coord_rotation = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[i]);
                double cose = lens->anglecos[i];
                double sine = lens->anglesin[i];
                //
                double x = true_coord.x*cose + true_coord.y*sine;
                double y = true_coord.y*cose - true_coord.x*sine;
                //:
                double R = sqrt(x*x*(1 - lens->ellipticity[i]/3.) + y*y*(1 + lens->ellipticity[i]/3.));
                result.x = (1 - lens->ellipticity[i]/3.)*lens->b0[i]*x/R;
                result.y = (1 + lens->ellipticity[i]/3.)*lens->b0[i]*y/R;
                //
                grad.x += result.x*cose - result.y*sine;
                grad.y += result.y*cose + result.x*sine;
        }
        return grad;
}
//
//
//
struct point module_potentialDerivatives_totalGradient_8_SOA(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos)
{
	asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
	// 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
	// 
	struct point grad, clumpgrad;
	grad.x = 0;
	grad.y = 0;
	//
	for(int i = shalos; i < shalos + nhalos; i++)
	{
		//IACA_START;
		//
		struct point true_coord, true_coord_rot; //, result;
		//double       R, angular_deviation;
		complex      zis;
		//
		//result.x = result.y = 0.;
		//
		//@@printf("image_x = %f image_y = %f\n",  pImage->x, pImage->y);	
		true_coord.x = pImage->x - lens->position_x[i];
		true_coord.y = pImage->y - lens->position_y[i];
		//printf("x = %f y = %f\n",  true_coord.x, true_coord.y);	
		/*positionning at the potential center*/
		// Change the origin of the coordinate system to the center of the clump
		true_coord_rot = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[i]);
		//
		double x   = true_coord_rot.x;
		double y   = true_coord_rot.y;
		//@@printf("x = %f y = %f\n",  x, y);	
		double eps = lens->ellipticity_potential[i];
		double rc  = lens->rcore[i];
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
		complex zci, znum, zden, zres;
		double norm;
		//	
		zci.re  = 0;
		zci.im  = -0.5*(1. - eps*eps)/sqe;
		//@@printf("zci = %f %f\n", zci.re, zci.im);	
		//
		znum.re = cx1*x;
		znum.im = 2.*sqe*sqrt(rc*rc + rem2) - y/cx1;
		//
		zden.re = x;
		zden.im = 2.*rc*sqe - y;
		norm    = (zden.re*zden.re + zden.im*zden.im);     // zis = znum/zden
		//@@printf("norm = %f\n", norm);
		//
		zis.re  = (znum.re*zden.re + znum.im*zden.im)/norm;
		zis.im  = (znum.im*zden.re - znum.re*zden.im)/norm;
		//@@printf("zis = %f %f\n", zis.re, zis.im);	
		norm    = zis.re;
		zis.re  = log(sqrt(norm*norm + zis.im*zis.im));  // ln(zis) = ln(|zis|)+i.Arg(zis)
		zis.im  = atan2(zis.im, norm);
		//@@printf("y,x = %f %f\n", zis.im, norm);	
		//  norm = zis.re;
		zres.re = zci.re*zis.re - zci.im*zis.im;   // Re( zci*ln(zis) )
		zres.im = zci.im*zis.re + zis.im*zci.re;   // Im( zci*ln(zis) )
		//
		//@@printf("zres: %f %f\n", zres.re, zres.im); 
		//
		zis.re  = zres.re;
		zis.im  = zres.im;
		//
		//zres.re = zis.re*b0;
		//zres.im = zis.im*b0;
		// rotation
		clumpgrad.x = zis.re;
		clumpgrad.y = zis.im;
		clumpgrad = rotateCoordinateSystem(clumpgrad, -lens->ellipticity_angle[i]);
		//
		clumpgrad.x = lens->b0[i]*clumpgrad.x;
		clumpgrad.y = lens->b0[i]*clumpgrad.y;
		//
		//clumpgrad.x = lens->b0[i]*zis.re;
		//clumpgrad.y = lens->b0[i]*zis.im;
		//nan check
		//if(clumpgrad.x == clumpgrad.x or clumpgrad.y == clumpgrad.y)
		//{
			// add the gradients
		grad.x += clumpgrad.x;
		grad.y += clumpgrad.y;
		//@@printf("grad: %f %f\n", grad.x, grad.y); 
		//@@std::cout << "grad.x = " << grad.x << " grad.y = " << grad.y << std::endl;
		//}
	}
	//IACA_END;
	//
	return(grad);
}
//
//
//
struct point module_potentialDerivatives_totalGradient_8_SOA_v2(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos)
{
        asm volatile("# module_potentialDerivatives_totalGradient_8_SOA_v2 begins");
        //std::cout << "# module_potentialDerivatives_totalGradient_8_SOA_v2 begins" << std::endl;
	//
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        //
        struct point grad, clumpgrad;
        grad.x = 0;
        grad.y = 0;
	//printf("%d %d\n", shalos, nhalos);
        //
        for(int i = shalos; i < shalos + nhalos; i++)
        {
                //IACA_START;
                //
                struct point true_coord, true_coord_rot; //, result;
                complex      zis;
                double b0   = lens->b0[i];
                //
                true_coord.x = pImage->x - lens->position_x[i];
                true_coord.y = pImage->y - lens->position_y[i];
                //
                //true_coord_rot = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[i]);
                //
                double cose = lens->anglecos[i];
                double sine = lens->anglesin[i];
		//
                double x = true_coord.x*cose + true_coord.y*sine;
                double y = true_coord.y*cose - true_coord.x*sine;
                //
                double eps = lens->ellipticity_potential[i];
                double rc  = lens->rcore[i];
                //
                double sqe  = sqrt(eps);
                //
                double cx1  = (1. - eps)/(1. + eps);
                double cxro = (1. + eps)*(1. + eps);
                double cyro = (1. - eps)*(1. - eps);
                //
                double rem2 = x*x/cxro + y*y/cyro;
                //
                complex zci, znum, zden, zres;
                double norm;
                //
                zci.re  = 0;
                zci.im  = -0.5*(1. - eps*eps)/sqe;
                //
		KERNEL(rc, zres)
		//
                grad.x += b0*(zres.re*cose - zres.im*sine);
                grad.y += b0*(zres.im*cose + zres.re*sine);
                //
        }
        //IACA_END;
        //
        return(grad);
}
//
//
//
struct point module_potentialDerivatives_totalGradient_8_SOA_v2_novec(const struct point *pImage, const struct Potential_SOA *lens
, int shalos, int nhalos)
{
        asm volatile("# module_potentialDerivatives_totalGradient_8_SOA_v2 begins");
        //std::cout << "# module_potentialDerivatives_totalGradient_8_SOA_v2 begins" << std::endl;
        //
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        //
        struct point grad, clumpgrad;
        grad.x = 0;
        grad.y = 0;
        //printf("%d %d\n", shalos, nhalos);
        //
#pragma novector
        for(int i = shalos; i < shalos + nhalos; i++)
        {
                //IACA_START;
                //
                struct point true_coord, true_coord_rot; //, result;
                complex      zis;
                double b0   = lens->b0[i];
                //
                true_coord.x = pImage->x - lens->position_x[i];
                true_coord.y = pImage->y - lens->position_y[i];
                //
                //true_coord_rot = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[i]);
                //
                double cose = lens->anglecos[i];
                double sine = lens->anglesin[i];
                //
                double x = true_coord.x*cose + true_coord.y*sine;
                double y = true_coord.y*cose - true_coord.x*sine;
                //
                double eps = lens->ellipticity_potential[i];
                double rc  = lens->rcore[i];
                //
                double sqe  = sqrt(eps);
                //
                double cx1  = (1. - eps)/(1. + eps);
                double cxro = (1. + eps)*(1. + eps);
                double cyro = (1. - eps)*(1. - eps);
                //
                double rem2 = x*x/cxro + y*y/cyro;
                //
                complex zci, znum, zden, zres;
                double norm;
                //
                zci.re  = 0;
                zci.im  = -0.5*(1. - eps*eps)/sqe;
                //
                KERNEL(rc, zres)
                //
                grad.x += b0*(zres.re*cose - zres.im*sine);
                grad.y += b0*(zres.im*cose + zres.re*sine);
                //
        }
        //IACA_END;
        //
        return(grad);
}
//
//
//
struct point module_potentialDerivatives_totalGradient_81_SOA(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos)
{
        asm volatile("# module_potentialDerivatives_totalGradient_81_SOA begins");
	//std::cout << "# module_potentialDerivatives_totalGradient_SOA begins" << std::endl;
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        // 
        struct point grad, clumpgrad;
        grad.x = 0;
        grad.y = 0;
        for(int i = shalos; i < shalos + nhalos; i++)
        {
                //IACA_START;
                //
                struct point true_coord, true_coord_rot; //, result;
                //double       R, angular_deviation;
                complex      zis;
                //
                //result.x = result.y = 0.;
                //
                true_coord.x = pImage->x - lens->position_x[i];
                true_coord.y = pImage->y - lens->position_y[i];
                /*positionning at the potential center*/
                // Change the origin of the coordinate system to the center of the clump
                true_coord_rot = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[i]);
                //
                double x    = true_coord_rot.x;
                double y    = true_coord_rot.y;
                double eps  = lens->ellipticity_potential[i];
                double rc   = lens->rcore[i];
                double rcut = lens->rcut[i];
		double b0   = lens->b0[i];
		double t05  = b0*rcut/(rcut - rc);
                //
                double sqe  = sqrt(eps);
                //
                double cx1  = (1. - eps)/(1. + eps);
                double cxro = (1. + eps)*(1. + eps);
                double cyro = (1. - eps)*(1. - eps);
                //
                double rem2 = x*x/cxro + y*y/cyro;
                //
                complex zci, znum, zden, zres_rc, zres_rcut;
                double norm;
                //      
                zci.re  = 0;
                zci.im  = -0.5*(1. - eps*eps)/sqe;
		// step 1
		{
			KERNEL(rc, zres_rc)
		}
		// step 2
		{
			KERNEL(rcut, zres_rcut)
                }
		zis.re  = t05*(zres_rc.re - zres_rcut.re);
		zis.im  = t05*(zres_rc.im - zres_rcut.im); 
                // rotation
                clumpgrad.x = zis.re;
                clumpgrad.y = zis.im;
                clumpgrad = rotateCoordinateSystem(clumpgrad, -lens->ellipticity_angle[i]);
                //
                grad.x += clumpgrad.x;
                grad.y += clumpgrad.y;
                //}
        }
        //IACA_END;
        //
        return(grad);
}
//
//
//
struct point module_potentialDerivatives_totalGradient_81_SOA_v2(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos)
{
        asm volatile("# module_potentialDerivatives_totalGradient_81_SOA begins");
        //std::cout << "# module_potentialDerivatives_totalGradient_81_SOA begins" << std::endl;
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        //
        struct point grad, clumpgrad;
        grad.x = 0;
        grad.y = 0;
        for(int i = shalos; i < shalos + nhalos; i++)
        {
                //IACA_START;
                //
                struct point true_coord, true_coord_rot; //, result;
                //double       R, angular_deviation;
                complex      zis;
                //
                //result.x = result.y = 0.;
                //
                true_coord.x = pImage->x - lens->position_x[i];
                true_coord.y = pImage->y - lens->position_y[i];
                /*positionning at the potential center*/
                // Change the origin of the coordinate system to the center of the clump
                double cose = lens->anglecos[i];
                double sine = lens->anglesin[i];
                double x = true_coord.x*cose + true_coord.y*sine;
                double y = true_coord.y*cose - true_coord.x*sine;
		//
                double eps  = lens->ellipticity_potential[i];
                double rc   = lens->rcore[i];
                double rcut = lens->rcut[i];
                double b0   = lens->b0[i];
                double t05  = b0*rcut/(rcut - rc);
                //
                double sqe  = sqrt(eps);
                //
                double cx1  = (1. - eps)/(1. + eps);
                double cxro = (1. + eps)*(1. + eps);
                double cyro = (1. - eps)*(1. - eps);
                //
                double rem2 = x*x/cxro + y*y/cyro;
                //
                complex zci, znum, zden, zres_rc, zres_rcut;
                double norm;
                //
                zci.re  = 0;
                zci.im  = -0.5*(1. - eps*eps)/sqe;
		//
                // step 1
                {
			KERNEL(rc, zres_rc)
                }
                // step 2
                {
			KERNEL(rcut, zres_rcut)
                }
                zis.re  = t05*(zres_rc.re - zres_rcut.re);
                zis.im  = t05*(zres_rc.im - zres_rcut.im);
                // rotation
		grad.x += (zis.re*cose - zis.im*sine);
                grad.y += (zis.im*cose + zis.re*sine);
                //
        }
        //
        return(grad);
}
//
//
//
struct point module_potentialDerivatives_totalGradient_81_SOA_v2_novec(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos)
{
        asm volatile("# module_potentialDerivatives_totalGradient_81_SOA begins");
        //std::cout << "# module_potentialDerivatives_totalGradient_81_SOA begins" << std::endl;
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        struct point grad, clumpgrad;
        grad.x = 0;
        grad.y = 0;
#pragma novector
        for(int i = shalos; i < shalos + nhalos; i++)
        {
                //IACA_START;
                //
                struct point true_coord, true_coord_rot; //, result;
                //double       R, angular_deviation;
                complex      zis;
                //
                //result.x = result.y = 0.;
                //
                true_coord.x = pImage->x - lens->position_x[i];
                true_coord.y = pImage->y - lens->position_y[i];
                /*positionning at the potential center*/
                // Change the origin of the coordinate system to the center of the clump
                double cose = lens->anglecos[i];
                double sine = lens->anglesin[i];
                double x = true_coord.x*cose + true_coord.y*sine;
                double y = true_coord.y*cose - true_coord.x*sine;
                //
                double eps  = lens->ellipticity_potential[i];
                double rc   = lens->rcore[i];
                double rcut = lens->rcut[i];
                double b0   = lens->b0[i];
                double t05  = b0*rcut/(rcut - rc);
                //
                double sqe  = sqrt(eps);
                //
                double cx1  = (1. - eps)/(1. + eps);
                double cxro = (1. + eps)*(1. + eps);
                double cyro = (1. - eps)*(1. - eps);
                //
                double rem2 = x*x/cxro + y*y/cyro;
                //
                complex zci, znum, zden, zres_rc, zres_rcut;
                double norm;
                //
                zci.re  = 0;
                zci.im  = -0.5*(1. - eps*eps)/sqe;
                //
                // step 1
                {
                        KERNEL(rc, zres_rc)
                }
                // step 2
                {
                        KERNEL(rcut, zres_rcut)
                }
                zis.re  = t05*(zres_rc.re - zres_rcut.re);
                zis.im  = t05*(zres_rc.im - zres_rcut.im);
                // rotation
                grad.x += (zis.re*cose - zis.im*sine);
                grad.y += (zis.im*cose + zis.re*sine);
                //
        }
        //
        return(grad);
}
//
//
//
typedef struct point (*halo_func_t) (const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos); 
halo_func_t halo_func[100] = 
{ 
0, 0, 0, 0, 0, module_potentialDerivatives_totalGradient_5_SOA_v2, 0, 0, module_potentialDerivatives_totalGradient_8_SOA_v2,  0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0,  module_potentialDerivatives_totalGradient_81_SOA_v2, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};
//
//
//
struct point module_potentialDerivatives_totalGradient_SOA(const struct point *pImage, const struct Potential_SOA *lens, int nhalos)
{
        struct point grad, clumpgrad;
        //
        grad.x = clumpgrad.x = 0;
        grad.y = clumpgrad.y = 0;
	//
	int shalos = 0;
	//
	//module_potentialDerivatives_totalGradient_81_SOA(pImage, lens, 0, nhalos);
	//return;
	/*
	int* p_type = &(lens->type)[0];
	int* lens_type = (int*) malloc(nhalos*sizeof(int));
	memcpy(lens_type, &(lens->type)[0], nhalos*sizeof(int));
	*/
	//quicksort(lens_type, nhalos);
	//
	while (shalos < nhalos)
	{
		int lens_type = lens->type[shalos];
		int count     = 1;
		while (lens->type[shalos + count] == lens_type) count++;
		//std::cerr << "type = " << lens_type << " " << count << " " << shalos << " " << " func = " << halo_func[lens_type] << std::endl;
		//	
		clumpgrad = (*halo_func[lens_type])(pImage, lens, shalos, count);
		//
		grad.x += clumpgrad.x;
		grad.y += clumpgrad.y;
		shalos += count;
	}	

        return(grad);
}



typedef struct point (*halo_func_t_novec) (const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);
halo_func_t_novec halo_func_novec[100] =
{
0, 0, 0, 0, 0, module_potentialDerivatives_totalGradient_5_SOA_v2_novec, 0, 0, module_potentialDerivatives_totalGradient_8_SOA_v2_novec,  0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0,  module_potentialDerivatives_totalGradient_81_SOA_v2_novec, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};
//
//
//
struct point module_potentialDerivatives_totalGradient_SOA_novec(const struct point *pImage, const struct Potential_SOA *lens, int
 nhalos)
{
        struct point grad, clumpgrad;
        //
        grad.x = clumpgrad.x = 0;
        grad.y = clumpgrad.y = 0;
        //
        int shalos = 0;
        //
        //module_potentialDerivatives_totalGradient_81_SOA(pImage, lens, 0, nhalos);
        //return;
        /*
        int* p_type = &(lens->type)[0];
        int* lens_type = (int*) malloc(nhalos*sizeof(int));
        memcpy(lens_type, &(lens->type)[0], nhalos*sizeof(int));
        */
        //quicksort(lens_type, nhalos);
        //
        while (shalos < nhalos)
        {
                int lens_type = lens->type[shalos];
                int count     = 1;
                while (lens->type[shalos + count] == lens_type) count++;
                //std::cerr << "type = " << lens_type << " " << count << " " << shalos << std::endl;
                //
                clumpgrad = (*halo_func_novec[lens_type])(pImage, lens, shalos, count);
                //
                grad.x += clumpgrad.x;
                grad.y += clumpgrad.y;
                shalos += count;
        }

        return(grad);
}



