/**
Lenstool-HPC: HPC based massmodeling software and Lens-map generation
Copyright (C) 2017  Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

@brief: Handvectorized gradient function for AVX2

*/

#include <iostream>
#include <string.h>
//#include <cuda_runtime.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <immintrin.h>
//
#include "structure_hpc.hpp"
#include "simd_math_avx.h"
#include "gradient.hpp"
//
//#include "iacaMarks.h"
//
//
//
void print256(__m256d __reg)
{
	printf("%f\n", ((double*) &__reg)[0]);
	printf("%f\n", ((double*) &__reg)[1]);
	printf("%f\n", ((double*) &__reg)[2]);
	printf("%f\n", ((double*) &__reg)[3]);
}
//
//
//

//
//
//
#ifdef _double

inline
static 
void rotateCoordinateSystem_avx(__m256d Q_x, __m256d Q_y, const __m256d P_x, const __m256d P_y, double* angles)
{
	//
	__m256d theta = _mm256_loadu_pd(angles);
	/*positionning at the potential center*/
	__m256d cos_theta = _mm256_cos_pd(theta);
	__m256d sin_theta = _mm256_sin_pd(theta);
	// rotation: 6 ops
	Q_x = P_x*cos_theta + P_y*sin_theta;
	Q_y = P_y*cos_theta - P_x*sin_theta;
}
//
//
//
struct point module_potentialDerivatives_totalGradient_5_SOA_AVX(const struct point *pImage, const struct Potential_SOA *lens, int
shalos, int nhalos)
{
        asm volatile("# module_potentialDerivatives_totalGradient_SIS_SOA begins");
        //
        struct point grad, clumpgrad;
        grad.x = 0;
        grad.y = 0;
        //
        // smearing the image coordinates on registers
        __m256d image_x  = _mm256_set1_pd(pImage->x);
        __m256d image_y  = _mm256_set1_pd(pImage->y);
        //
        __m256d __grad_x = _mm256_set1_pd(0.);
        __m256d __grad_y = _mm256_set1_pd(0.);
        //
        __m256d __result_x = _mm256_set1_pd(0.);
        __m256d __result_y = _mm256_set1_pd(0.);
        //
        __m256d one      = _mm256_set1_pd( 1. );
        //
        int i;
        int imax = shalos + nhalos;
//#pragma unroll
        for(i = shalos; i < imax - imax%4; i = i + 4)
        {
                //
                struct point true_coord, true_coord_rotation;
                //
                __m256d b0                = _mm256_loadu_pd(&lens->b0[i]);
                //
                //true_coord.x = pImage->x - lens->position_x[i];
                //true_coord.y = pImage->y - lens->position_y[i];
                //
                __m256d true_coord_x      = _mm256_sub_pd(image_x, _mm256_loadu_pd(&lens->position_x[i]));
                __m256d true_coord_y      = _mm256_sub_pd(image_y, _mm256_loadu_pd(&lens->position_y[i]));
                //
                __m256d theta     = _mm256_loadu_pd(&lens->ellipticity_angle[i]);
                /*positionning at the potential center*/
                __m256d cos_theta;
                __m256d sin_theta = _mm256_sincos_pd(&cos_theta, theta);
                //
                __m256d x         = _mm256_add_pd(_mm256_mul_pd(true_coord_x, cos_theta), _mm256_mul_pd(true_coord_y, sin_theta));
                __m256d y         = _mm256_sub_pd(_mm256_mul_pd(true_coord_y, cos_theta), _mm256_mul_pd(true_coord_x, sin_theta));
                //
                //true_coord_rotation = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[i]);
                //
                __m256d ell_pot               = _mm256_loadu_pd(&lens->ellipticity_potential[i]);
                //
                __m256d _R  = __SQRT(ADD(MUL(x, MUL(x, SUB(one, ell_pot))), MUL(y, MUL(y, ADD(one, ell_pot)))));
                //
                __result_x = MUL(MUL(SUB(one, ell_pot), b0), MUL(x, __INV(_R)));
                __result_y = MUL(MUL(ADD(one, ell_pot), b0), MUL(y, __INV(_R)));
                //
                __grad_x = __grad_x + _mm256_sub_pd(_mm256_mul_pd(__result_x, cos_theta), _mm256_mul_pd(__result_y, sin_theta));
                __grad_y = __grad_y + _mm256_add_pd(_mm256_mul_pd(__result_y, cos_theta), _mm256_mul_pd(__result_x, sin_theta));
                //
        }
        //
        grad.x  = ((double*) &__grad_x)[0];
        grad.x += ((double*) &__grad_x)[1];
        grad.x += ((double*) &__grad_x)[2];
        grad.x += ((double*) &__grad_x)[3];
        //
        grad.y  = ((double*) &__grad_y)[0];
        grad.y += ((double*) &__grad_y)[1];
        grad.y += ((double*) &__grad_y)[2];
        grad.y += ((double*) &__grad_y)[3];
        //
        // end of peeling
        //
        if (nhalos%4 > 0)
        {
                struct point grad_peel;
                grad_peel = module_potentialDerivatives_totalGradient_5_SOA(pImage, lens, i, nhalos%4);
                //
                grad.x += grad_peel.x;
                grad.y += grad_peel.y;
        }
        //


        return grad;
}
//
//
//
struct point module_potentialDerivatives_totalGradient_5_SOA_AVX_v2(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos)
{
        asm volatile("# module_potentialDerivatives_totalGradient_SIS_SOA begins");
        //
        struct point grad, clumpgrad;
        grad.x = 0;
        grad.y = 0;
	//
        // smearing the image coordinates on registers
        __m256d image_x  = _mm256_set1_pd(pImage->x);
        __m256d image_y  = _mm256_set1_pd(pImage->y);
        //
        __m256d __grad_x = _mm256_set1_pd(0.);
        __m256d __grad_y = _mm256_set1_pd(0.);
	//
	__m256d __result_x = _mm256_set1_pd(0.);
        __m256d __result_y = _mm256_set1_pd(0.);
	//
	__m256d one      = _mm256_set1_pd( 1. );
	//
	int i;
        int imax = shalos + nhalos;
//#pragma unroll
        for(i = shalos; i < imax - imax%4; i = i + 4)
        {
                //
                struct point true_coord, true_coord_rotation;
		//
		__m256d b0        = _mm256_loadu_pd(&lens->b0[i]);
		__m256d theta     = _mm256_loadu_pd(&lens->ellipticity_angle[i]);
		__m256d ell_pot   = _mm256_loadu_pd(&lens->ellipticity_potential[i]);
                //
                //true_coord.x = pImage->x - lens->position_x[i];
                //true_coord.y = pImage->y - lens->position_y[i];
		//
		__m256d true_coord_x      = _mm256_sub_pd(image_x, _mm256_loadu_pd(&lens->position_x[i]));
                __m256d true_coord_y      = _mm256_sub_pd(image_y, _mm256_loadu_pd(&lens->position_y[i]));
                //
                /*positionning at the potential center*/
		__m256d cos_theta = _mm256_loadu_pd(&lens->anglecos[i]);
                __m256d sin_theta = _mm256_loadu_pd(&lens->anglesin[i]);
                //
		__m256d x         = _mm256_add_pd(_mm256_mul_pd(true_coord_x, cos_theta), _mm256_mul_pd(true_coord_y, sin_theta));
                __m256d y         = _mm256_sub_pd(_mm256_mul_pd(true_coord_y, cos_theta), _mm256_mul_pd(true_coord_x, sin_theta));
		//
                //true_coord_rotation = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[i]);
		//
		//	
                __m256d _R  = __INV(__SQRT(ADD(MUL(x, MUL(x, SUB(one, ell_pot))), MUL(y, MUL(y, ADD(one, ell_pot))))));
		//
                __result_x = MUL(MUL(SUB(one, ell_pot), b0), MUL(x, _R));
                __result_y = MUL(MUL(ADD(one, ell_pot), b0), MUL(y, _R));
		//
		__grad_x = __grad_x + _mm256_sub_pd(_mm256_mul_pd(__result_x, cos_theta), _mm256_mul_pd(__result_y, sin_theta));
                __grad_y = __grad_y + _mm256_add_pd(_mm256_mul_pd(__result_y, cos_theta), _mm256_mul_pd(__result_x, sin_theta));
		//
        }
	//
	grad.x  = ((double*) &__grad_x)[0];
        grad.x += ((double*) &__grad_x)[1];
        grad.x += ((double*) &__grad_x)[2];
        grad.x += ((double*) &__grad_x)[3];
        //
        grad.y  = ((double*) &__grad_y)[0];
        grad.y += ((double*) &__grad_y)[1];
        grad.y += ((double*) &__grad_y)[2];
        grad.y += ((double*) &__grad_y)[3];
        //
        // end of peeling
        //
        if (nhalos%4 > 0)
        {
                struct point grad_peel;
                grad_peel = module_potentialDerivatives_totalGradient_5_SOA(pImage, lens, i, nhalos%4);
                //
                grad.x += grad_peel.x;
                grad.y += grad_peel.y;
        }
        //


        return grad;
}
//
//
//
#if 0
struct point module_potentialDerivatives_totalGradient_5_SOA_AVX(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos)
{
        asm volatile("# module_potentialDerivatives_totalGradient_SIS_SOA begins");
        printf("# module_potentialDerivatives_totalGradient_SIS_SOA_AVX begins\n");
        //
        struct point grad, result;
	//
        grad.x = 0;
        grad.y = 0;
	//
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
#endif
//
//
//
struct point module_potentialDerivatives_totalGradient_8_SOA_AVX(const struct point *pImage, const struct Potential_SOA *lens, const int shalos, const int nhalos)
{
	asm volatile("# module_potentialDerivatives_totalGradient_8_SOA_AVX begins");
	//
        struct point grad, clumpgrad;
        grad.x = 0;
        grad.y = 0;
        //
        // smearing the image coordinates on registers
        __m256d image_x  = _mm256_set1_pd(pImage->x);
        __m256d image_y  = _mm256_set1_pd(pImage->y);
	//
	__m256d __grad_x = _mm256_set1_pd(0.);
	__m256d __grad_y = _mm256_set1_pd(0.);
        //
        int i;
	int imax = shalos + nhalos;
//#pragma unroll
        for(i = shalos; i < imax - imax%4; i = i + 4)
        {
		//IACA_START;
                //
		__m256d two   = _mm256_set1_pd( 2. );
		__m256d one   = _mm256_set1_pd( 1. );
		__m256d zero  = _mm256_set1_pd( 0. );
		__m256d half  = _mm256_set1_pd( 0.5);
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
                __m256d theta     = _mm256_loadu_pd(&lens->ellipticity_angle[i]);
                /*positionning at the potential center*/
		__m256d cos_theta;
		__m256d sin_theta = _mm256_sincos_pd(&cos_theta, theta);		
		// rotation: 6 ops
		__m256d x         = _mm256_add_pd(_mm256_mul_pd(true_coord_x, cos_theta), _mm256_mul_pd(true_coord_y, sin_theta));
		__m256d y         = _mm256_sub_pd(_mm256_mul_pd(true_coord_y, cos_theta), _mm256_mul_pd(true_coord_x, sin_theta));
		//
                __m256d sqe   = __SQRT(eps);
                //
                __m256d cx1  = one_minus_eps*one_plus_eps_rcp;    // (1. - eps)/(1. + eps); 3 ops 
                __m256d cxro = one_plus_eps*one_plus_eps;         // (1. + eps)*(1. + eps); 3 ops
                __m256d cyro = one_minus_eps*one_minus_eps;       // (1. - eps)*(1. - eps); 3 ops
                __m256d rem2 = x*x*__INV(cxro) + y*y*__INV(cyro); // x*x/(cxro) + y*y/(cyro); ~5 ops
                //      
                __m256d zci_re  = zero;
                __m256d zci_im  = mhalf*(one - eps*eps)*__INV(sqe); // ~4 ops
                __m256d znum_re = cx1*x;
                __m256d znum_im = two*sqe*__SQRT(rc*rc + rem2) - y*__INV(cx1); // ~4 ops
                //
                __m256d zden_re = x;
                __m256d zden_im = _mm256_mul_pd(_mm256_set1_pd(2.), _mm256_mul_pd(rc, sqe)); 
		//
		zden_im = _mm256_sub_pd(zden_im, y);
                __m256d norm    = (zden_re*zden_re + zden_im*zden_im);     
                __m256d zis_re  = (znum_re*zden_re + znum_im*zden_im)*__INV(norm); // 3 ops
                __m256d zis_im  = (znum_im*zden_re - znum_re*zden_im)*__INV(norm); // 3 ops
		//
                norm    	= zis_re;
		//
                zis_re  = _mm256_log_pd(__SQRT(norm*norm + zis_im*zis_im));  // 3 ops// ln(zis) = ln(|zis|)+i.Arg(zis)
                zis_re  = _mm256_log_pd(__SQRT(norm*norm + zis_im*zis_im));  // 3 ops// ln(zis) = ln(|zis|)+i.Arg(zis)
                zis_im  = _mm256_atan2_pd(zis_im, norm); //
		//
                __m256d zres_re = zci_re*zis_re - zci_im*zis_im;   // Re( zci*ln(zis) ) 3 ops
                __m256d zres_im = zci_im*zis_re + zis_im*zci_re;   // Im( zci*ln(zis) ) 3 ops
		//
		zis_re = b0*zres_re;
		zis_im = b0*zres_im;
		//
                sin_theta = _mm256_sincos_pd(&cos_theta, zero - theta);
                // rotation: 6 ops
		__grad_x    = __grad_x + _mm256_add_pd(_mm256_mul_pd(zis_re, cos_theta), _mm256_mul_pd(zis_im, sin_theta));
		__grad_y    = __grad_y + _mm256_sub_pd(_mm256_mul_pd(zis_im, cos_theta), _mm256_mul_pd(zis_re, sin_theta));
		//
	}
	//
	grad.x  = ((double*) &__grad_x)[0];
	grad.x += ((double*) &__grad_x)[1];
	grad.x += ((double*) &__grad_x)[2];
	grad.x += ((double*) &__grad_x)[3];
	//
	grad.y  = ((double*) &__grad_y)[0];
	grad.y += ((double*) &__grad_y)[1];
	grad.y += ((double*) &__grad_y)[2];
	grad.y += ((double*) &__grad_y)[3];
	//
	// end of peeling
	//
	if (nhalos%4 > 0)
	{
		struct point grad_peel;	
		grad_peel = module_potentialDerivatives_totalGradient_8_SOA(pImage, lens, i, nhalos%4); 
		//
		grad.x += grad_peel.x;
		grad.y += grad_peel.y;
	}
	//
	return(grad);
}
//
//
//
struct point module_potentialDerivatives_totalGradient_8_SOA_AVX_v2(const struct point *pImage, const struct Potential_SOA *lens, const int shalos, const int nhalos)
{
        asm volatile("# module_potentialDerivatives_totalGradient_8_SOA_AVX begins");
        //
        struct point grad, clumpgrad;
        grad.x = 0;
        grad.y = 0;
        //
        // smearing the image coordinates on registers
        __m256d image_x  = _mm256_set1_pd(pImage->x);
        __m256d image_y  = _mm256_set1_pd(pImage->y);
        //
        __m256d __grad_x = _mm256_set1_pd(0.);
	__m256d __grad_y = _mm256_set1_pd(0.);
	//
	__m256d two      = _mm256_set1_pd( 2. );
	__m256d one      = _mm256_set1_pd( 1. );
	__m256d zero     = _mm256_set1_pd( 0. );
	__m256d half     = _mm256_set1_pd( 0.5);
	__m256d mhalf    = _mm256_set1_pd(-0.5);
	//
	int i;
	int imax = shalos + nhalos;
//#pragma unroll
        for(i = shalos; i < imax - imax%4; i = i + 4)
        {
		//printf("i = %d lens = %p\n", i, &lens->anglecos[i]);fflush(stdout);
                //IACA_START;
                //
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
                // 2 load
		__m256d cose 	= _mm256_loadu_pd(&lens->anglecos[i]);
		__m256d sine 	= _mm256_loadu_pd(&lens->anglesin[i]);
		//
		__m256d x	=  ADD(MUL(true_coord_x, cose), MUL(true_coord_y, sine));
		__m256d y	=  SUB(MUL(true_coord_y, cose), MUL(true_coord_x, sine));
                //
                __m256d sqe   = __SQRT(eps);
                //
                __m256d cx1  = one_minus_eps*one_plus_eps_rcp;    // (1. - eps)/(1. + eps); 3 ops
                __m256d cxro = one_plus_eps*one_plus_eps;         // (1. + eps)*(1. + eps); 3 ops
                __m256d cyro = one_minus_eps*one_minus_eps;       // (1. - eps)*(1. - eps); 3 ops
                __m256d rem2 = x*x*__INV(cxro) + y*y*__INV(cyro); // x*x/(cxro) + y*y/(cyro); ~5 ops
                //
                __m256d zci_re  = zero;
                __m256d zci_im  = mhalf*(one - eps*eps)*__INV(sqe); // ~4 ops
                __m256d znum_re = cx1*x;
                __m256d znum_im = two*sqe*__SQRT(rc*rc + rem2) - y*__INV(cx1); // ~4 ops
                //
                __m256d zden_re = x;
                __m256d zden_im = _mm256_mul_pd(_mm256_set1_pd(2.), _mm256_mul_pd(rc, sqe));
                //
                zden_im = _mm256_sub_pd(zden_im, y);
		//
                __m256d norm    = (zden_re*zden_re + zden_im*zden_im);
                __m256d zis_re  = (znum_re*zden_re + znum_im*zden_im)*__INV(norm); // 3 ops
                __m256d zis_im  = (znum_im*zden_re - znum_re*zden_im)*__INV(norm); // 3 ops
                //
                norm            = zis_re;
                //
                zis_re  = _mm256_log_pd(__SQRT(norm*norm + zis_im*zis_im));  // 3 ops// ln(zis) = ln(|zis|)+i.Arg(zis)
                zis_re  = _mm256_log_pd(__SQRT(norm*norm + zis_im*zis_im));  // 3 ops// ln(zis) = ln(|zis|)+i.Arg(zis)
                zis_im  = _mm256_atan2_pd(zis_im, norm); //
                //
                __m256d zres_re = zci_re*zis_re - zci_im*zis_im;   // Re( zci*ln(zis) ) 3 ops
                __m256d zres_im = zci_im*zis_re + zis_im*zci_re;   // Im( zci*ln(zis) ) 3 ops
                //
                zis_re = b0*zres_re;
                zis_im = b0*zres_im;
                //
                __grad_x    = __grad_x + SUB(MUL(zis_re, cose), MUL(zis_im, sine));
                __grad_y    = __grad_y + ADD(MUL(zis_im, cose), MUL(zis_re, sine));
                //
        }
        //
        grad.x  = ((double*) &__grad_x)[0];
        grad.x += ((double*) &__grad_x)[1];
        grad.x += ((double*) &__grad_x)[2];
        grad.x += ((double*) &__grad_x)[3];
        //
        grad.y  = ((double*) &__grad_y)[0];
        grad.y += ((double*) &__grad_y)[1];
        grad.y += ((double*) &__grad_y)[2];
        grad.y += ((double*) &__grad_y)[3];
        //
        // end of peeling
        //
        if (nhalos%4 > 0)
        {
                struct point grad_peel;
                grad_peel = module_potentialDerivatives_totalGradient_8_SOA(pImage, lens, i, nhalos%4);
                //
                grad.x += grad_peel.x;
                grad.y += grad_peel.y;
        }
        //
        return(grad);
}
//
//
//
struct point module_potentialDerivatives_totalGradient_81_SOA_AVX(const struct point *pImage, const struct Potential_SOA *lens, const int shalos, const int nhalos)
{
        struct point grad, clumpgrad;
        grad.x = 0;
        grad.y = 0;
        //
        // smearing the image coordinates on registers
        __m256d image_x  = _mm256_set1_pd(pImage->x);
        __m256d image_y  = _mm256_set1_pd(pImage->y);
        //
        __m256d __grad_x = _mm256_set1_pd(0.);
        __m256d __grad_y = _mm256_set1_pd(0.);
        //
        int i;
	int imax = shalos + nhalos;
//#pragma unroll
	for(i = shalos; i < imax - imax%4; i = i + 4)
        //for(i = 0; i < nhalos - nhalos%4; i = i + 4)
        {
                //IACA_START;
                //
                __m256d two   = _mm256_set1_pd( 2. );
                __m256d one   = _mm256_set1_pd( 1. );
                __m256d zero  = _mm256_set1_pd( 0. );
                __m256d half  = _mm256_set1_pd( 0.5);
                __m256d mhalf = _mm256_set1_pd(-0.5);
                //
                // 2 loads
                __m256d true_coord_x      = _mm256_sub_pd(image_x, _mm256_loadu_pd(&lens->position_x[i]));
                __m256d true_coord_y      = _mm256_sub_pd(image_y, _mm256_loadu_pd(&lens->position_y[i]));
                // 2 loads
                __m256d rc                = _mm256_loadu_pd(&lens->rcore[i]);
                __m256d rcut              = _mm256_loadu_pd(&lens->rcut[i]);
                __m256d b0                = _mm256_loadu_pd(&lens->b0[i]);
                __m256d t05 		  = b0*rcut*__INV(rcut - rc);
                //__m256d t05 		  = b0; //*rcut*__INV(rcut - rc);
                // 2 loads
                __m256d eps               = _mm256_loadu_pd(&lens->ellipticity_potential[i]);
                //
                __m256d one_minus_eps     = _mm256_sub_pd(one, eps);
                __m256d one_plus_eps      = _mm256_add_pd(one, eps);
                __m256d one_plus_eps_rcp  = __INV(one_plus_eps);
                // 1 load
                __m256d theta     = _mm256_loadu_pd(&lens->ellipticity_angle[i]);
                /*positionning at the potential center*/
                //__m256d cos_theta = _mm256_cos_pd(theta);
                //__m256d sin_theta = _mm256_sin_pd(theta);
                __m256d cos_theta;
                __m256d sin_theta = _mm256_sincos_pd(&cos_theta, theta);
                // rotation: 6 ops
                __m256d x         = _mm256_add_pd(_mm256_mul_pd(true_coord_x, cos_theta), _mm256_mul_pd(true_coord_y, sin_theta));
                __m256d y         = _mm256_sub_pd(_mm256_mul_pd(true_coord_y, cos_theta), _mm256_mul_pd(true_coord_x, sin_theta));
                //
                __m256d sqe   = __SQRT(eps);
                //
                __m256d cx1  = one_minus_eps*one_plus_eps_rcp;    // (1. - eps)/(1. + eps); 3 ops
                __m256d cxro = one_plus_eps*one_plus_eps;         // (1. + eps)*(1. + eps); 3 ops
                __m256d cyro = one_minus_eps*one_minus_eps;       // (1. - eps)*(1. - eps); 3 ops
                __m256d rem2 = x*x*__INV(cxro) + y*y*__INV(cyro); // x*x/(cxro) + y*y/(cyro); ~5 ops
		//
                __m256d zci_re  = zero;
                __m256d zci_im  = mhalf*(one - eps*eps)*__INV(sqe); // ~4 ops
		//
		__m256d znum_re, znum_im;
		__m256d zden_re = zero, zden_im = zero;
		__m256d norm;
		//
		__m256d zis_re  = zero, zis_im  = zero;
		__m256d zres_rc_re = zero, zres_rc_im = zero;
		__m256d zres_rcut_re = zero, zres_rcut_im = zero;
                //
		// step 1
		//
		{
			//  2.*sqe*sqrt(rc*rc + rem2) - y/cx1, 7 ops
			znum_re = cx1*x;
			znum_im = two*sqe*__SQRT(rc*rc + rem2) - y*__INV(cx1); // ~4 ops
			//
			zden_re = x;
			zden_im = _mm256_mul_pd(_mm256_set1_pd(2.), _mm256_mul_pd(rc, sqe));
			//
			zden_im = _mm256_sub_pd(zden_im, y);
			//norm    = (zden.re*zden.re + zden.im*zden.im); 3 ops
			norm    = (zden_re*zden_re + zden_im*zden_im);
			zis_re  = (znum_re*zden_re + znum_im*zden_im)*__INV(norm); // 3 ops
			zis_im  = (znum_im*zden_re - znum_re*zden_im)*__INV(norm); // 3 ops
			//
			norm            = zis_re;
			//
			zis_re  = _mm256_log_pd(__SQRT(norm*norm + zis_im*zis_im));  // 3 ops// ln(zis) = ln(|zis|)+i.Arg(zis)
			zis_im  = _mm256_atan2_pd(zis_im, norm); //
			//
			zres_rc_re = zci_re*zis_re - zci_im*zis_im;   // Re( zci*ln(zis) ) 3 ops
			zres_rc_im = zci_im*zis_re + zis_im*zci_re;   // Im( zci*ln(zis) ) 3 ops
		}
                //
		// step 2
		//
		{
			//  2.*sqe*sqrt(rcut*rcut + rem2) - y/cx1, 7 ops
			znum_re = cx1*x;
			znum_im = two*sqe*__SQRT(rcut*rcut + rem2) - y*__INV(cx1); // ~4 ops
			//
			zden_re = x;
			zden_im = _mm256_mul_pd(_mm256_set1_pd(2.), _mm256_mul_pd(rcut, sqe));
			zden_im = _mm256_sub_pd(zden_im, y);
			//
			norm    = (zden_re*zden_re + zden_im*zden_im); // 3 ops
			//
			zis_re  = (znum_re*zden_re + znum_im*zden_im)*__INV(norm); // 3 ops
			zis_im  = (znum_im*zden_re - znum_re*zden_im)*__INV(norm); // 3 ops
			//
			norm            = zis_re;
			//
			zis_re  = _mm256_log_pd(__SQRT(norm*norm + zis_im*zis_im));  // 3 ops// ln(zis) = ln(|zis|)+i.Arg(zis)
			zis_im  = _mm256_atan2_pd(zis_im, norm); //
			//
			zres_rcut_re = (zci_re*zis_re - zci_im*zis_im);   // Re( zci*ln(zis) ) 3 ops
			zres_rcut_im = (zci_im*zis_re + zis_im*zci_re);   // Im( zci*ln(zis) ) 3 ops		
		}
		//
		// postlogue
		//
                zis_re    = t05*(zres_rc_re - zres_rcut_re);
                zis_im    = t05*(zres_rc_im - zres_rcut_im);
                //
                //cos_theta   = _mm256_cos_pd(zero - theta);
                //sin_theta   = _mm256_sin_pd(zero - theta);
		sin_theta = _mm256_sincos_pd(&cos_theta, zero - theta);
                // rotation: 6 ops
                __grad_x   = __grad_x + _mm256_add_pd(_mm256_mul_pd(zis_re, cos_theta), _mm256_mul_pd(zis_im, sin_theta));
                __grad_y   = __grad_y + _mm256_sub_pd(_mm256_mul_pd(zis_im, cos_theta), _mm256_mul_pd(zis_re, sin_theta));
        }
        //
        //
        //
        grad.x  = ((double*) &__grad_x)[0];
        grad.x += ((double*) &__grad_x)[1];
        grad.x += ((double*) &__grad_x)[2];
        grad.x += ((double*) &__grad_x)[3];
        //
        grad.y  = ((double*) &__grad_y)[0];
        grad.y += ((double*) &__grad_y)[1];
        grad.y += ((double*) &__grad_y)[2];
        grad.y += ((double*) &__grad_y)[3];
        //
        // end of peeling
        /*
        if (nhalos%4 > 0)
        {
                struct point grad_peel;
                grad_peel = module_potentialDerivatives_totalGradient_SOA(pImage, lens, i, nhalos);
                //grad_peel = rotateCoordinateSystem(grad_peel, -theta);
                //
                grad.x += grad_peel.x;
                grad.y += grad_peel.y;
        }
        */
	for (; i < nhalos; ++i)
	{
		struct point true_coord;
		true_coord.x   = pImage->x - lens->position_x[i];
		true_coord.y   = pImage->y - lens->position_y[i];
		//printf("x, y = %f, %f\n", lens->position.x, lens->position.y);
		struct point true_coord_rot = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[i]);
		//	
		complex zis, zis_cut;
		double t05 = lens->b0[i]*lens->rcut[i]/(lens->rcut[i] - lens->rcore[i]);
		zis     = piemd_1derivatives_ci05(true_coord_rot.x, true_coord_rot.y, lens->ellipticity_potential[i], lens->rcore[i]);
		//
		zis_cut = piemd_1derivatives_ci05(true_coord_rot.x, true_coord_rot.y, lens->ellipticity_potential[i], lens->rcut[i]);

		struct point clumpgrad;
		clumpgrad.x = t05*(zis.re - zis_cut.re);
		clumpgrad.y = t05*(zis.im - zis_cut.im);
		clumpgrad = rotateCoordinateSystem(clumpgrad, -lens->ellipticity_angle[i]);
		//
		grad.x += clumpgrad.x;
		grad.y += clumpgrad.y;
	} 
        //IACA_END;
        //
        return(grad);
}


typedef struct point (*halo_func_avx_t) (const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);
halo_func_avx_t halo_func_avx[100] =
{
0, 0, 0, 0, 0, module_potentialDerivatives_totalGradient_5_SOA_AVX, 0, 0, module_potentialDerivatives_totalGradient_8_SOA_AVX, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0,   module_potentialDerivatives_totalGradient_81_SOA_AVX, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};
//
//
//
struct point module_potentialDerivatives_totalGradient_SOA_AVX(const struct point *pImage, const struct Potential_SOA *lens, int nhalos)
{
        struct point grad, clumpgrad;
        //
        grad.x = clumpgrad.x = 0;
        grad.y = clumpgrad.y = 0;
        //
        int shalos = 0;
        //
        /*
        int* p_type = &(lens->type)[0];
        int* lens_type = (int*) malloc(nhalos*sizeof(int));
        memcpy(lens_type, &(lens->type)[0], nhalos*sizeof(int));
        */
        //
        while (shalos < nhalos)
        {
                int lens_type = lens->type[shalos];
                int count     = 1;
                while (lens->type[shalos + count] == lens_type) count++;
                //std::cout << "type = " << lens_type << " " << count << " " << shalos << std::endl;
                //
                clumpgrad = (*halo_func_avx[lens_type])(pImage, lens, shalos, count);
                //
                grad.x += clumpgrad.x;
                grad.y += clumpgrad.y;
                shalos += count;
        }

        return(grad);
}

#endif

