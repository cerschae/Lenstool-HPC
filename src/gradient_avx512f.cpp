#include <iostream>
#include <string.h>
#include <cuda_runtime.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <immintrin.h>
//
#include "structure_hpc.h"
#include "simd_math_avx512f.h"
#include "gradient.hpp"
//#include "iacaMarks.h"
//
//
//
struct point module_potentialDerivatives_totalGradient_SOA_AVX512(const struct point *pImage, const struct Potential_SOA *lens, const int nhalos)
{
        struct point grad, clumpgrad;
                grad.x = 0;
                grad.y = 0;
                //
                // smearing the image coordinates on registers
                __m512d image_x  = _mm512_set1_pd(pImage->x);
                __m512d image_y  = _mm512_set1_pd(pImage->y);
                //
                __m512d __grad_x = _mm512_set1_pd(0.);
                __m512d __grad_y = _mm512_set1_pd(0.);
                //
                int i;
#pragma unroll
                for(i = 0; i < nhalos - nhalos%8; i = i + 8)
                {
                        //IACA_START;
                        //
                        __m512d two   = _mm512_set1_pd( 2. );
                        __m512d one   = _mm512_set1_pd( 1. );
                        __m512d zero  = _mm512_set1_pd( 0. );
                        __m512d half  = _mm512_set1_pd( 0.5);
                        __m512d mhalf = _mm512_set1_pd(-0.5);
                        //
                        // 2 loads
                        __m512d true_coord_x      = _mm512_sub_pd(image_x, _mm512_loadu_pd(&lens->position_x[i]));
                        __m512d true_coord_y      = _mm512_sub_pd(image_y, _mm512_loadu_pd(&lens->position_y[i]));
                        // 2 loafs
                        __m512d rc                = _mm512_loadu_pd(&lens->rcore[i]);
                        __m512d b0                = _mm512_loadu_pd(&lens->b0[i]);
                        // 2 loads
                        __m512d eps               = _mm512_loadu_pd(&lens->ellipticity_potential[i]);
                        //
                        __m512d one_minus_eps     = _mm512_sub_pd(one, eps);
                        __m512d one_plus_eps      = _mm512_add_pd(one, eps);
                        __m512d one_plus_eps_rcp  = __INV(one_plus_eps);
                        // 1 load
                        __m512d theta     = _mm512_loadu_pd(&lens->ellipticity_angle[i]);
                        /*positionning at the potential center*/
                        __m512d cos_theta = _mm512_cos_pd(theta);
                        __m512d sin_theta = _mm512_sin_pd(theta);
                        // rotation: 6 ops
                        __m512d x         = _mm512_add_pd(_mm512_mul_pd(true_coord_x, cos_theta), _mm512_mul_pd(true_coord_y, sin_theta));
                        __m512d y         = _mm512_sub_pd(_mm512_mul_pd(true_coord_y, cos_theta), _mm512_mul_pd(true_coord_x, sin_theta));
                        //
                        __m512d sqe   = __SQRT(eps);
                        //
                        __m512d cx1  = _mm512_mul_pd(one_minus_eps, one_plus_eps_rcp);    // (1. - eps)/(1. + eps); 3 ops
                        //__m512d cx1  = one_minus_eps*one_plus_eps_rcp;    // (1. - eps)/(1. + eps); 3 ops
                        __m512d cxro = one_plus_eps*one_plus_eps;         // (1. + eps)*(1. + eps); 3 ops
                        __m512d cyro = one_minus_eps*one_minus_eps;       // (1. - eps)*(1. - eps); 3 ops
                        __m512d rem2 = x*x*__INV(cxro) + y*y*__INV(cyro); // x*x/(cxro) + y*y/(cyro); ~5 ops
                        //
                        __m512d zci_re  = zero;
                        __m512d zci_im  = mhalf*(one - eps*eps)*__INV(sqe); // ~4 ops
                        //  2.*sqe*sqrt(rc*rc + rem2) - y/cx1, 7 ops
                        __m512d znum_re = cx1*x;
                        __m512d znum_im = two*sqe*__SQRT(rc*rc + rem2) - y*__INV(cx1); // ~4 ops
                        //
                        __m512d zden_re = x;
                        __m512d zden_im = _mm512_mul_pd(_mm512_set1_pd(2.), _mm512_mul_pd(rc, sqe));
                        //
                       //
                        zden_im = _mm512_sub_pd(zden_im, y);
                        //norm    = (zden.re*zden.re + zden.im*zden.im); 3 ops
                        __m512d norm    = (zden_re*zden_re + zden_im*zden_im);
                        __m512d zis_re  = (znum_re*zden_re + znum_im*zden_im)*__INV(norm); // 3 ops
                        __m512d zis_im  = (znum_im*zden_re - znum_re*zden_im)*__INV(norm); // 3 ops
                        //
                        norm            = zis_re;
                        //
                        zis_re  = _mm512_log_pd(__SQRT(norm*norm + zis_im*zis_im));  // 3 ops// ln(zis) = ln(|zis|)+i.Arg(zis)
                        zis_re  = _mm512_log_pd(__SQRT(norm*norm + zis_im*zis_im));  // 3 ops// ln(zis) = ln(|zis|)+i.Arg(zis)
                        zis_im  = _mm512_atan2_pd(zis_im, norm); //
                        //
                        __m512d zres_re = zci_re*zis_re - zci_im*zis_im;   // Re( zci*ln(zis) ) 3 ops
                        __m512d zres_im = zci_im*zis_re + zis_im*zci_re;   // Im( zci*ln(zis) ) 3 ops
                        //
                        zis_re = b0*zres_re;
                        zis_im = b0*zres_im;
                        //
                        cos_theta = _mm512_cos_pd(zero - theta);
                        sin_theta = _mm512_sin_pd(zero - theta);
                        // rotation: 6 ops
                        __grad_x    = __grad_x + _mm512_add_pd(_mm512_mul_pd(zis_re, cos_theta), _mm512_mul_pd(zis_im, sin_theta));
                        __grad_y    = __grad_y + _mm512_sub_pd(_mm512_mul_pd(zis_im, cos_theta), _mm512_mul_pd(zis_re, sin_theta));
                        //
                }
                //
                grad.x  = ((double*) &__grad_x)[0];
                grad.x += ((double*) &__grad_x)[1];
                grad.x += ((double*) &__grad_x)[2];
                grad.x += ((double*) &__grad_x)[3];
                grad.x += ((double*) &__grad_x)[4];
                grad.x += ((double*) &__grad_x)[5];
                grad.x += ((double*) &__grad_x)[6];
                grad.x += ((double*) &__grad_x)[7];
                //
                grad.y  = ((double*) &__grad_y)[0];
                grad.y += ((double*) &__grad_y)[1];
                grad.y += ((double*) &__grad_y)[2];
                grad.y += ((double*) &__grad_y)[3];
                grad.y += ((double*) &__grad_y)[4];
                grad.y += ((double*) &__grad_y)[5];
                grad.y += ((double*) &__grad_y)[6];
                grad.y += ((double*) &__grad_y)[7];
                //
                // end of peeling
                //
                if (nhalos%8 > 0)
                {
                        struct point grad_peel;
                        grad_peel = module_potentialDerivatives_totalGradient_SOA(pImage, lens, i, nhalos);
                        //
                        grad.x += grad_peel.x;
                        grad.y += grad_peel.y;
                }
                //
                //IACA_END;
                //
                return(grad);
}


struct point module_potentialDerivatives_totalGradient_81_SOA_AVX512(const struct point *pImage, const struct Potential_SOA *lens, const int nhalos)
{
	//_MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
        struct point grad, clumpgrad;
	grad.x = 0;
	grad.y = 0;
	//
	// smearing the image coordinates on registers
	__m512d image_x  = _mm512_set1_pd(pImage->x);
	__m512d image_y  = _mm512_set1_pd(pImage->y);
	//
	__m512d __grad_x = _mm512_set1_pd(0.);
	__m512d __grad_y = _mm512_set1_pd(0.);
	//
	int i;
#pragma unroll
	for(i = 0; i < nhalos - nhalos%8; i = i + 8)
	{
		//IACA_START;
		//
		__m512d two   = _mm512_set1_pd( 2. );
		__m512d one   = _mm512_set1_pd( 1. );
		__m512d zero  = _mm512_set1_pd( 0. );
		__m512d half  = _mm512_set1_pd( 0.5);
		__m512d mhalf = _mm512_set1_pd(-0.5);
		//
		// 2 loads
		__m512d true_coord_x      = _mm512_sub_pd(image_x, _mm512_loadu_pd(&lens->position_x[i]));
		__m512d true_coord_y      = _mm512_sub_pd(image_y, _mm512_loadu_pd(&lens->position_y[i]));
		// 3 loads
		__m512d rc                = _mm512_loadu_pd(&lens->rcore[i]);
		__m512d rcut              = _mm512_loadu_pd(&lens->rcut[i]);
		__m512d b0                = _mm512_loadu_pd(&lens->b0[i]);
		__m512d t05               = b0*rcut*__INV(rcut - rc); 
		// 1 loads
		__m512d eps               = _mm512_loadu_pd(&lens->ellipticity_potential[i]);
		//
		__m512d one_minus_eps     = _mm512_sub_pd(one, eps);
		__m512d one_plus_eps      = _mm512_add_pd(one, eps);
		__m512d one_plus_eps_rcp  = __INV(one_plus_eps);
		// 1 load
		__m512d theta     = _mm512_loadu_pd(&lens->ellipticity_angle[i]);
		/*positionning at the potential center*/
		//__m512d cos_theta = _mm512_cos_pd(theta);
		//__m512d sin_theta = _mm512_sin_pd(theta);
		__m512d cos_theta;
                __m512d sin_theta = _mm512_sincos_pd(&cos_theta, theta);
		// rotation: 6 ops
		__m512d x         = _mm512_add_pd(_mm512_mul_pd(true_coord_x, cos_theta), _mm512_mul_pd(true_coord_y, sin_theta));
		__m512d y         = _mm512_sub_pd(_mm512_mul_pd(true_coord_y, cos_theta), _mm512_mul_pd(true_coord_x, sin_theta));
		//
		__m512d sqe   = __SQRT(eps);
		//
		__m512d cx1  = _mm512_mul_pd(one_minus_eps,one_plus_eps_rcp);    // (1. - eps)/(1. + eps); 3 ops
		//__m512d cx1  = one_minus_eps*one_plus_eps_rcp;    // (1. - eps)/(1. + eps); 3 ops
		__m512d cxro = one_plus_eps*one_plus_eps;         // (1. + eps)*(1. + eps); 3 ops
		__m512d cyro = one_minus_eps*one_minus_eps;       // (1. - eps)*(1. - eps); 3 ops
		__m512d rem2 = x*x*__INV(cxro) + y*y*__INV(cyro); // x*x/(cxro) + y*y/(cyro); ~5 ops
		//
		__m512d zci_re  = zero;
		__m512d zci_im  = mhalf*(one - eps*eps)*__INV(sqe); // ~4 ops
		//
		//  2.*sqe*sqrt(rc*rc + rem2) - y/cx1, 7 ops
		//
		__m512d znum_re = zero, znum_im = zero;
		__m512d zden_re = zero, zden_im = zero;
		__m512d norm;
		//
		__m512d zis_re       = zero, zis_im       = zero;
		__m512d zres_rc_re   = zero, zres_rc_im   = zero;
		__m512d zres_rcut_re = zero, zres_rcut_im = zero;
		//
		// part 1
		//
		{
			znum_re = cx1*x;
			znum_im = two*sqe*__SQRT(rc*rc + rem2) - y*__INV(cx1); // ~4 ops
			//
			zden_re = x;
			zden_im = _mm512_mul_pd(_mm512_set1_pd(2.), _mm512_mul_pd(rc, sqe));
			zden_im = _mm512_sub_pd(zden_im, y);
			//norm    = (zden.re*zden.re + zden.im*zden.im); 3 ops
			norm    = (zden_re*zden_re + zden_im*zden_im);
			zis_re  = (znum_re*zden_re + znum_im*zden_im)*__INV(norm); // 3 ops
			zis_im  = (znum_im*zden_re - znum_re*zden_im)*__INV(norm); // 3 ops
			//
			norm            = zis_re;
			//
			zis_re  = _mm512_log_pd(__SQRT(norm*norm + zis_im*zis_im));  // 3 ops// ln(zis) = ln(|zis|)+i.Arg(zis)
			zis_im  = _mm512_atan2_pd(zis_im, norm); //
			//
			zres_rc_re = zci_re*zis_re - zci_im*zis_im;   // Re( zci*ln(zis) ) 3 ops
			zres_rc_im = zci_im*zis_re + zis_im*zci_re;   // Im( zci*ln(zis) ) 3 ops
		}
		//
		// part 2
		//
		{
			znum_re = cx1*x;
			znum_im = two*sqe*__SQRT(rcut*rcut + rem2) - y*__INV(cx1); // ~4 ops
			//
			zden_re = x;
			zden_im = _mm512_mul_pd(_mm512_set1_pd(2.), _mm512_mul_pd(rcut, sqe));
			zden_im = _mm512_sub_pd(zden_im, y);
			//norm    = (zden.re*zden.re + zden.im*zden.im); 3 ops
			//
			norm    = (zden_re*zden_re + zden_im*zden_im);
			//
			zis_re  = (znum_re*zden_re + znum_im*zden_im)*__INV(norm); // 3 ops
			zis_im  = (znum_im*zden_re - znum_re*zden_im)*__INV(norm); // 3 ops
			//
			norm            = zis_re;
			//
			zis_re  = _mm512_log_pd(__SQRT(norm*norm + zis_im*zis_im));  // 3 ops// ln(zis) = ln(|zis|)+i.Arg(zis)
			zis_im  = _mm512_atan2_pd(zis_im, norm); //
			//
			zres_rcut_re = (zci_re*zis_re - zci_im*zis_im);   // Re( zci*ln(zis) ) 3 ops
			zres_rcut_im = (zci_im*zis_re + zis_im*zci_re);   // Im( zci*ln(zis) ) 3 ops
		}
		//
		//
		//
		zis_re = t05*(zres_rc_re - zres_rcut_re);
		zis_im = t05*(zres_rc_im - zres_rcut_im);
		//
		//cos_theta = _mm512_cos_pd(zero - theta);
		//sin_theta = _mm512_sin_pd(zero - theta);
		cos_theta;
                sin_theta = _mm512_sincos_pd(&cos_theta, zero - theta);
		// rotation: 6 ops
		__grad_x    = __grad_x + _mm512_add_pd(_mm512_mul_pd(zis_re, cos_theta), _mm512_mul_pd(zis_im, sin_theta));
		__grad_y    = __grad_y + _mm512_sub_pd(_mm512_mul_pd(zis_im, cos_theta), _mm512_mul_pd(zis_re, sin_theta));
		//
	}
	//
	grad.x  = ((double*) &__grad_x)[0];
	grad.x += ((double*) &__grad_x)[1];
	grad.x += ((double*) &__grad_x)[2];
	grad.x += ((double*) &__grad_x)[3];
	grad.x += ((double*) &__grad_x)[4];
	grad.x += ((double*) &__grad_x)[5];
	grad.x += ((double*) &__grad_x)[6];
	grad.x += ((double*) &__grad_x)[7];
	//
	grad.y  = ((double*) &__grad_y)[0];
	grad.y += ((double*) &__grad_y)[1];
	grad.y += ((double*) &__grad_y)[2];
	grad.y += ((double*) &__grad_y)[3];
	grad.y += ((double*) &__grad_y)[4];
	grad.y += ((double*) &__grad_y)[5];
	grad.y += ((double*) &__grad_y)[6];
	grad.y += ((double*) &__grad_y)[7];
	//
	// end of peeling
	// 
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
	//
	//IACA_END;
	//
	return(grad);
}
