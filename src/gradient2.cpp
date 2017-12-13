#include <iostream>
#include <iomanip>
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
#include "structure_hpc.hpp"
#include "gradient2.hpp"

//lenstool fct
/*
* derivates of I0.5 KK
* Parameters :
* - (x,y) is the computation position of the potential
* - eps is the ellepticity (a-b)/(a+b)
* - rc is the core radius
* - b0 asymptotic Einstein radius E0. (6pi*vdisp^2/c^2)
*
* Return a the 4 second derivatives of the PIEMD potential
*/
void mdci05_hpc(type_t x, type_t y, type_t eps, type_t rc, type_t b0, struct matrix *res)
{
    type_t   ci, sqe, cx1, cxro, cyro, wrem;
    type_t  didyre, didyim, didxre;// didxim;
    type_t  cx1inv, den1, num2, den2;

    sqe = sqrt(eps);
    cx1 = (1. - eps) / (1. + eps);
    cx1inv = 1. / cx1;
    cxro = (1. + eps) * (1. + eps);     /* rem^2=x^2/(1+e^2) + y^2/(1-e^2) Eq 2.3.6*/
    cyro = (1. - eps) * (1. - eps);
    ci = 0.5 * (1. - eps * eps) / sqe;
    wrem = sqrt(rc * rc + x * x / cxro + y * y / cyro); /*wrem^2=w^2+rem^2 with w core radius*/
    den1 = 2.*sqe * wrem - y * cx1inv;
    den1 = cx1 * cx1 * x * x + den1 * den1;
    num2 = 2.*rc * sqe - y;
    den2 = x * x + num2 * num2;

    didxre = ci * ( cx1 * (2.*sqe * x * x / cxro / wrem - 2.*sqe * wrem + y * cx1inv) / den1 + num2 / den2 );
    didyre = ci * ( (2 * sqe * x * y * cx1 / cyro / wrem - x) / den1 + x / den2 );

    didyim = ci * ( (2 * sqe * wrem * cx1inv - y * cx1inv * cx1inv - 4 * eps * y / cyro +
                     2 * sqe * y * y / cyro / wrem * cx1inv) / den1 - num2 / den2 );

    res->a = b0 * didxre;
    res->b = res->d = b0 * didyre; //(didyre+didxim)/2.;
    res->c = b0 * didyim;

//  return(res);


}

void printmat(matrix A){
	std::cerr << A.a << " " << A.b << " " << A.c << " " << A.d <<  std::endl;
}

//
// SOA versions, vectorizable
//
//
struct matrix module_potentialDerivatives_totalGradient2_81_SOA_v2(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos)
{
        asm volatile("# module_potentialDerivatives_totalGradient2_81_SOA begins");
        //std::cout << "# module_potentialDerivatives_totalGradient_81_SOA begins" << std::endl;
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        //
        type_t t05;
        struct matrix grad2, clump, clumpcore, clumpcut;
        grad2.a = clump.a = 0;
        grad2.b = clump.b = 0;
        grad2.c = clump.c = 0;
        grad2.d = clump.d = 0;
        for(int i = shalos; i < shalos + nhalos; i++)
        {
			struct point true_coord, true_coord_rotation;
			//True coord
			true_coord.x = pImage->x - lens->position_x[i];
			true_coord.y = pImage->y - lens->position_y[i];
			//
			//std::cerr << "start" << std::endl;
			//Rotation
			type_t cose = lens->anglecos[i];
			type_t sine = lens->anglesin[i];
			//
			type_t x = true_coord.x*cose + true_coord.y*sine;
			type_t y = true_coord.y*cose - true_coord.x*sine;
			// 81 comput
            t05 = lens->rcut[i] / (lens->rcut[i] - lens->rcore[i]);
            mdci05_hpc(x, y, lens->ellipticity_potential[i], lens->rcore[i], lens->b0[i], &clumpcore);
            mdci05_hpc(x, y, lens->ellipticity_potential[i], lens->rcut[i], lens->b0[i], &clumpcut);
            //
            clumpcore.a = t05 * (clumpcore.a - clumpcut.a);
            clumpcore.b = t05 * (clumpcore.b - clumpcut.b);
            clumpcore.c = t05 * (clumpcore.c - clumpcut.c);
            clumpcore.d = t05 * (clumpcore.d - clumpcut.d);
            //printmat(clumpcore);
            //rotation matrix  1
            clumpcut.a = clumpcore.a * cose + clumpcore.b * sine;
            clumpcut.b = clumpcore.a * -sine + clumpcore.b * cose;
            clumpcut.c = clumpcore.d * -sine + clumpcore.c * cose;
            clumpcut.d = clumpcore.d * cose + clumpcore.c * sine;
            //printmat(clumpcut);
            //rotation matrix  2
            clump.a = cose * clumpcut.a + sine * clumpcut.d;
            clump.b = cose * clumpcut.b + sine * clumpcut.c;
            clump.c = -sine * clumpcut.b + cose * clumpcut.c;
            clump.d = -sine * clumpcut.a + cose * clumpcut.d;
            //printmat(clump);

            grad2.a += clump.a;
            grad2.b += clump.b;
            grad2.c += clump.c;
            grad2.d += clump.d;
        }
        //
        return(grad2);
}
//
//This natrix handles the calling of the gradient functions for the different type without losing time to switch or if condition
//
typedef struct matrix (*halo_g2_func_t) (const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);
halo_g2_func_t halo_g2_func[100] =
{ 
0, 0, 0, 0, 0, 0, 0, 0, 0,  0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0,  module_potentialDerivatives_totalGradient2_81_SOA_v2, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};
//
//
//
struct matrix module_potentialDerivatives_totalGradient2_SOA(const struct point *pImage, const struct Potential_SOA *lens, int nhalos)
{
        struct matrix grad2, clumpgrad;
        //
        grad2.a = clumpgrad.a = 0;
        grad2.b = clumpgrad.b = 0;
        grad2.c = clumpgrad.c = 0;
        grad2.d = clumpgrad.d = 0;
	//
	int shalos = 0;

	while (shalos < nhalos)
	{
		int lens_type = lens->type[shalos];
		int count     = 1;
		while (lens->type[shalos + count] == lens_type and shalos + count < nhalos) count++;
		//std::cerr << "type = " << lens_type << " " << count << " " << shalos << " " << nhalos << " " <<  " func = " << halo_g2_func[lens_type] << std::endl;
		//	
		clumpgrad = (*halo_g2_func[lens_type])(pImage, lens, shalos, count);
		//
		grad2.a += clumpgrad.a;
		grad2.b += clumpgrad.b;
		grad2.c += clumpgrad.c;
		grad2.d += clumpgrad.d;
		shalos += count;
	}	

        return(grad2);
}
//



