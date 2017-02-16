
#include "gradient_GPU.cuh"

#include <iostream>
#include <string.h>
#include <cuda_runtime.h>
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
#include "utils.hpp"



//
// SOA versions, vectorizable
//
__device__ struct point module_potentialDerivatives_totalGradient_5_SOA_GPU(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos)
{
        //asm volatile("# module_potentialDerivatives_totalGradient_SIS_SOA begins");
  //
  struct point grad, clumpgrad;
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
    true_coord_rotation = rotateCoordinateSystem_GPU(true_coord, lens->ellipticity_angle[i]);
    double R = sqrt(true_coord_rotation.x*true_coord_rotation.x*(1 - lens->ellipticity_potential[i])+true_coord_rotation.y*true_coord_rotation.y*(1 + lens->ellipticity_potential[i]));
    //
    grad.x += (1 - lens->ellipticity[i]/3.)*lens->b0[i]*true_coord_rotation.x/R;
    grad.y += (1 + lens->ellipticity[i]/3.)*lens->b0[i]*true_coord_rotation.y/R;
  }
  return grad;
}

//
//
//
__device__ struct point module_potentialDerivatives_totalGradient_8_SOA_GPU(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos)
{
  //asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
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
    true_coord.x = pImage->x - lens->position_x[i];
    true_coord.y = pImage->y - lens->position_y[i];
    double cosi,sinu;
    //cosi = cos(lens->ellipticity_angle[i]);
    //sinu = sin(lens->ellipticity_angle[i]);
    cosi = lens->anglecos[i];
    sinu = lens->anglesin[i];
    //positionning at the potential center
    // Change the origin of the coordinate system to the center of the clump
    //true_coord_rot = rotateCoordinateSystem_GPU(true_coord, lens->ellipticity_angle[i]);
    true_coord_rot = rotateCoordinateSystem_GPU_2(true_coord, cosi,sinu);
    //
    double x   = true_coord_rot.x;
    double y   = true_coord_rot.y;
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
    // rotation
    clumpgrad.x = zis.re;
    clumpgrad.y = zis.im;
    //clumpgrad = rotateCoordinateSystem_GPU(true_coord, -lens->ellipticity_angle[i]);
    clumpgrad = rotateCoordinateSystem_GPU_2(clumpgrad, cosi,-sinu);
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
    //}
  }
  //IACA_END;
  //
  return(grad);
}


__device__ struct point module_potentialDerivatives_totalGradient_81_SOA_GPU(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos)
{
        //asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
  //std::cout << "# module_potentialDerivatives_totalGradient_SOA begins" << std::endl;
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
                true_coord.x = pImage->x - lens->position_x[i];
                true_coord.y = pImage->y - lens->position_y[i];
                //positionning at the potential center
                // Change the origin of the coordinate system to the center of the clump
                true_coord_rot = rotateCoordinateSystem_GPU(true_coord, lens->ellipticity_angle[i]);
                //
                double x    = true_coord_rot.x;
                double y    = true_coord_rot.y;
                double eps  = lens->ellipticity_potential[i];
                double rc   = lens->rcore[i];
                double rcut = lens->rcut[i];
    double b0   = lens->b0[i];
    double t05  = b0*rcut/(rcut - rc);
    //printf("b0 = %f, rcut = %f, rc = %f, t05 = %f\n", b0, rcut, rc, t05);
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
                complex zci, znum, zden, zres_rc, zres_rcut;
                double norm;
                //
                zci.re  = 0;
                zci.im  = -0.5*(1. - eps*eps)/sqe;
                //
    // step 1
    //
    {
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
      zres_rc.re = zci.re*zis.re - zci.im*zis.im;   // Re( zci*ln(zis) )
      zres_rc.im = zci.im*zis.re + zis.im*zci.re;   // Im( zci*ln(zis) )
    }
    //
    // step 2
    //
    {
                        znum.re = cx1*x;
                        znum.im = 2.*sqe*sqrt(rcut*rcut + rem2) - y/cx1;
                        //
                        zden.re = x;
                        zden.im = 2.*rcut*sqe - y;
                        norm    = (zden.re*zden.re + zden.im*zden.im);     // zis = znum/zden
                        //
                        zis.re  = (znum.re*zden.re + znum.im*zden.im)/norm;
                        zis.im  = (znum.im*zden.re - znum.re*zden.im)/norm;
                        norm    = zis.re;
                        zis.re  = log(sqrt(norm*norm + zis.im*zis.im));  // ln(zis) = ln(|zis|)+i.Arg(zis)
                        zis.im  = atan2(zis.im, norm);
                        //  norm = zis.re;
                        zres_rcut.re = zci.re*zis.re - zci.im*zis.im;   // Re( zci*ln(zis) )
                        zres_rcut.im = zci.im*zis.re + zis.im*zci.re;   // Im( zci*ln(zis) )
                }
    zis.re  = t05*(zres_rc.re - zres_rcut.re);
    zis.im  = t05*(zres_rc.im - zres_rcut.im);
    //printf("%f %f\n", zis.re, zis.im);
                //
                //zres.re = zis.re*b0;
                //zres.im = zis.im*b0;
                // rotation
                clumpgrad.x = zis.re;
                clumpgrad.y = zis.im;
                clumpgrad = rotateCoordinateSystem_GPU(clumpgrad, -lens->ellipticity_angle[i]);
                //
                //clumpgrad.x = lens->b0[i]*clumpgrad.x;
                //clumpgrad.y = lens->b0[i]*clumpgrad.y;
                //
                //clumpgrad.x = lens->b0[i]*zis.re;
                //clumpgrad.y = lens->b0[i]*zis.im;
                //nan check
                //if(clumpgrad.x == clumpgrad.x or clumpgrad.y == clumpgrad.y)
                //{
                        // add the gradients
                grad.x += clumpgrad.x;
                grad.y += clumpgrad.y;
    //printf("grad = %f %f\n", clumpgrad.x, clumpgrad.y);
                //}
        }
        //IACA_END;
        //
        return(grad);
}
//
//
//

typedef struct point (*halo_func_GPU_t) (const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);

__constant__ halo_func_GPU_t halo_func_GPU[100] =
{
0, 0, 0, 0, 0, module_potentialDerivatives_totalGradient_5_SOA_GPU, 0, 0, module_potentialDerivatives_totalGradient_8_SOA_GPU,  0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0,  module_potentialDerivatives_totalGradient_81_SOA_GPU, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};
//
//
//
__device__ struct point module_potentialDerivatives_totalGradient_SOA_GPU(const struct point *pImage, const struct Potential_SOA *lens, int nhalos)
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

  //quicksort(lens_type, nhalos);
  //*/

  //printf ("%f \n",lens->position_x[shalos]);


  while (shalos < nhalos)
  {

    int lens_type = lens->type[shalos];
    int count     = 1;
    while (lens->type[shalos + count] == lens_type) count++;
    //std::cerr << "type = " << lens_type << " " << count << " " << shalos << std::endl;
    //printf ("%d %d %d \n",lens_type,count,shalos);
    //
    clumpgrad = (*halo_func_GPU[lens_type])(pImage, lens, shalos, count);
    //
    grad.x += clumpgrad.x;
    grad.y += clumpgrad.y;
    shalos += count;
  }

        return(grad);
}

__device__ inline struct point rotateCoordinateSystem_GPU(struct point P, double theta)
{
  struct  point   Q;

  Q.x = P.x*cos(theta) + P.y*sin(theta);
  Q.y = P.y*cos(theta) - P.x*sin(theta);

  return(Q);
}

__device__ inline struct point rotateCoordinateSystem_GPU_2(struct point P, double cosi, double sinu)
{
  struct  point   Q;

  Q.x = P.x*cosi + P.y*sinu;
  Q.y = P.y*cosi - P.x*sinu;

  return(Q);
}
