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

@brief: Function for single order computation on GPUs

*/

#ifndef GRADIENT_GPU_CUH_
#define GRADIENT_GPU_CUH_

//#include "cudafunctions.cuh"
#include <fstream>
#include <structure_hpc.hpp>
//
__global__
void
module_potentialDerivatives_totalGradient_SOA_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nhalos, int nbgridcells);
//
__global__
void
module_potentialDerivatives_totalGradient_SOA_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart = 0, int jstart = 0);

#if 0
__device__ struct point module_potentialDerivatives_totalGradient_SOA_GPU(const struct point *pImage, const struct Potential_SOA *lens, int nhalos);

__device__ inline struct point rotateCoordinateSystem_GPU(struct point P, double theta);
__device__ inline struct point rotateCoordinateSystem_GPU_2(struct point P, double cosi, double sinu);
//
//__device__
__device__ point module_potentialDerivatives_totalGradient_5_SOA_GPU(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);
__device__ point module_potentialDerivatives_totalGradient_8_SOA_GPU(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);
__device__ point module_potentialDerivatives_totalGradient_81_SOA_GPU(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);

__global__ void module_potentialDerivatives_totalGradient_SOA_GPU(double *grid_grad_x, double *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int nhalos);

__global__
void
module_potentialDerivatives_totalGradient_8_SOA_GPU_cur(double *grid_grad_x, double *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int shalos, int nhalos);

__global__
void
module_potentialDerivatives_totalGradient_8_SOA_GPU_SM2(double *grid_grad_x, double *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int shalos, int nhalos);

__global__
void
module_potentialDerivatives_totalGradient_8_SOA_GPU_SM3(double *grid_grad_x, double *grid_grad_y, const struct Potential_SOA
 *lens, const struct grid_param *frame, int nbgridcells, int shalos, int nhalos);


__global__
void
module_potentialDerivatives_totalGradient_8_SOA_GPU_SM4(double *grid_grad_x, double *grid_grad_y, const struct Potential_SOA lens, const struct grid_param *frame, int nbgridcells, int shalos, int nhalos/*, double* dtimer*/);

__global__
void
module_potentialDerivatives_totalGradient_8_SOA_GPU_v2(double *grid_grad_x, double *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int i, int nhalos);




__global__
void
gradient_grid_kernel(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcells, const struct Potential_SOA *lens);

__global__
void
gradient_grid_kernel_v2(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcells, const struct Potential_SOA *lens);

__global__
void
gradient_grid_piemd_GPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcell, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut);

__global__
void
gradient_grid_sis_GPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcell, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut);

__global__
void
gradient_grid_piemd_GPU_multiple(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcell, int Ndevice, int indexactual, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut);

__global__
void
gradient_grid_kernel_multiple(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, int nbgridcells, const struct Potential_SOA *lens, int indexactual, int ncells);

__device__
static
double
atomicAdd_double(double* address, double val);
#endif

#define KERNEL_8 \
double x = true_coord_x*cosi[i] + true_coord_y*sinu[i]; \
double y = true_coord_y*cosi[i] - true_coord_x*sinu[i]; \
double rem2 = x*x*inv_onepeps + y*y*inv_onemeps; \
double norm; \
complex      zis; \
complex      znum; \
complex      zden; \
complex      zres; \
znum.re = cx1*x; \
znum.im = 2.*sqe*sqrt(rc*rc + rem2) - y/cx1; \
zden.re = x; \
zden.im = 2.*rc*sqe - y; \
norm    = (x*x + zden.im*zden.im); \
zis.re  = (znum.re*x + znum.im*zden.im)/norm; \
zis.im  = (znum.im*x - znum.re*zden.im)/norm; \
norm    = zis.re; \
zis.re  = log(sqrt(norm*norm + zis.im*zis.im)); \
zis.im  = atan2(zis.im, norm); \
zres.re = - zci_im*zis.im; \
zres.im =   zci_im*zis.re; \
grad_x += b0[i]*(zres.re*cosi[i] - zres.im*sinu[i]); \
grad_y += b0[i]*(zres.im*cosi[i] + zres.re*sinu[i]); \
\


#define KERNEL_8_reg(X) \
double x = true_coord_x*cosi[i] + true_coord_y*sinu[i]; \
double y = true_coord_y*cosi[i] - true_coord_x*sinu[i]; \
double rem2 = x*x*inv_onepeps + y*y*inv_onemeps; \
double norm; \
double    zis_re, zis_im; \
double    znum_re, znum_im; \
double  /*zden_re,*/ zden_im; \
double    zres_re, zres_im; \
znum_re = cx1*x; \
znum_im = 2.*sqe*sqrt(rc*rc + rem2) - y/cx1; \
/*zden_re = x;*/ \
zden_im = 2.*rc*sqe - y; \
norm    = (x*x + zden_im*zden_im); \
zis_re  = (znum_re*x + znum_im*zden_im)/norm; \
zis_im  = (znum_im*x - znum_re*zden_im)/norm; \
norm    = zis_re; \
zis_re  = log(sqrt(norm*norm + zis_im*zis_im)); \
zis_im  = atan2(zis_im, norm); \
zres_re = - zci_im*zis_im; \
zres_im =   zci_im*zis_re; \
grad_x = b0[i]*(zres_re*cosi[i] - zres_im*sinu[i]); \
grad_y = b0[i]*(zres_im*cosi[i] + zres_re*sinu[i]); \
\


#endif /* GRADIENT_GPU_CUH_ */
