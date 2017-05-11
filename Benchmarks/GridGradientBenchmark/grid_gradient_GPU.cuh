/*
 * gradientgpu.cuh
 *
 *  Created on: Nov 29, 2016
 *      Author: cerschae
 */

#ifndef GRID_GRADIENT_GPU_CUH_
#define GRID_GRADIENT_GPU_CUH_

#include "cudafunctions.cuh"
#include "gradient_GPU.cuh"
#include <structure_hpc.h>

void gradient_grid_GPU_sorted(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens,int nbgridcells);
void gradient_grid_pinned(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int *Nlens,int nbgridcell);
void gradient_grid_pinned_multiple(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int *Nlens,int nbgridcell);
void gradient_grid_GPU_sub(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos,int nbgridcells, int indexactual, int Ncells );
void gradient_grid_GPU_multiple(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells);

__global__ void gradient_grid_kernel(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcells, const struct Potential_SOA *lens);
__global__ void gradient_grid_kernel_v2(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcells, const struct Potential_SOA *lens);
__global__ void gradient_grid_piemd_GPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcell, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut);
__global__ void gradient_grid_sis_GPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcell, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut);
__global__ void gradient_grid_piemd_GPU_multiple(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcell, int Ndevice, int indexactual, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut);
__global__ void gradient_grid_kernel_multiple(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, int nbgridcells, const struct Potential_SOA *lens, int indexactual, int ncells);

__device__ static double atomicAdd_double(double* address, double val);

#endif /* GRADIENTGPU_CUH_ */


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
grad_y += b0[i]*(zres.im*cosi[i] + zres.re*sinu[i]); 


#define KERNEL_8_reg(X) \
double x = (true_coord_x + X*dx)*cosi[i] +          true_coord_y*sinu[i]; \
double y =          true_coord_y*cosi[i] - (true_coord_x + X*dx)*sinu[i]; \
double rem2 = x*x*inv_onepeps + y*y*inv_onemeps; \
double norm; \
double      zis_re, zis_im; \
double      znum_re, znum_im; \
double      /*zden_re,*/ zden_im; \
double      zres_re, zres_im; \
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
grad_x += b0[i]*(zres_re*cosi[i] - zres_im*sinu[i]); \
grad_y += b0[i]*(zres_im*cosi[i] + zres_re*sinu[i]); \
\
