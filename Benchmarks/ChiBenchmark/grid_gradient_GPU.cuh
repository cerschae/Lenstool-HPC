/*
 * gradientgpu.cuh
 *
 *  Created on: Nov 29, 2016
 *      Author: cerschae
 */

#ifndef GRID_GRADIENT_GPU_CUH_
#define GRID_GRADIENT_GPU_CUH_

#include "cudafunctions.cuh"
#include <structure_hpc.h>

void gradient_grid_GPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int *Nlens,int nbgridcells);
void gradient_grid_pinned(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int *Nlens);
void gradient_grid_pinned_multiple(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int *Nlens);

struct point grad_halo_piemd_SOA_GPU(const struct point *pImage, int Nlens, const struct Potential_SOA *lens);

__global__ void piemd_GPU(double *grad_x_gpu, double * grad_y_gpu, point *pImage, int Nlens, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut);
__global__ void gradient_grid_piemd_GPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut);
__global__ void gradient_grid_sis_GPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut);
__global__ void gradient_grid_piemd_GPU_multiple(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, int Ndevice, int indexactual, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut);

__device__ struct point rotateCoordinateSystem_GPU(struct point P, double theta);
__device__ static double atomicAdd_double(double* address, double val);

#endif /* GRADIENTGPU_CUH_ */
