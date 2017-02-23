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
void gradient_grid_GPU_sub(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos,int nbgridcells, int indexactual, int Ncells );
void gradient_grid_GPU_multiple(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells);

__global__ void gradient_grid_kernel(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcells, const struct Potential_SOA *lens);
__global__ void gradient_grid_kernel_multiple(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, int nbgridcells, const struct Potential_SOA *lens, int indexactual, int ncells);

__device__ static double atomicAdd_double(double* address, double val);


__host__ void allocatekernelmemory(double *grid_grad_x_kernel, double *grid_grad_y_kernel, struct grid_param *frame_kernel, struct Potential_SOA *lens_kernel,const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int Ncells);
void freekernelmemory(double *grid_grad_x_kernel, double *grid_grad_y_kernel, struct grid_param *frame_kernel, struct Potential_SOA *lens_kernel);


#endif /* GRADIENTGPU_CUH_ */
