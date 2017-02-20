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


__global__ void gradient_grid_kernel(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcells, const struct Potential_SOA *lens);
__global__ void gradient_grid_piemd_GPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcell, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut);
__global__ void gradient_grid_sis_GPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcell, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut);
__global__ void gradient_grid_piemd_GPU_multiple(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcell, int Ndevice, int indexactual, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut);
__global__ void gradient_grid_kernel_multiple(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, int nbgridcells, const struct Potential_SOA *lens, int indexactual, int ncells);

__device__ static double atomicAdd_double(double* address, double val);

#endif /* GRADIENTGPU_CUH_ */
