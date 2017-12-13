/*
 * gradientgpu.cuh
 *
 *  Created on: Nov 29, 2016
 *      Author: cerschae
 */

#ifndef GRID_GRADIENT2_GPU_CUH_
#define GRID_GRADIENT2_GPU_CUH_

#include <structure_hpc.hpp>
//#include "gradient2_GPU.cuh"


//static
//extern "C"
//void gradient2_grid_GPU_sorted(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, int nbgridcells);
void gradient2_grid_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, int nbgridcells);
//
void gradient2_grid_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart = 0, int jstart = 0);
//
void module_potentialDerivatives_totalGradient2_SOA_CPU_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct grid_param *frame, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos);

#endif /* GRADIENTGPU_CUH_ */
