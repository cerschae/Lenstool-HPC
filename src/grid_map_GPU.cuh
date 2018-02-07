/*
 * gradientgpu.cuh
 *
 *  Created on: Nov 29, 2016
 *      Author: cerschae
 */

#ifndef GRID_MAP_GPU_CUH_
#define GRID_MAP_GPU_CUH_

#include <structure_hpc.hpp>
#include "grid_gradient2_GPU.cuh"
#include "module_cosmodistances.hpp"

//type for mapping functions
typedef void (*map_gpu_function_t) (type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t z,int mode_amp,int nhalos,int nbgridcells_x, int nbgridcells_y);
// Amplification
//map_gpu_function_t amplif_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t z,int mode_amp,int nhalos,int nbgridcells_x, int nbgridcells_y);
void amplif_grid_CPU_GPU(type_t *,type_t *,type_t *,type_t *,type_t *, type_t , type_t ,int ,int ,int , int );

//General mapping function
void map_grid_GPU(map_gpu_function_t mapfunction, type_t *amplif,const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens,int nbgridcells, int mode_amp, type_t z  );
void map_grid_GPU(map_gpu_function_t mapfunction, type_t *amplif,const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, int mode_amp, type_t z, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart = 0, int jstart = 0);
//
#endif /* GRADIENTGPU_CUH_ */

