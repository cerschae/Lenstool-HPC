/**
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)
* @date   July 2017
* @version 0,1
*
*/

#ifndef GRID_MAP_DPL_GPU_CUH_
#define GRID_MAP_DPL_GPU_CUH_

//#include <fstream>

//#include "gradient2_GPU.cuh"
#include "grid_gradient_GPU.cuh"
#include "module_cosmodistances.hpp"
#include <structure_hpc.hpp>
//#include <cuda.h>

////type for mapping functions
typedef void (*map_gpu_function_t) (type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t dos, type_t dl, type_t h, type_t z,int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame);

////Map function selection
map_gpu_function_t select_map_dpl_function(const struct runmode_param* runmode);

////Map functions
void dpl_grid_CPU_GPU(type_t *map, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame);
////General mapping function
void map_grid_dpl_GPU(map_gpu_function_t mapfunction, type_t *dplx, type_t *dply, const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens,int nbgridcells, int mode_amp, type_t z  );
void map_grid_dpl_GPU(map_gpu_function_t mapfunction, type_t *dplx, type_t *dply, const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, int mode_amp, type_t z, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart = 0, int jstart = 0);

//
#endif /* GRADIENTGPU_CUH_ */

