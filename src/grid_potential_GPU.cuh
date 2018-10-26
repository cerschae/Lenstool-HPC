/**
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)
* @date   July 2017
* @version 0,1
*
*/
#ifndef GRID_POT_GPU_CUH_
#define GRID_POT_GPU_CUH_

#include <structure_hpc.hpp>
//#include "gradient2_GPU.cuh"


//static
//extern "C"
void potential_grid_GPU(type_t *grid_potential, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart);
void potential_grid_GPU(type_t *grid_potential, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos ,int nbgridcells);

void
module_potential_SOA_CPU_GPU(type_t *grid_potential, const struct grid_param *frame, const struct Potential_SOA *lens_gpu, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart);

#endif
