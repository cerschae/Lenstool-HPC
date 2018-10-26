/**
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)
* @date   July 2017
* @version 0,1
*
*/
#ifndef POTENTIAL_GPU_CUH_
#define POTENTIAL_GPU_CUH_

//#include "cudafunctions.cuh"
#include <fstream>
#include <structure_hpc.hpp>
#include <lt_math_GPU.cuh>
//#include <cuda.h>


//
//__global__ void module_potentialDerivatives_totalGradient2_SOA_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct Potential_SOA *lens, const struct grid_param *frame, int nhalos, int nbgridcells);
//
//__global__ void module_potentialDerivatives_totalGradient2_SOA_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct Potential_SOA *lens, const struct grid_param *frame, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart = 0, int jstart = 0);

__global__
void
module_potential_totalPotential_SOA_GPU(type_t *potential_GPU, const struct Potential_SOA *lens, const struct grid_param *frame, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart);
#endif /* POTENTIAL_GPU_CUH_ */
