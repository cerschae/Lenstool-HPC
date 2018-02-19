/*
/**
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)
* @date   July 2017
* @version 0,1
*
*/
 */

#ifndef GRID_GRADIENT_GPU_CUH_
#define GRID_GRADIENT_GPU_CUH_

//#include "cudafunctions.cuh"
//#include "gradient_GPU.cuh"
#include <structure_hpc.hpp>

//gradient_grid_GPU_sorted(double*, double*, grid_param const*, Potential_SOA const*, int, int)

//static
//extern "C"
void gradient_grid_GPU_sorted(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, int nbgridcells);
void gradient_grid_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, int nbgridcells);
//
void gradient_grid_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart = 0, int jstart = 0);
//
void module_potentialDerivatives_totalGradient_SOA_CPU_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos);
#if 0

//
//void 
//gradient_grid_GPU_sorted(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, int nbgridcells);

void
module_potentialDerivatives_totalGradient_SOA_CPU_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_cpu, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos);
//
void 
gradient_grid_pinned(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int *Nlens,int nbgridcell);

void 
gradient_grid_pinned_multiple(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int *Nlens,int nbgridcell);

void 
gradient_grid_GPU_sub(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos,int nbgridcells, int indexactual, int Ncells );

#endif

#endif /* GRADIENTGPU_CUH_ */

