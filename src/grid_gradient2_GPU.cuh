/**
Lenstool-HPC: HPC based massmodeling software and Lens-map generation
Copyright (C) 2017  Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

@brief: Function for second order derivative computation over a grid

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

void
module_potentialDerivatives_totalGradient2_SOA_CPU_GPU(type_t *grid_grad2_a, type_t *grid_grad_b, type_t *grid_grad2_c, type_t *grid_grad_d,  const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart);


void
module_potentialDerivatives_Kmap_SOA_CPU_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct grid_param *frame, const struct Potential_SOA *lens_gpu, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart);

#endif /* GRADIENTGPU_CUH_ */
