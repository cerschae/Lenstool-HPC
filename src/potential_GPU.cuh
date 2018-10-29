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

@brief: Function for single deflection potential computation

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
