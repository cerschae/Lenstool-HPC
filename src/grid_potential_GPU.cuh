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

@brief: Function for deflection computation over a grid

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
