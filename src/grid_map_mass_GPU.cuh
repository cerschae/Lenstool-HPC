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

@brief: Function for mass computation over a grid

*/

#ifndef GRID_MAP_MASS_GPU_CUH_
#define GRID_MAP_MASS_GPU_CUH_

//#include <fstream>

//#include "gradient2_GPU.cuh"
#include "grid_gradient2_GPU.cuh"
#include "module_cosmodistances.hpp"
#include <structure_hpc.hpp>
//#include <cuda.h>

////type for mapping functions
typedef void (*map_gpu_function_t) (type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t dos, type_t dl, type_t h, type_t z,int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame);

////Map function selection
map_gpu_function_t select_map_mass_function(const struct runmode_param* runmode);

////Map functions
//Mass
void mass_1_grid_CPU_GPU(type_t *,type_t *,type_t *,type_t *,type_t *, type_t , type_t, type_t , type_t , type_t ,int ,int ,const struct grid_param*);
void mass_2_grid_CPU_GPU(type_t *,type_t *,type_t *,type_t *,type_t *, type_t , type_t, type_t , type_t , type_t ,int ,int ,const struct grid_param*);
void mass_3_grid_CPU_GPU(type_t *,type_t *,type_t *,type_t *,type_t *, type_t , type_t, type_t , type_t , type_t ,int ,int ,const struct grid_param*);
void mass_4_grid_CPU_GPU(type_t *,type_t *,type_t *,type_t *,type_t *, type_t , type_t, type_t , type_t , type_t ,int ,int ,const struct grid_param*);

////General mapping function
void map_grid_mass_GPU(map_gpu_function_t mapfunction, type_t *map, const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos ,int nbgridcells,int mode_amp, type_t zm, type_t zs );
void map_grid_mass_GPU(map_gpu_function_t mapfunction, type_t *map,const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int mode_amp,  type_t zm, type_t zs, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart);

//
#endif /* GRADIENTGPU_CUH_ */

