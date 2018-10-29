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

#ifndef GRID_GRADIENT2_CPU_HPP_
#define GRID_GRADIENT2_CPU_HPP_


#include "structure_hpc.hpp"
#include <iostream>
#include <string.h>
//#include <cuda_runtime.h>
#include <math.h>
#include <gradient_avx.hpp>
#include <gradient2.hpp>
#include <sys/time.h>
#include <fstream>
//#include <immintrin.h>



static void gradient2_grid_general_CPU(matrix *grid_grad2, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart);
static void gradient2_grid_general_CPU(matrix *grid_grad2, const struct grid_param *frame, int Nlens, int nbgridcells,  const struct Potential_SOA *lens, int istart, int jstart);
//
void gradient2_grid_CPU(matrix *grid_grad2, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, int nbgridcells, int istart = 0, int jstart = 0);
void gradient2_grid_CPU(matrix *grid_grad2, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, int nbgridcells_x, int nbgridcells_y, int istart = 0, int jstart = 0);
void gradient2_grid_CPU(matrix *grid_grad2, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart = 0, int jstart = 0);



#endif /* GRID_GRADIENT_CPU_HPP_ */
