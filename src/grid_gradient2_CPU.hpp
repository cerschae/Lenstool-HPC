/*
 * grid_gradient_CPU.hpp
 *
 *  Created on: Jan 12, 2017
 *      Author: cerschae
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