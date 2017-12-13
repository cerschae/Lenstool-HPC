/*
 * grid_gradient_CPU.hpp
 *
 *  Created on: Jan 12, 2017
 *      Author: cerschae
 */

#ifndef GRID_AMPLIF_CPU_HPP_
#define GRID_AMPLIF_CPU_HPP_


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



void amplif_grid_CPU(type_t *amplif, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells,int mode_amp, int istart = 0, int jstart = 0);
//
void amplif_5_grid_general_CPU(type_t *amplif, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart);

#endif /* GRID_GRADIENT_CPU_HPP_ */