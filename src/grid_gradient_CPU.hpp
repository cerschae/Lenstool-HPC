/*
 * grid_gradient_CPU.hpp
 *
 *  Created on: Jan 12, 2017
 *      Author: cerschae
 */

#ifndef GRID_GRADIENT_CPU_HPP_
#define GRID_GRADIENT_CPU_HPP_


#include "structure_hpc.h"
#include <iostream>
#include <string.h>
//#include <cuda_runtime.h>
#include <math.h>
#include <gradient_avx.hpp>
#include <sys/time.h>
#include <fstream>
//#include <immintrin.h>

void gradient_grid_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, int nbgridcells);

static void gradient_grid_general_CPU_old(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, int nbgridcells,  const struct Potential_SOA *lens);
static void gradient_grid_general_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, int nbgridcells,  const struct Potential_SOA *lens);


#endif /* GRID_GRADIENT_CPU_HPP_ */
