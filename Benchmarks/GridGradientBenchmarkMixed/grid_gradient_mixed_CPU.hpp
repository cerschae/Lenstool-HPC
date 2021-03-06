/*
 * grid_gradient_CPU.hpp
 *
 *  Created on: Jan 12, 2017
 *      Author: cerschae
 */

#ifndef GRID_GRADIENT_CPU_HPP_
#define GRID_GRADIENT_CPU_HPP_


#include "structure_hpc.hpp"
#include <iostream>
#include <string.h>
//#include <cuda_runtime.h>
#include <math.h>
#include <gradient_avx.hpp>
#include <sys/time.h>
#include <fstream>
//#include <immintrin.h>
//
//
struct point_double
{
        type_t x;
        type_t y;
};
//
//
//
void gradient_grid_CPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, int nbgridcells, int istart = 0, int jstart = 0);
//
static void gradient_grid_general_CPU_old(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, int Nlens, int nbgridcells,  const struct Potential_SOA *lens);
//static void gradient_grid_general_mixed_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, int nbgridcells,  const struct Potential_SOA *lens, int istart, int jstart);
void gradient_grid_general_mixed_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, int nbgridcells, const struct Potential_SOA *lens, int istart, int jstart);
void gradient_grid_general_double_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, int nbgridcells, const struct Potential_SOA *lens, int istart, int jstart);




#endif /* GRID_GRADIENT_CPU_HPP_ */
