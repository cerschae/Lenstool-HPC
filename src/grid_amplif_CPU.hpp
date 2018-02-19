/**
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)
* @date   July 2017
* @version 0,1
*
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



void amplif_grid_CPU(type_t *amplif,const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells,int mode_amp,type_t z, int istart = 0, int jstart = 0);
//
void amplif_5_grid_general_CPU(type_t *amplif,const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y,type_t z, int istart, int jstart);

#endif /* GRID_GRADIENT_CPU_HPP_ */
