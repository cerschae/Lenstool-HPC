#pragma once
#ifndef __CHI_GPU_CUH__
#define __CHI_GPU_CUH__

#include <structure_hpc.h>
//#include <gradient_avx.hpp>
//#include <grid_srcplane_conversion.hpp>
//#include <grid_gradient_CPU.hpp>
#include <grid_gradient_GPU.cuh>
#include <cuda_runtime.h>
//#ifdef __AVX512F__
//#include "gradient_avx512f.hpp"
//#endif

void chi_bruteforce_SOA_GPU_grid_gradient(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images);
void chi_bruteforce_GPU_CPU(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images);

#endif
