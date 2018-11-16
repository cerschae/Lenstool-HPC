#pragma once
#ifndef __CHI_CPU_HPP__
#define __CHI_CPU_HPP__

#include <structure_hpc.hpp>
#include <gradient_avx.hpp>
#include <grid_srcplane_conversion.hpp>
#include <grid_gradient_CPU.hpp>
#ifdef __AVX512F__
#include "gradient_avx512f.hpp"
#endif

void mychi_bruteforce_SOA_CPU_grid_gradient_barycentersource(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images);

#endif
