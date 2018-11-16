#pragma once

#include <structure_hpc.hpp>
#include <gradient_avx.hpp>
#include <grid_srcplane_conversion.hpp>
#include <grid_gradient_CPU.hpp>
#ifdef __AVX512F__
#include "gradient_avx512f.hpp"
#endif



void PotentialSOAAllocation_GPU(Potential_SOA **lens_SOA, const int nhalos);
void PotentialSODeallocation_GPU(Potential_SOA *lens_SOA);

