#pragma once

#include <structure_hpc.hpp>
#include <gradient_avx.hpp>
#include <grid_gradient_CPU.hpp>
#ifdef __AVX512F__
#include "gradient_avx512f.hpp"
#endif



void PotentialSOAAllocation(Potential_SOA **lens_SOA, const int nhalos);
void PotentialSODeallocation(Potential_SOA *lens_SOA);
void module_readParameters_PotentialSOA_local(std::string infile, Potential_SOA *lens_SOA, int nhalos, int n_tot_halos, cosmo_param cosmology);

