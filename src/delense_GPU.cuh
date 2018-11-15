#pragma once

#ifdef __WITH_MPI
#include <mpi.h>
#endif

#include <structure_hpc.hpp>
#include <grid_srcplane_conversion.hpp>
#include <grid_gradient_CPU.hpp>
#ifdef __AVX512F__
#include "gradient_avx512f.hpp"
#else
#include <gradient_avx.hpp>
#endif
#include "delense_CPU_utils.hpp"


//
//void delense_barycenter(struct point*** image_pos, int* locimagesfound, int* numimg, const runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, const struct galaxy* sources, double* grid_gradient_x, double* grid_gradient_y);
void delense_barycenter(struct point* image_pos, int* locimagesfound, int* numimg, const runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, const struct galaxy* sources, double* grid_gradient_x, double* grid_gradient_y);
//
//void delense(int nsets, int nimagestot, struct point* image_pos[nsets][nimagestot][MAXIMPERSOURCE], int* locimagesfound[nsets][nimagestot], int* numimg, const runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, const struct galaxy sources[nsets][nimagestot], double* grid_gradient_x, double* grid_gradient_y);
//void delense(int nsets, int nimagestot, struct point image_pos[nsets][nimagestot][MAXIMPERSOURCE], int locimagesfound[nsets][nimagestot], int* numimg, const runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, const struct galaxy sources[nsets][nimagestot], double* grid_gradient_x, double* grid_gradient_y);
//
void delense(struct point* image_pos, int* locimagesfound, int* numimg, const runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, const struct galaxy* sources, double* grid_gradient_x, double* grid_gradient_y);
//
void delense_barycenter_GPU(struct point* image_pos, int* locimagesfound, int* numimg, const runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, const struct galaxy* sources, double* grid_gradient_x, double* grid_gradient_y);
//
void delense_barycenter_GPU_new(struct point* image_pos, int* locimagesfound, int* numimg, const runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, const struct galaxy* sources, double* grid_gradient_x, double* grid_gradient_y);
//
