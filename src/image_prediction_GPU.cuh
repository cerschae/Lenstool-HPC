#pragma once
#ifndef __IMG_PRED_GPU__
#define __IMG_PRED_GPU__

#include <structure_hpc.hpp>

#include "delense_CPU_utils.hpp"
#include "delense_CPU.hpp"
#include "delense_GPU.cuh"
#include "grid_gradient_GPU.cuh"

void image_prediction(struct galaxy *predicted_images, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nImagesSet, struct galaxy *images);

//void mychi_bruteforce_SOA_GPU_grid_gradient_barycentersource(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images);

#endif
