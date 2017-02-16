/*
 * gradient_GPU.cuh
 *
 *  Created on: Feb 1, 2017
 *      Author: cerschae
 */

#ifndef GRADIENT_GPU_CUH_
#define GRADIENT_GPU_CUH_

#include "cudafunctions.cuh"
#include <fstream>
#include <structure_hpc.h>

__device__ struct point module_potentialDerivatives_totalGradient_SOA_GPU(const struct point *pImage, const struct Potential_SOA *lens, int nhalos);

__device__ inline struct point rotateCoordinateSystem_GPU(struct point P, double theta);
__device__ inline struct point rotateCoordinateSystem_GPU_2(struct point P, double cosi, double sinu);

#endif /* GRADIENT_GPU_CUH_ */
