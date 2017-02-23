/**
* @file   cudafunctions.cuh
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch)
* @date   July 2015
* @version 0,1
* @brief  Header file for cudafunctions.cu
*
*/







//Header guard
#ifndef CUDAFUNCTIONS_CUH
#define CUDAFUNCTIONS_CUH




// Include
//===========================================================================================================
#include <strings.h>
#include "structure.h"
#include <iostream>
#include <cuda_runtime.h>

void cudasafe( cudaError_t error, std::string message);


//End header guard
#endif
