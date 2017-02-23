/**
* @file   cudafunctions.cpp
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch)
* @date   July 2015
* @version 0,1
* @brief  Functions for safe calling of cuda function and timing functions
*
*/






/// Include
///==========================================================================================================
#include<stdio.h>
#include<string>
#include<strings.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<time.h>
#include"cudafunctions.cuh"
#include <cuda_runtime.h>
#include<typeinfo>



/** @brief Safe Cuda call function. Checks for error, prints and then stops the execution.
*/
void cudasafe( cudaError_t error, std::string message)
{
  if(error!=cudaSuccess) {
     fprintf(stderr,"ERROR: %s : %s \n",message.c_str(),cudaGetErrorString(error));
     exit(-1);
  }
}
