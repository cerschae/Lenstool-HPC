/**
* @file   module_writeFits.h
* @Author Markus Rexroth, EPFL
* @date   July 2015
* @version 0,1
* @brief  Header file for module_writeFits.c
*/







//Header guard
#ifndef MODULE_WRITEFITS_CUH
#define MODULE_WRITEFITS_CUH




// Include
//===========================================================================================================
#include <stdio.h>
#include <structure_hpc.hpp>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>




// Function declarations
//===========================================================================================================
void module_writeFits(std::string path, std::string filename, int ii, type_t *map, const struct runmode_param* runmode, const struct grid_param* frame, type_t ra, type_t dec );

// We are calling C functions here, but we are calling them with a cuda/c++ compiler, thus we need to declare them as extern "C"
#ifdef __cplusplus
extern "C" {
#endif
int module_writeFits_Image(char *filename, double *ima, int nx,int ny, double xmin,double xmax,double ymin,double ymax);
#ifdef __cplusplus
}
extern "C" {
#endif
int module_writeFits_ImageAbsoluteCoordinates(char *filename, double *ima, int nx,int ny, double xmin,double xmax,double ymin,double ymax, double ra,double dec);
#ifdef __cplusplus
}
extern "C" {
#endif
int module_writeFits_cube(char *filename, double ***cube, int nx,int ny, int nz, double xmin,double xmax,double ymin,double ymax, double lmin, double lmax);
#ifdef __cplusplus
}
#endif

#endif
