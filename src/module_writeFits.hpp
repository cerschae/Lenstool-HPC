/**
Lenstool-HPC: HPC based massmodeling software and Lens-map generation
Copyright (C) 2017  Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

@brief: Function for writing fitsfiles

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
void module_writeFits(std::string path, std::string filename, int ii, type_t *map, const struct runmode_param* runmode, int size, const struct grid_param* frame, type_t ra, type_t dec );
void module_writeFits(std::string path, std::string filename, type_t *map, const struct runmode_param* runmode, int size, const struct grid_param* frame, type_t ra, type_t dec );

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
