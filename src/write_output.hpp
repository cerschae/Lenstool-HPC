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

@brief: Function for writing files

*/





//Header guard
#ifndef WRITE_OUTPUT_H
#define WRITE_OUTPUT_H




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
void write_output_config(std::string infile, char* path);
void write_output_images(std::string  path, const struct runmode_param *runmode, int nimages[], galaxy images[] );

#endif
