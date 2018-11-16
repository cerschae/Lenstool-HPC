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

@brief: Functions for writing files

*/

/// Include
///==========================================================================================================
#include <cstring>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "write_output.hpp"
#include <iomanip>


/// Functions
///==========================================================================================================

/** @brief copy the input file and save the copy in the output folder
 *
 * @param infile path to config file
 * @param path of result folder
 *
 */

void write_output_config(std::string infile, char* path){
	std::string line1, input_file(path);
	input_file.append("/configuration.par");
	std::ofstream COPY(input_file.c_str());  // Output stream
	std::ifstream IN(infile.c_str(),std::ios::in); // open file
    	if ( IN )
    	{
		if ( COPY )
		{
			while(std::getline(IN,line1))
			COPY << line1 << std::endl;
		}
	}
	else
	{
		fprintf(stderr, "ERROR: Input parameter file %s not found\n", infile.c_str());
        	exit(-1);
	}

}

/** @brief writing the predicted images to predicted_images.txt
 * print the position of the images of the source on the lens plane
 *
 * @param path of result folder
 * @param runmode runmode parameters
 * @param nimages number of images per source
 * @param images images
 */

void write_output_images(char* path, const struct runmode_param *runmode, int nimages[], galaxy images[] ){
	std::string line1, input_file(path);
	input_file.append("/predicted_images.txt");
	std::ofstream OUT(input_file.c_str());  // Output stream
	int index =0;
    if ( OUT ){

		OUT << "#Images of sources lensed unto the lens plane " << std::endl;
		OUT << "#Corresponding Source, Center_x, Center_y, Shape_a, Shape_b, Shape_ellipticity, Redshift" << std::endl;
		OUT << "#[],\t[arcsec],\t[arcsec],\t[arcsec],\t[arcsec],\t[rad],\t[]" << std::endl;
		OUT << std::setprecision(7);
		for(int i=0; i < runmode->nsets; i++){
			//printf("nimages %d ", nimages[i]);
			for(int j = 0; j < nimages[i]; j++){
				index = i * MAXIMPERSOURCE + j;
				OUT << i+1 << " " << images[index].center.x << " " << images[index].center.y << " " << images[index].shape.a << " " << images[index].shape.b << " " << images[index].shape.theta << " " << images[index].redshift <<  " " << images[index].mag << std::endl;
			}
		}
	}
}

