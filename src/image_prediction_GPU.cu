#include <iostream>
#include <string.h>
//#include <cuda_runtime.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <unistd.h>
#include <assert.h>


#ifndef __xlC__
#warning "gnu compilers"
#include <immintrin.h>
#endif

//#include "simd_math.h"
#include "image_prediction_GPU.cuh"

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __WITH_MPI
#warning "MPI enabled"
#include <mpi.h>
#include "mpi_check.h"
#include "chi_comm.hpp"
#endif



//#include "chi_computation.hpp"
#ifdef __WITH_GPU
#include <cuda_runtime.h>
#include "cudafunctions.cuh"
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


//returns predicted images for every constraint (images) given
void image_prediction(struct galaxy *predicted_images, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nImagesSet,struct galaxy *images){

	struct galaxy sources[runmode->nsets];
	point predicted_images_pos[runmode->nsets][MAXIMPERSOURCE];

	double *grid_gradient_x, *grid_gradient_y;

	int grid_size     = runmode->nbgridcells;
	const double dx   = (frame->xmax - frame->xmin)/(runmode->nbgridcells - 1);
	const double dy   = (frame->ymax - frame->ymin)/(runmode->nbgridcells - 1);
	int index             = 0;

//==============Calculating the gradient over the whole grid =====================

	//Allocating memory for gradient
#ifdef __WITH_GPU
	cudaMallocManaged(&grid_gradient_x, grid_size*grid_size*sizeof(double));
	cudasafe(cudaGetLastError(), "grid_gradient_x allocation");
	cudaMallocManaged(&grid_gradient_y, grid_size*grid_size*sizeof(double));
	cudasafe(cudaGetLastError(), "grid_gradient_y allocation");
#else
	grid_gradient_x   = (double *)malloc((int) grid_size*grid_size*sizeof(double));
	grid_gradient_y   = (double *)malloc((int) grid_size*grid_size*sizeof(double));
#endif

	//Computing Gradient
#ifdef __WITH_GPU
	//@@printf("Using GPUs...\n"); fflush(stdout);
#ifdef __WITH_UM
	//@@printf("Unified memory...\n"); fflush(stdout);
	gradient_grid_GPU_UM(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, dx, dy, grid_size, grid_size, 0, 0);
#else
	gradient_grid_GPU(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, dx, dy, grid_size, grid_size, 0, 0);
#endif
	//gradient_grid_GPU(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, grid_size);
#else
	int zero = 0;
	gradient_grid_CPU(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, runmode->nbgridcells);
#endif

	int numsets = 0;
	for( int  source_id = 0; source_id < runmode->nsets; source_id ++)
		numsets += nImagesSet[source_id];
	//
	for( int  source_id = 0; source_id < runmode->nsets; source_id ++)
	{
		// number of images in the image plane for the specific image (1,3,5...)
		unsigned short int nimages = nImagesSet[source_id];
		//@@printf("@@ source_id = %d, nimages = %d\n",  source_id, nimages_strongLensing[source_id]);
		//Init sources
		sources[source_id].center.x = sources[source_id].center.y = 0.;
		//____________________________ image (constrains) loop ________________________________
		for(unsigned short int image_id = 0; image_id < nimages; image_id++)
		{
			//
			struct galaxy sources_temp;
			struct point Grad = module_potentialDerivatives_totalGradient_SOA(&images[index + image_id].center, lens, runmode->nhalos);
			//
			// find the position of the constrain in the source plane
			sources_temp.center.x = images[index + image_id].center.x - images[index + image_id].dr*Grad.x;
			sources_temp.center.y = images[index + image_id].center.y - images[index + image_id].dr*Grad.y;

			//________ Adding up for barycenter comp _________
			sources[source_id].center.x += sources_temp.center.x;
			sources[source_id].center.y += sources_temp.center.y;
			//time += myseconds();
			sources[source_id].redshift = 0.;
			sources[source_id].shape.a = sources[source_id].shape.b = sources[source_id].shape.theta = (type_t) 0.;
			sources[source_id].mag = 0.;
			//
		}
		//________ Dividing by nimages for final barycenter result _________
		sources[source_id].center.x /= nimages;
		sources[source_id].center.y /= nimages;

		sources[source_id].redshift = images[index].redshift;
		sources[source_id].dr       = images[index].dr;
		sources[source_id].dls      = images[index].dls;
		sources[source_id].dos      = images[index].dos;
		index += nimages;
	}

	int locimagesfound[runmode->nsets];
	int numimg    = 0;
	delense_barycenter(&predicted_images_pos[0][0], &locimagesfound[0], &numimg, runmode, lens, frame, nImagesSet, sources, grid_gradient_x, grid_gradient_y);


	std::cout << "Images:" << std::endl;
	for(int  source_id = 0; source_id < runmode->nsets; source_id++){
		std::cout << "Sources:" << source_id << ": " << sources[source_id].center.x << " " << sources[source_id].center.y << std::endl;
		unsigned short int nimages = nImagesSet[source_id];
		for(unsigned short int image_id = 0; image_id < nimages; image_id++){
			std::cout << "Images:" << predicted_images_pos[source_id][image_id].x << " " << predicted_images_pos[source_id][image_id].y << " " <<  std::endl;
			predicted_images[source_id * MAXIMPERSOURCE + image_id].center.x = predicted_images_pos[source_id][image_id].x;
			predicted_images[source_id * MAXIMPERSOURCE + image_id].center.y = predicted_images_pos[source_id][image_id].y;
			predicted_images[source_id * MAXIMPERSOURCE + image_id].redshift = sources[source_id].redshift;
			predicted_images[source_id * MAXIMPERSOURCE + image_id].shape.a = predicted_images[source_id * MAXIMPERSOURCE + image_id].shape.b = predicted_images[source_id * MAXIMPERSOURCE + image_id].shape.theta = (type_t) 0.;
			predicted_images[source_id * MAXIMPERSOURCE + image_id].mag = 0.;

		}
	}


}
