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
#include "chi_CPU.hpp"
//#ifdef __WITH_GPU
#include "grid_gradient_GPU.cuh"
////#endif
//#include "gradient_GPU.cuh"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __WITH_MPI
#warning "MPI enabled"
#include <mpi.h>
#include "mpi_check.h"
#include "chi_comm.hpp"
#endif

#include "delense_CPU_utils.hpp"
#include "delense_CPU.hpp"
#include "delense_GPU.cuh"

#include "chi_computation.hpp"
#ifdef __WITH_GPU
#include <cuda_runtime.h>
#include "cudafunctions.cuh"
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


//
//
//
double myseconds()
{
	struct timeval  tp;
	struct timezone tzp;
	int i;

	i = gettimeofday(&tp,&tzp);
	return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

void mychi_bruteforce_SOA_CPU_grid_gradient_barycentersource(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images)
{
	//
	//double dx, dy; //, x_pos, y_pos;        //pixelsize
	//
	//	double im_dist[MAXIMPERSOURCE]; // distance to a real image for an "theoretical" image found from a theoretical source
	//	int im_index;       		// index of the closest real image for an image found from a theoretical source
	//	int second_closest_id;   	// index of the second closest real image for an image found from a theoretical source
	//	int thread_found_image = 0;   	// at each iteration in the lens plane, turn to 1 whether the thread find an image
	// 	struct point im_position, temp; // position of the image found from a theoretical source + temp variable for comparison
	//	struct triplet Tsup, Tinf, Tsupsource, Tinfsource;// triangles for the computation of the images created by the theoretical sources
	//	
	unsigned int nsets = runmode->nsets;
	//
	struct galaxy sources[nsets]; // theoretical sources (common for a set)
	int    nimagesfound  [nsets][MAXIMPERSOURCE]; // number of images found from the theoretical sources
	int    locimagesfound[nsets];
	//int*    locimagesfound = (int*) malloc(nsets*sizeof(int));
	struct point image_pos [nsets][MAXIMPERSOURCE];
	//struct point* image_pos = (struct point*) malloc(nsets*MAXIMPERSOURCE*sizeof(struct point));	
	//double image_dist    [nsets][runmode->nimagestot][MAXIMPERSOURCE];
	struct point tim     [MAXIMPERSOURCE]; // theoretical images (computed from sources)
	//
	int world_size = 1;
	int world_rank = 0;
#ifdef __WITH_MPI
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	unsigned int verbose = (world_rank == 0);
	//
	int grid_size     = runmode->nbgridcells;
	int loc_grid_size = runmode->nbgridcells/world_size;
	//
	double y_len      = fabs(frame->ymax - frame->ymin);
	int    y_len_loc  = runmode->nbgridcells/world_size;
	int    y_pos_loc  = (int) world_rank*y_len_loc;
	int    y_bound    = y_len_loc;
	//
	if ((world_rank + 1) != world_size) y_bound++;
	//
	//
	const double dx   = (frame->xmax - frame->xmin)/(runmode->nbgridcells - 1);
	const double dy   = (frame->ymax - frame->ymin)/(runmode->nbgridcells - 1);
	//
	double *grid_gradient_x, *grid_gradient_y;
	//
	//grid_gradient_x   = (double *)malloc((int) grid_size*loc_grid_size*sizeof(double));
	//grid_gradient_y   = (double *)malloc((int) grid_size*loc_grid_size*sizeof(double));
	//grid_gradient_x   = (double *)malloc((int) grid_size*y_bound*sizeof(double));
	//grid_gradient_y   = (double *)malloc((int) grid_size*y_bound*sizeof(double));
#ifdef __WITH_GPU
        cudaMallocManaged(&grid_gradient_x, grid_size*y_bound*sizeof(double));
	cudasafe(cudaGetLastError(), "grid_gradient_x allocation");
        cudaMallocManaged(&grid_gradient_y, grid_size*y_bound*sizeof(double));
	cudasafe(cudaGetLastError(), "grid_gradient_y allocation");
#else
        grid_gradient_x   = (double *)malloc((int) grid_size*y_bound*sizeof(double));
        grid_gradient_y   = (double *)malloc((int) grid_size*y_bound*sizeof(double));
#endif
	//
	const int grid_dim = runmode->nbgridcells;
	//Packaging the image to sourceplane conversion
	double time = -myseconds();
	double grad_time = -myseconds();
	//gradient_grid_CPU(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, grid_dim);
#ifdef __WITH_GPU
	//@@printf("Using GPUs...\n"); fflush(stdout);
#ifdef __WITH_UM
	//@@printf("Unified memory...\n"); fflush(stdout);
	gradient_grid_GPU_UM(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, dx, dy, grid_size, y_bound, 0, y_pos_loc);
#else
	gradient_grid_GPU(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, dx, dy, grid_size, y_bound, 0, y_pos_loc);
#endif
	//gradient_grid_GPU(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, grid_size);
#else
	int zero = 0;
	gradient_grid_CPU(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, dx, dy, grid_size, y_bound, zero, y_pos_loc);
#endif
	//gradient_grid_CPU(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, grid_size, loc_grid_size);
#ifdef __WITH_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	grad_time += myseconds();
	//
	int index             = 0;       // index tracks the image within the total image array in the image plane
	*chi  		      = 0;
	//
	double chi2_time      = 0.;
	double delense_time   = 0.;
	double image_time     = 0.;
	//
	int images_found      = 0;
	long int images_total = 0;
	//
	//printf("@@nsets = %d nx = %d ny = %d, xmin = %f, dx = %f, ymin = %f, dy = %f\n", nsets, runmode->nbgridcells, runmode->nbgridcells, frame->xmin, dx, frame->ymin, dy );
	//
	int numsets = 0;
	//
	for( int  source_id = 0; source_id < nsets; source_id ++)
		numsets += nimages_strongLensing[source_id];
	//
	time = -myseconds();
	//
	// nsets     : number of images in the source plane
	// nimagestot: number of images in the image plane
	//
	//
	// Image Lensing
	//
	image_time -= myseconds();
	//
	//unsigned int nsets = nsets;
	//
	for( int  source_id = 0; source_id < nsets; source_id ++)
	{
		// number of images in the image plane for the specific image (1,3,5...)
		unsigned short int nimages = nimages_strongLensing[source_id];
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
//#if 1
		index += nimages;
	}
#ifdef __WITH_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	image_time += myseconds();
	//
	// Delensing
	//
	delense_time -= myseconds();
	index 	      = 0;
	int numimg    = 0;
#ifndef __WITH_GPU
	delense_barycenter(&image_pos[0][0], &locimagesfound[0], &numimg, runmode, lens, frame, nimages_strongLensing, sources, grid_gradient_x, grid_gradient_y);	
#else
	//delense_barycenter_GPU(&image_pos[0][0], &locimagesfound[0], &numimg, runmode, lens, frame, nimages_strongLensing, sources, grid_gradient_x, grid_gradient_y);
	//struct point* ip[MAXIMPERSOURCE] = image_pos;	
	delense_barycenter_GPU(&image_pos[0][0], &locimagesfound[0], &numimg, runmode, lens, frame, nimages_strongLensing, sources, grid_gradient_x, grid_gradient_y);
	printf("done. numimg = %d\n", numimg);fflush(stdout);
#endif
	//
	// Communications
	//
#ifdef __WITH_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	delense_time += myseconds();
	//
	double comm_time = -myseconds();
	//
	int          numimagesfound    [nsets];
	struct point imagesposition    [nsets][MAXIMPERSOURCE];
	//
	memset(&numimagesfound, 0, nsets*sizeof(int));
	memset(&imagesposition, 0, nsets*MAXIMPERSOURCE*sizeof(point));
	//
#ifdef __WITH_MPI
	delense_comm(&numimagesfound[0], &imagesposition[0][0],  &numimg, nimages_strongLensing, &locimagesfound[0], &image_pos[0][0], nsets, world_rank, world_size);
#endif
	//
	comm_time += myseconds();
	//
	// image extraction to compute
	//
	chi2_time = -myseconds();
	if (verbose)
		chi_computation(chi, &images_found, &numimagesfound[0], &imagesposition[0][0], nimages_strongLensing, images, nsets);

#ifdef __WITH_GPU
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	//
	chi2_time += myseconds();
	time      += myseconds();
	//
	if (verbose)
	{
		//
		//		int nthreads = 1;
		//
		//#pragma omp parallel
		//		nthreads = omp_get_num_threads();
		//
		printf("	overall time  = %f s.\n", time);
		printf("		- grad    time = %f s.\n", grad_time);
		printf("		- image   time = %f s.\n", image_time);
		printf("		- delense time = %f s.\n", delense_time);
		printf("		- comm    time = %f s.\n", comm_time);
		printf("		- chi2    time = %f s.\n", chi2_time);
		//
		printf("	images found: %d out of %ld\n", images_found, images_total);
	}
	//
#ifdef __WITH_GPU
        cudaFree(grid_gradient_x);
        cudaFree(grid_gradient_y);
#else
        free(grid_gradient_x);
        free(grid_gradient_y);
#endif
}
