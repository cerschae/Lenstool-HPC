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
#include<mpi.h>
#include"mpi_check.h"
#endif


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


double myseconds()
{
        struct timeval  tp;
        struct timezone tzp;
        int i;

        i = gettimeofday(&tp,&tzp);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}


void mychi_bruteforce_SOA_CPU_grid_gradient(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images)
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
	struct galaxy sources[runmode->nsets][runmode->nimagestot]; // theoretical sources (common for a set)
	int    nimagesfound  [runmode->nsets][runmode->nimagestot][MAXIMPERSOURCE]; // number of images found from the theoretical sources
	int    locimagesfound[runmode->nsets][runmode->nimagestot];
	struct point image_pos [runmode->nsets][runmode->nimagestot][MAXIMPERSOURCE];
	//double image_dist    [runmode->nsets][runmode->nimagestot][MAXIMPERSOURCE];
	struct point tim     [runmode->nimagestot][MAXIMPERSOURCE]; // theoretical images (computed from sources)
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
	grid_gradient_x   = (double *)malloc((int) grid_size*y_bound*sizeof(double));
	grid_gradient_y   = (double *)malloc((int) grid_size*y_bound*sizeof(double));
	//
	//uais
        //if (verbose) printf("@@%d: nsets = %d nimagestot = %d maximgpersource = %d, grid_size = %d, loc_grid_size = %d, y_pos_loc = %d\n", world_rank, runmode->nsets, runmode->nimagestot, MAXIMPERSOURCE, grid_size, loc_grid_size, y_pos_loc);
	//
	const int grid_dim = runmode->nbgridcells;
	//Packaging the image to sourceplane conversion
	double time = -myseconds();
	//gradient_grid_CPU(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, grid_dim);
#ifdef __WITH_GPU
	gradient_grid_GPU(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, dx, dy, grid_size, y_bound, 0, y_pos_loc);
	//gradient_grid_GPU(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, grid_size);
#else
	gradient_grid_CPU(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, dx, dy, grid_size, y_bound, 0, y_pos_loc);
#endif

#if 0
	if (world_rank == 0)
	for (int jj = 0; jj < loc_grid_size; ++jj)
		for (int ii = 0; ii < runmode->nbgridcells; ++ii)
			printf("%d %d: %f %f\n", ii, jj, grid_gradient_x[jj*grid_size + ii], grid_gradient_y[jj*grid_size + ii]);
#endif

	//gradient_grid_CPU(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, grid_size, loc_grid_size);
	time += myseconds();
	if (verbose) printf("	Gridgrad time = %f s.\n", time);
	//
	int index          = 0;       // index tracks the image within the total image array in the image plane
	*chi  		   = 0;
	//
	double chi2_time   = 0.;
	double loop_time   = 0.;
	double image_time  = 0.;
	//
	int images_found   = 0;
	long int images_total   = 0;
	//
        //printf("@@nsets = %d nx = %d ny = %d, xmin = %f, dx = %f, ymin = %f, dy = %f\n", runmode->nsets, runmode->nbgridcells, runmode->nbgridcells, frame->xmin, dx, frame->ymin, dy );
	//
	int numsets = 0;
	//
	for( int  source_id = 0; source_id < runmode->nsets; source_id ++)
		numsets += nimages_strongLensing[source_id];
	//printf("@@Total numsets = %d\n", numsets);
	//
	time = -myseconds();
	//
	// nsets     : number of images in the source plane
	// nimagestot: number of images in the image plane
	//
	image_time -= myseconds();
	//
	unsigned int nsets = runmode->nsets;
	for( int  source_id = 0; source_id < runmode->nsets; source_id ++)
	{
		// number of images in the image plane for the specific image (1,3,5...)
		unsigned short int nimages = nimages_strongLensing[source_id];
		//@@printf("@@ source_id = %d, nimages = %d\n",  source_id, nimages_strongLensing[source_id]);
		//____________________________ image (constrains) loop ________________________________
		for(unsigned short int image_id = 0; image_id < nimages; image_id++)
		{
			//printf("@@  nimages = %d\n",  nimages_strongLensing[source_id]);
			//________ computation of theoretical sources _________
			// output: sources[source_id].center
			//printf("Image = %f %f\n", images[index + image_id].center.x, images[index + image_id].center.y);
			mychi_transformImageToSourcePlane_SOA(runmode->nhalos, 
							&images[index + image_id].center,
							 images[index + image_id].dr, 
							 lens, 
							&sources[source_id][image_id].center);
			//
			// void mychi_transformImageToSourcePlane_SOA(const int Nlens, const struct point *image_point, double dlsds, const struct Potential_SOA *lens, struct point *source_point)
			//
			struct point Grad = module_potentialDerivatives_totalGradient_SOA(&images[index + image_id].center, lens, runmode->nhalos);
			//printf("	image %d, %d (%d) = (%.15f, %.15f) -> (%.15f, %.15f)\n", source_id, image_id, nimages, images[index + image_id].center.x, images[index + image_id].center.y, sources[source_id].center.x, sources[source_id].center.y );
			//printf("Grad.x = %f, %f\n", Grad.x, grid_gradient_x[images[index + image_id].center.x/dx]); 
			//
			// find the position of the constrain in the source plane
			sources[source_id][image_id].center.x = images[index + image_id].center.x - images[index + image_id].dr*Grad.x;
                        sources[source_id][image_id].center.y = images[index + image_id].center.y - images[index + image_id].dr*Grad.y;
			//
			//time += myseconds(); 
			//
			sources[source_id][image_id].redshift = images[index + image_id].redshift;
			//
			sources[source_id][image_id].dr       = images[index + image_id].dr;
			sources[source_id][image_id].dls      = images[index + image_id].dls;
			sources[source_id][image_id].dos      = images[index + image_id].dos;
#if 1
		}
		index += nimages; 
	}
#ifdef __WITH_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
	image_time += myseconds();
	//
	// main loop
	//
	//double y_len     = fabs(frame->ymax - frame->ymin);
	//int    y_len_loc = runmode->nbgridcells/world_size; 
	//int    y_pos_loc = (int) world_rank*y_len_loc;
	//printf("%d out of %d: y_id = %d to %d\n", world_rank, world_size, y_pos_loc, y_pos_loc + y_len_loc - 1);
	//fflush(stdout);
        //MPI_Barrier(MPI_COMM_WORLD);
	//
	loop_time    -= myseconds();
	index 	      = 0;
	int numimg    = 0;
	//
        for( int  source_id = 0; source_id < runmode->nsets; source_id ++)
        {
		// number of images in the image plane for the specific image (1,3,5...)
		unsigned short int nimages = nimages_strongLensing[source_id];
		//printf("@@ source_id = %d, nimages = %d\n",  source_id, nimages_strongLensing[source_id]);
		//____________________________ image (constrains) loop ________________________________
		for(unsigned short int image_id = 0; image_id < nimages; image_id++)
		{
#endif
			//
			//struct point image_pos [MAXIMPERSOURCE];
			//
			//MPI_Barrier(MPI_COMM_WORLD);
			//if (verbose) printf("source = %d, image = %d\n", source_id, image_id);
			//if (verbose) fflush(stdout);	
			int loc_images_found = 0;
#ifdef __WITH_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif	
#pragma omp parallel 
#pragma omp for reduction(+: images_total) 
			//for (int y_id = 0; y_id < (runmode->nbgridcells - 1); ++y_id )
			//for (int y_id = 0; y_id < (y_len_loc - 1); ++y_id)
			for (int y_id = 0; y_id < (y_bound - 1); ++y_id)
			//for (int y_id = world_rank*y_pos_loc; y_id < (world_rank*y_pos_loc + y_len_loc - 1); ++y_id)
			{
				//for (int y_id = 0; (y_id < runmode->nbgridcells - 1) /*&& (loc_images_found != nimages)*/; ++y_id )
				for (int x_id = 0; x_id < runmode->nbgridcells - 1 ; ++x_id)
				{
					//int yy_pos = MIN(y_pos_loc + y_id, runmode->nbgridcells - 1);
					images_total++;
					//
					double x_pos = frame->xmin + (            x_id)*dx;
					double y_pos = frame->ymin + (y_pos_loc + y_id)*dy;
					//printf("%d: x_id = %d xpos = %f, y_id = %d ypos = %f\n", world_rank, x_id, x_pos, y_pos_loc + y_id, y_pos);
					//
					// Define the upper + lower triangle, both together = square = pixel
					// 
					struct triplet Tsup, Tinf;
					//
					Tsup.a.x = x_pos;
					Tsup.b.x = x_pos;
					Tsup.c.x = x_pos + dx;
					Tinf.a.x = x_pos + dx;
					Tinf.b.x = x_pos + dx;
					Tinf.c.x = x_pos;
					//
					Tsup.a.y = y_pos;
					Tsup.b.y = y_pos + dy;
					Tsup.c.y = y_pos;
					Tinf.a.y = y_pos + dy;
					Tinf.b.y = y_pos;
					Tinf.c.y = y_pos + dy;
					//
					// Lens to Sourceplane conversion of triangles
					// 
					// double time = -myseconds();
					//
					struct triplet Timage;
					struct triplet Tsource;
					//
					//printf("	- Tsup = %f %f\n", Tsup.a.x, Tsup.a.y);
					mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_upper(&Tsup, sources[source_id][image_id].dr, &Tsource, grid_gradient_x, grid_gradient_y, y_id*grid_dim + x_id, grid_dim);
					//mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_upper(&Tsup, sources[source_id][image_id].dr, &Tsource, grid_gradient_x, grid_gradient_y, (y_pos_loc + y_id)*grid_dim + x_id, grid_dim);
					//@@if (world_rank == 1)
					//
					int thread_found_image = 0;
					//
					if (mychi_insideborder(&sources[source_id][image_id].center, &Tsource) == 1)
					{
						//printf("	-> %d found: y_id = %d, x_id = %d : pixel %f %f -> %f %f\n", world_rank, y_pos_loc + y_id, x_id, Tsup.a.x, Tsup.a.y, Tsource.a.x, Tsource.a.y);
						thread_found_image = 1;
						Timage = Tsup;
					}
					else 
					{
						//printf("	- Tinf = %f %f\n", Tinf.a.x, Tinf.a.y);
						mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_lower(&Tinf, sources[source_id][image_id].dr, &Tsource, grid_gradient_x, grid_gradient_y, y_id*grid_dim + x_id, grid_dim);
						//mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_lower(&Tinf, sources[source_id][image_id].dr, &Tsource, grid_gradient_x, grid_gradient_y, (y_id + y_pos_loc)*grid_dim + x_id, grid_dim);
						//@@if (world_rank == 1)
						if (mychi_inside(&sources[source_id][image_id].center, &Tsource) == 1)
						{
							//printf("	-> %d found: y_id = %d, x_id = %d : pixel %f %f -> %f %f\n", world_rank, y_pos_loc + y_id, x_id, Tinf.a.x, Tinf.a.y, Tsource.a.x, Tsource.a.y);
							thread_found_image = 1;
							Timage = Tinf;
						}
					}
#if 1
					if (thread_found_image)
					{
#pragma omp critical
						{
							image_pos[source_id][image_id][loc_images_found] = mychi_barycenter(&Timage); // get the barycenter of the triangle
							locimagesfound[source_id][image_id]++;
							loc_images_found++;
							numimg ++;
						}
					}
#endif
				}
			}
//		if(locimagesfound[source_id][image_id] != 0)
//		{
//			numimg += locimagesfound[source_id][image_id];
			//printf("	-> %d: %d %d: number of images = %d, total = %d\n", world_rank, source_id, image_id, locimagesfound[source_id][image_id], totimg);
//			fflush(stdout);
//		}
		MPI_Barrier(MPI_COMM_WORLD);
#if 1
		}
		index += nimages_strongLensing[source_id];
	}
	//printf("%d: %d images found out of %d total images\n", world_rank, numimg, images_total);
	//MPI_Barrier(MPI_COMM_WORLD);
	//
	double comm_time = -myseconds();
	int          numimagesfound    [runmode->nsets][runmode->nimagestot];
	memset(&numimagesfound, 0, runmode->nsets*runmode->nimagestot*sizeof(int));
	struct point imagesposition    [runmode->nsets][runmode->nimagestot][MAXIMPERSOURCE];
	memset(&imagesposition, 0, runmode->nsets*runmode->nimagestot*MAXIMPERSOURCE*sizeof(point));
	//
	int          numimagesfound_tmp[runmode->nsets][runmode->nimagestot];
	struct point imagesposition_tmp[runmode->nsets][runmode->nimagestot][MAXIMPERSOURCE];
	//
	memset(numimagesfound_tmp, 0, runmode->nsets*runmode->nimagestot*sizeof(int));
	memset(imagesposition_tmp, 0, runmode->nsets*runmode->nimagestot*sizeof(int));
	/*
	//if (verbose)
	{	
		int image_sum = 0;
		for (int ii = 0; ii < runmode->nsets; ++ii)
			for (int jj = 0; jj < runmode->nimagestot; ++jj)
			{
				image_sum += locimagesfound[ii][jj];	
			}
		printf("%d: num images found = %d\n", world_rank, image_sum); 
	}
	MPI_Barrier(MPI_COMM_WORLD);
	*/
#ifdef __WITH_MPI
	int total = 0;
	//MPI_Reduce(&locimagesfound, &imagesfound, runmode->nsets*runmode->nimagestot, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Reduce(&locimagesfound, &imagesfound, runmode->nsets*runmode->nimagestot, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	//
	if (!verbose)
	{
		MPI_CHECK(MPI_Send( &numimg        , 1                                                , MPI_INT   , 0, 666 + world_rank, MPI_COMM_WORLD ));
		if (numimg != 0)
		{
			MPI_CHECK(MPI_Send( &locimagesfound, runmode->nsets*runmode->nimagestot               , MPI_INT   , 0, 666 + world_rank, MPI_COMM_WORLD ));
			//MPI_CHECK(MPI_Send( &image_pos,      runmode->nsets*runmode->nimagestot*MAXIMPERSOURCE, MPI_points, 0, 667 + world_rank, MPI_COMM_WORLD ));
			MPI_CHECK(MPI_Send( &image_pos,      runmode->nsets*runmode->nimagestot*MAXIMPERSOURCE*2, MPI_DOUBLE, 0, 667 + world_rank, MPI_COMM_WORLD ));
		}
	}
	//
	if (verbose)
	{
		int image_sum = 0;
		//
		for (int ipe = 0; ipe < world_size; ++ipe)
		{
			MPI_Status status;
			//
			if (ipe == 0)
			{
				memcpy(&numimagesfound_tmp, &locimagesfound, runmode->nsets*runmode->nimagestot*sizeof(int));
				memcpy(&imagesposition_tmp, &image_pos,      runmode->nsets*runmode->nimagestot*MAXIMPERSOURCE*sizeof(point)); 
			}
			else
			{
				MPI_CHECK(MPI_Recv(&numimg            , 1                                                  , MPI_INT   , ipe, 666 + ipe, MPI_COMM_WORLD, &status));	
				if (numimg != 0)
				{
					MPI_CHECK(MPI_Recv(&numimagesfound_tmp, runmode->nsets*runmode->nimagestot                 , MPI_INT   , ipe, 666 + ipe, MPI_COMM_WORLD, &status));
					//MPI_CHECK(MPI_Recv(&imagesposition_tmp, runmode->nsets*runmode->nimagestot*MAXIMPERSOURCE, MPI_points, ipe, 667 + ipe, MPI_COMM_WORLD, &status));
					MPI_CHECK(MPI_Recv(&imagesposition_tmp, runmode->nsets*runmode->nimagestot*MAXIMPERSOURCE*2, MPI_DOUBLE, ipe, 667 + ipe, MPI_COMM_WORLD, &status));
				}
			}
			//
			//MPI_Reduce(&imagesfound_tmp, &total, ipe, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			//
			if (numimg != 0)
			for (int jj = 0; jj < runmode->nimagestot; ++jj)
			{
				for (int ii = 0; ii < runmode->nsets; ++ii)
				{
					//int img_len = numimagesfound[ii][jj];
					int img_len = numimagesfound_tmp[ii][jj];
					//printf("%d: %d %d, img_len = %d\n", ipe, ii, jj, img_len);
					image_sum  += img_len;
					if (img_len != 0)
					{
						//int loc_length = numimagesfound[ii][jj];
						int loc_length = numimagesfound[ii][jj];
						//printf("%d: %d %d, inserting %d images at position %d\n", ipe, ii, jj, img_len, loc_length);
						memcpy(&imagesposition[ii][jj][loc_length], &imagesposition_tmp[ii][jj], img_len*sizeof(point));
						numimagesfound[ii][jj] += img_len;
						numimg += img_len;

						/*
						if (verbose)
						{
							//for( int  source_id = 0; source_id < runmode->nsets; source_id ++)
							{
								// number of images in the image plane for the specific image (1,3,5...)
								//unsigned short int nimages = nimages_strongLensing[source_id];
								//____________________________ image (constrains) loop ________________________________
								//for(unsigned short int image_id = 0; image_id < nimages; image_id++)
								{
									int img_len = numimagesfound[ii][jj];
									for (int ij = 0; ij < img_len; ++ij)
										printf("	-> images = %f %f\n", imagesposition[ii][jj][ij].x, imagesposition[ii][jj][ij].y);
										//printf("ipe %d: source = %d image = %d: putting %d images to position %d number of images = %d, %f %f, total = %d\n", ipe, ii, jj, numimagesfound[ii][jj], imagesposition[ii][jj][ij].x, imagesposition[ii][jj][ij].y, totimg);
								}
							}
						}
						*/


						//printf("%d: %d %d, img_len = %d, loc_length = %d\n", ipe, ii, jj, img_len, numimagesfound[ii][jj]);
						//for (int ij = 0; ij < img_len; ++ij)
						//	printf("	%f %f\n", imagesposition[ii][jj][loc_length + ij].x, imagesposition[ii][jj][loc_length + ij].y);
					}
					/*
					   for (int ij = 0; ii < img_len; ++ij)
					   {
					   int loc_length = numimagesfound_tmp[ii][jj];
					   printf("	%d: point = %f %f\n", ij, imagesposition[ii][jj][loc_length + ij].x, imagesposition[ii][jj][loc_length + ij].y);
					   }
					   */
				}
			}
		}
		//printf("%d: num images found = %d, size of point = %d\n", ipe, image_sum, sizeof(point));
	}
	//
	//
	/*
	if (0*verbose)
		for( int  source_id = 0; source_id < runmode->nsets; source_id ++)
		{
			// number of images in the image plane for the specific image (1,3,5...)
			unsigned short int nimages = nimages_strongLensing[source_id];
			//____________________________ image (constrains) loop ________________________________
			for(unsigned short int image_id = 0; image_id < nimages; image_id++)
			{
				int img_len = numimagesfound[source_id][image_id];
				for (int ij = 0; ij < img_len; ++ij)
					printf("***** %d: %d %d: number of images = %d, %f %f, total = %d\n", world_rank, source_id, image_id, numimagesfound[source_id][image_id], imagesposition[source_id][image_id][ij].x, imagesposition[source_id][image_id][ij].y, totimg);
			}
		}
	*/
	//
	//
	//MPI_Barrier(MPI_COMM_WORLD);
	comm_time += myseconds();
	//
	//memcpy(imagesfound, locimagesfound, runmode->nsets*runmode->nimagestot*sizeof(int));
	//	
	/*
	   if (verbose) printf("--> total images found = %ld\n", images_total);
	   int total_images_found[runmode->nsets][runmode->nimagestot];		
	//if (verbose)
	memset(&total_images_found, 0, runmode->nsets*runmode->nimagestot*sizeof(int)); 
	//MPI_Reduce(&locimagesfound[0], &images_found, runmode->nsets*runmode->nimagestot*sizeof(int), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
	for (int ii = 0; ii < runmode->nsets; ++ii)
	for (int jj = 0; jj < runmode->nimagestot; ++jj)
	total += total_images_found[ii][jj];	
	//
	//printf("--> total images found = %d\n", total, loc_images_found);
	printf("--> total images found = %d\n", total);
	images_total = total;
	*/
#endif
	//	
	loop_time += myseconds();
	//
	// image extraction to compute 
	//
	chi2_time = -myseconds();
	index     =  0;
	//
	if (verbose)
		for( int  source_id = 0; source_id < runmode->nsets; source_id ++)
		{
			unsigned short int nimages = nimages_strongLensing[source_id];
			//
			//______________________Initialisation______________________
			//
			for (unsigned short int i = 0; i < nimages; ++i)
				for (unsigned short int j = 0; j < nimages; ++j)
					nimagesfound[source_id][i][j] = 0;
			//
			for( unsigned short int image_id = 0; image_id < nimages; image_id++)
			{
#endif
				//
				double image_dist[MAXIMPERSOURCE];
				//
				//printf("	Image %d, number of sources found %d\n", image_id, loc_images_found);
				//
				struct point image_position;
				int	     image_index;
				//
				for (int ii = 0; ii < /*loc*/numimagesfound[source_id][image_id]; ++ii)
				{	
					//
					//
					int image_index = 0;
					//
					//image_position = imagesposition_tmp[source_id][image_id][ii]; 
					image_position = imagesposition[source_id][image_id][ii]; 
					//image_position = image_pos[source_id][image_id][ii]; 
					//image_position = image_pos[ii]; 
					//printf("	* source %d, image %d, %d = %f %f\n", source_id, image_id, ii, image_position.x, image_position.y);
					//printf("	* source %d, image %d, %d = %f %f\n", source_id, image_id, ii, imagesposition/*_tmp*/[source_id][image_id][ii].x, imagesposition/*_tmp*/[source_id][image_id][ii].y); 
					image_dist[0] = mychi_dist(image_position, images[index + 0].center);  // get the distance to the real image
					//printf("        *** image %d = %f %f, distance = %f\n", index, images[index + 0].center.x, images[index + 0].center.y, image_dist[0]);
					for(int i = 1; i < nimages_strongLensing[source_id]; i++)
					{  // get the distance to each real image and keep the index of the closest real image

						image_dist[i] = mychi_dist(image_position, images[index + i].center);
						//printf("	*** image %d = %f %f, distance = %f\n", index + i, images[index + i].center.x, images[index + i].center.y, image_dist[i]);
						if (image_dist[i] < image_dist[image_index])
						{
							image_index = i;
						}
					}
					//printf("	-> %f %f %d\n", image_position.x, image_position.y, image_index);
					//
					// we should exit loops here
					//
					// p1_time += myseconds();
					//
					//if (thread_found_image == 1)
					//printf("%d %d %d: found image\n", x_id, y_id, image_id); 
					int skip_image = 0;
					// Sometimes due to the numerical errors at the centerpoint, 
					// for SIE potentials an additional image will appear at the center of the Potential.
					// This is due to the fact that it is not possible to simulate an infinity value 
					// at the center correctly, Check that sis correspond to Nlens[0]
					for (int iterator = 0; iterator < runmode->Nlens[0]; ++iterator)
					{
						if ( fabs(image_position.x - lens[0].position_x[iterator]) <= dx/2. and fabs(image_position.y  - lens[0].position_y[iterator]) <= dx/2.)
						{
							skip_image = 1;
							printf("WARNING: You are using SIE potentials. An image to close to one of the potential centers has been classified as numerical error and removed \n");
						}
					}
					//printf("%d %d %d %d\n", x_id, y_id,thread_found_image, skip_image);
					if (!skip_image)
					{
						//#pragma omp atomic
						images_found++;
						struct point temp;
						//printf("                        source %d, image %d, index %d, Images found: %d\n", source_id, image_id, image_index, nimagesfound[source_id][image_id][image_index]);
						//checking whether a closest image has already been found
						if (nimagesfound[source_id][image_id][image_index] == 0)
						{ // if no image found up to now

							//image position is allocated to theoretical image
							//#pragma omp critical
							tim[image_id][image_index] = image_position;  
							//#pragma omp atomic
							nimagesfound[source_id][image_id][image_index]++;
						}
						else if (nimagesfound[source_id][image_id][image_index] > 0)
						{ // if we have already found an image
							// If the new image we found is closer than the previous image
							//printf("	tim2: %f %f\n", image_dist[image_index], mychi_dist(images[index + image_index].center, tim[image_id][image_index]));
							if (image_dist[image_index] < mychi_dist(images[index + image_index].center, tim[image_id][image_index]))
							{
								temp = tim[image_id][image_index]; // we store the position of the old image in temp
								//#pragma omp critical
								tim[image_id][image_index] = image_position; // we link the observed image with the image we just found
								//printf("tim2 %d %d = %f %f\n", image_id, image_index, image_position.x, image_position.y);
							}
							else
							{
								temp = image_position; // we store the position of the image we just found in temp
							}
							// initialising second_closest_id to the highest value
							// Loop over all images in the set except the closest one
							// and initialize to the furthest away image
							int second_closest_id = 0;
							for (int i = 1; i < nimages_strongLensing[source_id] && (i != image_index); i++)
							{
								if(image_dist[i] > image_dist[second_closest_id]) second_closest_id=i;
							}
							///////////////////////////////////////////////////////////////
							// Loop over all images in the set that are not yet allocated to a theoretical image
							// and allocate the closest one
							// we search for an observed image not already linked (nimagesfound=0)
							for(int i = 0; i < nimages_strongLensing[source_id] && nimagesfound[source_id][image_id][i] == 0; i++) 
							{
								if(image_dist[i] < image_dist[second_closest_id])
								{
									second_closest_id = i;
									// im_index value changes only if we found a not linked yet image
									image_index = i; 
									//printf("tim3 %d %d = %f %f\n", image_id, image_index, temp.x, temp.y);
									//#pragma omp critical
									tim[image_id][image_index] = temp; // if we found an observed and not already linked image, we allocate the theoretical image temp
								}
							}
							// increasing the total number of images found (If we find more than 1 theoretical image linked to 1 real image, 
							// these theoretical
							//#pragma omp atomic
							nimagesfound[source_id][image_id][image_index]++; 
							// images are included in this number)
						}
					}
					//#pragma omp atomic
					//loc_images_found++;
					//thread_found_image  = 0; // for next iteration
				}
				//
			}
			//#pragma omp barrier
			//____________________________ end of image loop
			//
			//____________________________ computing the local chi square
			//
			double chiimage;
			//
			int _nimages = nimages_strongLensing[source_id];
			//
			for( int iter = 0; iter < _nimages*_nimages; iter++)
			{
				int i=iter/nimages_strongLensing[source_id];
				int j=iter%nimages_strongLensing[source_id];

				//printf("nimagesfound %d %d = %d\n", i, j, nimagesfound[i][j]);
				if(i != j)
				{
					// In the current method, we get the source in the source plane by ray tracing image in nimagesfound[i][i]. If we ray trace back,
					// we arrive again at the same position and thus the chi2 from it is 0. Thus we do not calculate the chi2 (-> if i!=j)
					if(nimagesfound[source_id][i][j] > 0)
					{
						double pow1 = images[index + j].center.x - tim[i][j].x;
						double pow2 = images[index + j].center.y - tim[i][j].y;
						//
						//chiimage = pow(images[index + j].center.x - tim[i][j].x, 2) + pow(images[index + j].center.y - tim[i][j].y, 2);  // compute the chi2
						chiimage = pow1*pow1 + pow2*pow2;  // compute the chi2
						//printf("%d %d = %.15f\n", i, j, chiimage);
						*chi    += chiimage;
					}
					else 
						if(nimagesfound[source_id][i][j] == 0)
						{
							// If we do not find a correpsonding image, we add a big value to the chi2 to disfavor the model
							*chi += 100.*nimages_strongLensing[source_id];
						}
				}
			}
			//
			//____________________________ end of computing the local chi square
			//
			//printf("%d: chi = %.15f\n", source_id, *chi);
			/*
			   for (int i=0; i < nimages_strongLensing[source_id]; ++i){
			   for (int j=0; j < nimages_strongLensing[source_id]; ++j){
			   printf(" %d",nimagesfound[i][j]);
			   }
			   printf("\n");
			   }*/

			//Incrementing Index: Images already treated by previous source_id
			index += nimages_strongLensing[source_id];
		}
	MPI_Barrier(MPI_COMM_WORLD);
	//
	chi2_time += myseconds();
	time      += myseconds();
	//
	if (verbose)
	{
		//
		int nthreads = 1;
		//
#pragma omp parallel
		nthreads = omp_get_num_threads();
		//
		printf("	overall time  = %f s. using %d threads\n", time, nthreads);
		printf("		- image  time = %f s.\n", image_time);
		printf("		- loop   time = %f s.\n", loop_time);
		printf("		- comm   time = %f s.\n", comm_time);
		printf("		- chi2   time = %f s.\n", chi2_time);
		//
		printf("	images found: %d out of %ld\n", images_found, images_total);
	}
	//
	free(grid_gradient_x);
	free(grid_gradient_y);
}


/** @brief Tranform a point from image to source plane. Result stored in sourcepoint argument
 *
 * Tranform a point from image to source plane using lensequation
 *
 * @param image_point	 image position
 * @param dlsds		 dls/ds
 * @param nhalos	 number of halos
 * @param potential_param	 gravitational potential information
 * @param source_point	 address where source information will be stored
 *
 *
 */
void mychi_transformImageToSourcePlane(const runmode_param *runmode, const struct point *image_point, double dlsds, const struct Potential *lens, struct point *source_point)
{   // dlsds is the distance between lens and source divided by the distance observer-source
	struct point Grad;  // gradient

	Grad = module_potentialDerivatives_totalGradient(runmode->nhalos, image_point, lens);
	//Grad = module_potentialDerivatives_totalGradient_SOA(image_point, lens, runmode->Nlens);

	source_point->x = image_point->x - dlsds*Grad.x;
	source_point->y = image_point->y - dlsds*Grad.y;
	//printf("dlsds %f", dlsds);
}


void mychi_transformImageToSourcePlane_SOA(const int Nlens, const struct point *image_point, double dlsds, const struct Potential_SOA *lens, struct point *source_point)
{   // dlsds is the distance between lens and source divided by the distance observer-source
	struct point Grad;  // gradient
	Grad = module_potentialDerivatives_totalGradient_SOA(image_point, lens, Nlens);
	//
	source_point->x = image_point->x - dlsds*Grad.x;
	source_point->y = image_point->y - dlsds*Grad.y;
	//printf("dlsds %f", dlsds);
}
//
//
//
	inline
void mychi_transformImageToSourcePlane_SOA_Packed( const struct point *image_point, double dlsds, struct point *source_point, double *grad_x, double * grad_y, int grad_id)
{

	source_point->x = image_point->x - dlsds*grad_x[grad_id];
	source_point->y = image_point->y - dlsds*grad_y[grad_id];
	int world_rank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	//if (world_rank == 1) 
	//printf("	%d: %f %f = %f %f - dlsds = %f grad id = %d grad = (%f %f)\n", world_rank, source_point->x, source_point->y, image_point->x, image_point->y, dlsds, grad_id, grad_x[grad_id], grad_y[grad_id]);
	//printf("dlsds %f", dlsds);
}


/** @brief Tranform a triangle from image to source plane. Result stored in S triangle argument
 *
 * Return a triplet of points in the source plane corresponding to the triplet
 * of images. dlsds is the lens efficiency at the source redshift.
 * I is the triangle in the image plane (input), S is the same triangle in the source plane (output)
 *
 * @param I	 triangle in image plane
 * @param dlsds	 dls/ds
 * @param nhalos	 number of halos
 * @param potential_param	 gravitational potential information
 * @param S	 address where triangle source information will be stored
 *
 *
 */

	inline
void mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_upper( struct triplet *I, double dlsds, struct triplet *S, double *grad_x, double * grad_y, int grad_id, int nbgridcell)
{
	mychi_transformImageToSourcePlane_SOA_Packed( &I->a, dlsds,   &S->a, grad_x, grad_y, grad_id             );
	mychi_transformImageToSourcePlane_SOA_Packed( &I->b, dlsds,   &S->b, grad_x, grad_y, grad_id + nbgridcell);
	mychi_transformImageToSourcePlane_SOA_Packed( &I->c, dlsds,   &S->c, grad_x, grad_y, grad_id +          1);
}

	inline
void mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_lower( struct triplet *I, double dlsds, struct triplet *S, double *grad_x, double * grad_y, int grad_id, int nbgridcell)
{
	mychi_transformImageToSourcePlane_SOA_Packed( &I->a, dlsds,   &S->a, grad_x, grad_y, grad_id + nbgridcell + 1);
	mychi_transformImageToSourcePlane_SOA_Packed( &I->b, dlsds,   &S->b, grad_x, grad_y, grad_id +              1);
	mychi_transformImageToSourcePlane_SOA_Packed( &I->c, dlsds,   &S->c, grad_x, grad_y, grad_id + nbgridcell    );
}



/** @brief Return the scalar triple product (a*b).c of the 3 vectors A[x,y,1], B[x,y,1], C[x,y,1].
 * If 2 of the 3 vectors are equal, colinear or form an orthogonal basis,
 * the triple product is 0.
 * This is also the determinant of the matrix
 *   | Ax  Bx  Cx |
 *   | Ay  By  Cy |
 *   |  1   1   1 |
 */
//inline
double mychi_determinant(const struct point *A,
		const struct point *B,
		const struct point *C)
{
	return( B->x*C->y - B->y*C->x +
			A->x*B->y - A->y*B->x +
			A->y*C->x - A->x*C->y );
}


/** @brief Return 1 if P is inside the triangle T, 0 otherwise.
 * Return 1 if P is inside the triangle T, 0 otherwise.
 * @param P  a point
 * @param T  a triplet of points.
 *
 *
 */
	inline
int mychi_inside(const struct point *P, struct triplet *T)
{
	double  s, s1, s2, d;

	d  = mychi_determinant(&T->a, &T->b, &T->c);
	s  = mychi_determinant(&T->a, &T->b, P)*d;
	if (s < 0.) return 0;
	s1 = mychi_determinant(&T->b, &T->c, P)*d;
	if (s1 < 0.) return 0;
	s2 = mychi_determinant(&T->c, &T->a, P)*d;
	if (s2 < 0.) return 0;
	return 1;

	return((s > 0.) && (s1 > 0.) && (s2 > 0.));  // If all determinants are positive,
	// the point must be inside the triangle
}


/*
   int
   mychi_inside2(const struct point *A, const struct point *B, const struct point *C)
   {

// Compute vectors        
v0 = C - A;
v1 = B - A;
v2 = P - A;

// Compute dot products
dot00 = dot(v0, v0);
dot01 = dot(v0, v1);
dot02 = dot(v0, v2);
dot11 = dot(v1, v1);
dot12 = dot(v1, v2);

// Compute barycentric coordinates
invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
u = (dot11 * dot02 - dot01 * dot12) * invDenom;
v = (dot00 * dot12 - dot01 * dot02) * invDenom;

// Check if point is in triangle
return (u >= 0) && (v >= 0) && (u + v < 1);
}
*/

/** @brief Return 1 if P is inside the triangle T or on its border, 0 otherwise.
 *
 * Return 1 if P is inside the triangle T or on its border, 0 otherwise.
 * @param  P  a point
 * @param  T  a triplet of points.
 *
 *
 */
	inline
int mychi_insideborder(const struct point *P, struct triplet *T)
{
	double  s, s1, s2, d;

	d  = mychi_determinant(&T->a, &T->b, &T->c);
	s  = mychi_determinant(&T->a, &T->b, P)*d;
	if (s < 0.) return 0;
	s1 = mychi_determinant(&T->b, &T->c, P)*d;
	if (s1 < 0.) return 0;
	s2 = mychi_determinant(&T->c, &T->a, P)*d;
	if (s2 < 0.) return 0;
	return 1;
	return((s >= 0.) && (s1 >= 0.) && (s2 >= 0.));  // If all determinants are positive or 0,
	// the point must be inside the triangle or on its border

}

/** @brief Barycentre of a triplet/triangle
 *
 * A is a structure triplet that contains 3 structures point a,b and c
 * Return value B is a point
 *
 *
 */
struct  point   mychi_barycenter(struct triplet *A)
{
	struct  point   B;

	B.x = (A->a.x + A->b.x + A->c.x) / 3.;
	B.y = (A->a.y + A->b.y + A->c.y) / 3.;
	return(B);
}

/** @brief Euclidean distance between 2 points
 *
 * Euclidean distance between 2 points
 *
 */
double mychi_dist(struct point A, struct point B)
{
	double  x, y;
	x = A.x - B.x;
	y = A.y - B.y;
	return(sqrt(x*x + y*y));
}

