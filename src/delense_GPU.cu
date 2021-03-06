// This file is part of lenstoolHPC
// authors: gilles.fourestey@epfl.ch

#include "delense_CPU_utils.hpp"
#include "delense_GPU_utils.cuh"
#include "structure_hpc.hpp"
#include "cudafunctions.cuh"
#ifdef __WITH_MPI
#include <mpi.h>
#endif

extern double myseconds();

#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 8


#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))


inline
void
moveToGPU(void** d_mem, void* h_mem, int size)
{
#ifdef __WITH_UM
	cudasafe(cudaMallocManaged((void**) d_mem, size), "memory allocation (UM)");
	cudasafe(cudaMemcpy(*d_mem, h_mem, size, cudaMemcpyHostToDevice), "H2D memory copy (UM)");
#else
	cudasafe(cudaMalloc((void**) d_mem, size), "memory allocation");
	cudasafe(cudaMemcpy(*d_mem, h_mem, size, cudaMemcpyHostToDevice), "H2D memory copy");
#endif
}

inline
void
moveFromGPU(void* h_mem, void* d_mem, int size)
{
#ifdef __WITH_UM
	//cudasafe(cudaMemcpy(h_mem, d_mem, size, cudaMemcpyDeviceToHost), "D2H memory copy (UM)");
	memcpy(h_mem, d_mem, size);
#else
	cudasafe(cudaMemcpy(h_mem, d_mem, size, cudaMemcpyHostToDevice), "D2H memory copy");
#endif
}
//
//
//
__global__
void
delense_GPU(struct point* image_pos, const runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const struct galaxy* sources, const int y_pos_loc, const int y_bound, const int source_id, double* grid_gradient_x, double* grid_gradient_y)
{
#define INDEX2D_BAR(y, x)         (MAXIMPERSOURCE*y + x)
        //
	int nbgridcells_x = runmode->nbgridcells;
	int nbgridcells_y = runmode->nbgridcells;
	//
        int col = blockIdx.x*blockDim.x + threadIdx.x;
        if (col > nbgridcells_x - 2) return;
	//
	int num_blocks_y  = blockDim.y;
	int row_y	  = threadIdx.y;
        int block_size_y  = y_bound/num_blocks_y;
        //
        const double dx   = (frame->xmax - frame->xmin)/(runmode->nbgridcells - 1);
        const double dy   = (frame->ymax - frame->ymin)/(runmode->nbgridcells - 1);
	//
	int img_pos       = 0;
	//
	int row_s  =     (row_y + 0)*block_size_y;
        int row_e  = MIN((row_y + 1)*block_size_y, nbgridcells_y);
	//
        //
	//for (int row = y_pos_loc; row < y_pos_loc + y_bound - 1; ++row)
	int yp = y_pos_loc + row_s; 
	//
	//printf("block id = %d, from %d to %d\n", row_y, yp, yp + block_size_y);
	//
	for (int row = 0; row < block_size_y - 0; ++row)
	{
		//
		int y_id = row;
		int x_id = col;
		//
		// int images_total;
		//
		double x_pos = frame->xmin + ( x_id +  0 )*dx;
		double y_pos = frame->ymin + ( y_id + yp )*dy;
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
		int g_loc = (y_id + row_s)*nbgridcells_x + x_id;
		//
		mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_upper_GPU(&Tsup, sources[source_id].dr, &Tsource, grid_gradient_x, grid_gradient_y, g_loc, nbgridcells_x);
		//
		//int thread_found_image = 0;
		//
		if (mychi_insideborder_GPU(&sources[source_id].center, &Tsource) == 1)
		{
//			printf("col = %d, %d, row = %d, %d, img_pos = %d\n", col, col*num_blocks_y*MAXIMPERSOURCE, row_y, row_y*MAXIMPERSOURCE, img_pos);
			//printf("found source = %d, row = %d, index 1 = %d, index 2 = %d, img_pos = %d\n", source_id, row + yp, col*num_blocks_y*MAXIMPERSOURCE, row_y*MAXIMPERSOURCE, img_pos);
			image_pos[col*num_blocks_y*MAXIMPERSOURCE + row_y*MAXIMPERSOURCE + img_pos + 0] = mychi_barycenter_GPU(&Tsup);
			img_pos++;
		}
		else
		{
			mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_lower_GPU(&Tinf, sources[source_id].dr, &Tsource, grid_gradient_x, grid_gradient_y, g_loc, nbgridcells_x);
			if (mychi_inside_GPU(&sources[source_id].center, &Tsource) == 1)
			{
//				printf("col = %d, %d, row = %d, %d, img_pos = %d\n", col, col*num_blocks_y*MAXIMPERSOURCE, row_y, row_y*MAXIMPERSOURCE, img_pos);
				//printf("found source = %d, row = %d, index 1 = %d, index 2 = %d, img_pos = %d\n", source_id, row, col*num_blocks_y*MAXIMPERSOURCE, row_y*MAXIMPERSOURCE, img_pos);
				image_pos[col*num_blocks_y*MAXIMPERSOURCE + row_y*MAXIMPERSOURCE + img_pos + 0] = mychi_barycenter_GPU(&Tinf);
				img_pos++;
			}
		}
	}
}
//
//
//
__global__
void
delense_GPU_old(struct point* image_pos, const runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const struct galaxy* sources, const int y_pos_loc, const int y_bound, const int source_id, double* grid_gradient_x, double* grid_gradient_y)
{
#define INDEX2D_BAR(y, x)         (MAXIMPERSOURCE*y + x)
        //
        int nbgridcells_x = runmode->nbgridcells;
        int nbgridcells_y = runmode->nbgridcells;
        int row = blockIdx.x*blockDim.x + threadIdx.x;
        if (row > y_bound) return;
        int num_blocks_y = gridDim.y;
        int block_size_y = nbgridcells_y/num_blocks_y;
        int x_pos_loc    = threadIdx.y*block_size_y;
        //
        //int col = blockIdx.y*blockDim.y + threadIdx.y;
        //if (col > runmode->nbgridcells) return;
        //
        const double dx = (frame->xmax - frame->xmin)/(runmode->nbgridcells - 1);
        const double dy = (frame->ymax - frame->ymin)/(runmode->nbgridcells - 1);
        //
        int img_pos = 0;
        //
        for (int col = 0; col < runmode->nbgridcells; ++col)
        {
                //
                int y_id = row;
                int x_id = col;
                //
                //int images_total;
                //
                double x_pos = frame->xmin + (x_pos_loc + x_id)*dx;
                double y_pos = frame->ymin + (y_pos_loc + y_id)*dy;
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
                mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_upper_GPU(&Tsup, sources[source_id].dr, &Tsource, grid_gradient_x, grid_gradient_y, y_id*runmode->nbgridcells + x_id, runmode->nbgridcells);
                //
                int thread_found_image = 0;
                //
                if (mychi_insideborder_GPU(&sources[source_id].center, &Tsource) == 1)
                {
                        image_pos[row*MAXIMPERSOURCE + img_pos + 0] = mychi_barycenter_GPU(&Tsup);
                        //printf("%d %d: T = %f %f\n", row, img_pos, image_pos[row*MAXIMPERSOURCE + img_pos + 0].x, image_pos[row*MAXIMPERSOURCE + img_pos + 0].y);
                        img_pos++;
                }
                else
                {
                        mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_lower_GPU(&Tinf, sources[source_id].dr, &Tsource, grid_gradient_x, grid_gradient_y, y_id*runmode->nbgridcells + x_id, runmode->nbgridcells);
                        if (mychi_inside_GPU(&sources[source_id].center, &Tsource) == 1)
                        {
                                image_pos[row*MAXIMPERSOURCE + img_pos + 0] = mychi_barycenter_GPU(&Tinf);
                                //printf("%d %d: T = %f %f\n", row, img_pos, image_pos[row*MAXIMPERSOURCE + img_pos + 0].x, image_pos[row*MAXIMPERSOURCE + img_pos + 0].y);
                                img_pos++;
                        }
                }
        }
}
//
//
//
#if 1
void delense_barycenter_GPU(struct point* image_pos, int* locimagesfound, int* numimg, const runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, const struct galaxy* sources, double* grid_gradient_x, double* grid_gradient_y)
{
#define INDEX2D_BAR(y, x)         (MAXIMPERSOURCE*y + x)
	//const unsigned int nimagestot  = runmode->nimagestot;
	const unsigned int nsets         = runmode->nsets;
	const unsigned int nbgridcells   = runmode->nbgridcells;
	const unsigned int nbgridcells_x = runmode->nbgridcells;
	const unsigned int nbgridcells_y = runmode->nbgridcells;
	//
	int world_size = 1;
	int world_rank = 0;
	//
#ifdef __WITH_MPI
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif
	//
	unsigned int verbose = (world_rank == 0);
	//
	int numimg_loc = 0;
	//
	for (int ii = 0; ii < nsets; ++ii)
		numimg_loc += nimages_strongLensing[ii];
	//
	//const int nbgridcells   = runmode->nbgridcells;
	const int grid_size     = nbgridcells;
	const int loc_grid_size = nbgridcells/world_size;
	//
	double y_len      = fabs(frame->ymax - frame->ymin);
	int    y_len_loc  = nbgridcells/world_size;
	int    y_pos_loc  = (int) world_rank*y_len_loc;
	int    y_bound    = y_len_loc;
	//
	if ((world_rank + 1) != world_size) y_bound += 1;
	//
	printf("%d: numimg_loc = %d, y_len_loc = %d, y_pos_loc = %d, y_bound = %d\n", world_rank, numimg_loc, y_len_loc, y_pos_loc, y_bound);
	//
	const double dx   = (frame->xmax - frame->xmin)/(nbgridcells - 1);
	const double dy   = (frame->ymax - frame->ymin)/(nbgridcells - 1);
	//
	int images_total  = 0;
	int index         = 0;
	//
	int block_size_x  = BLOCK_SIZE_X;
	int block_size_y  = BLOCK_SIZE_Y;
	//
	int GRID_SIZE_X   = (nbgridcells_y + block_size_y - 1)/block_size_y; // number of blocks
	int GRID_SIZE_Y   = 1; 
	//
	dim3 threads(block_size_x, block_size_y);
	dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);
	//
	double alloc_time = -myseconds();
	struct galaxy* sources_gpu;
	moveToGPU((void**) &sources_gpu, (void*) sources, sizeof(galaxy)*nsets); 
	//
	struct grid_param* frame_gpu;
	moveToGPU((void**) &frame_gpu  , (void*) frame  , sizeof(grid_param));
	//
	struct runmode_param* runmode_gpu;
	moveToGPU((void**) &runmode_gpu, (void*) runmode, sizeof(runmode_param));
	//
	struct point* image_pos_gpu;
	cudaMallocManaged((void**) &image_pos_gpu, grid_size*block_size_y*MAXIMPERSOURCE*sizeof(struct point));
	//
	//printf("image_pos size = %d\n", grid_size*block_size_y*MAXIMPERSOURCE);
	//
#if 1
	//
	memset(locimagesfound, 0, nsets*sizeof(int)); 
	//
	cudaDeviceSynchronize();
	cudasafe(cudaGetLastError(), "before image_pos_gpu allocation ");
	alloc_time += myseconds();
	//
#endif
	//
	double time_gpu, time_cpu;
	//
	time_gpu = 0.;
	time_cpu = 0.;
	//
	for( int  source_id = 0; source_id < nsets; source_id ++)
	{
		int loc_images_found = 0;
		//
		int x_pos_loc = 0;
		cudasafe(cudaGetLastError(), "before delense_GPU");
		time_gpu -= myseconds();
		//
		delense_GPU<<<grid, threads>>> (image_pos_gpu, runmode_gpu, lens, frame_gpu, sources_gpu, y_pos_loc, y_bound, source_id, grid_gradient_x, grid_gradient_y);
		//
		cudaDeviceSynchronize();
		cudasafe(cudaGetLastError(), "after  delense_GPU");
		time_gpu += myseconds();
		fflush(stdout);
#if 1
		time_cpu -= myseconds();
		int numimg_cpu = 0;
#pragma omp parallel for reduction(+ : numimg_cpu)
		for (int ii = 0; ii < block_size_y*nbgridcells_x*MAXIMPERSOURCE; ++ii)
		{
				struct point T = image_pos_gpu[ii];
				if((T.x != 0.) || (T.y != 0.))
				{
#pragma omp critical
					{
						image_pos[source_id*MAXIMPERSOURCE + locimagesfound[source_id]] = T; 
						locimagesfound[source_id]++;
						numimg_cpu++;
					}
					image_pos_gpu[ii].x = image_pos_gpu[ii].y = 0.;
				}
		}
		//memset(image_pos_gpu, 0, nbgridcells*MAXIMPERSOURCE*sizeof(struct point));
		//memset(image_pos_gpu, 0, loc_grid_size*MAXIMPERSOURCE*sizeof(struct point));
		*numimg = *numimg + numimg_cpu; 
		time_cpu += myseconds();
		//memset(image_pos_gpu, 0, 2*nbgridcells*nbgridcells*sizeof(struct point));
#endif
	}
	cudaFree(image_pos_gpu);
	printf("num img found = %d, gpu time = %f, cpu time = %f\n", *numimg, time_gpu, time_cpu);
}
#endif
