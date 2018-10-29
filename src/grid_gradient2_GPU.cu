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

@brief: Function for second order derivative computation over a grid

*/


#include <fstream>
#include "grid_gradient2_GPU.cuh"
#include "gradient2_GPU.cuh"
#include <structure_hpc.hpp>

#define BLOCK_SIZE_X 32
#define BLOCK_SIZE_Y 16

//#define ROT

#define _SHARED_MEM

#ifdef _SHARED_MEM
#define SHARED __shared__
#warning "shared memory"
extern __shared__ type_t shared[];
#else
#define SHARED 
#endif

#define Nx 1
#define Ny 0


#define cudasafe 

extern "C" 
{
	type_t myseconds();
}

__global__ void module_potentialDerivatives_totalGradient2_SOA_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d,  const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int nhalos);

////
void
module_potentialDerivatives_totalGradient2_SOA_CPU_GPU(type_t *grid_grad2_a, type_t *grid_grad_b, type_t *grid_grad2_c, type_t *grid_grad_d, const struct grid_param *frame, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos);
//
void gradient2_grid_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart);
//
//void
//module_potentialDerivatives_totalGradient_SOA_CPU_GPU_v2(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_cpu, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos);
//
//
//
void gradient2_grid_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos ,int nbgridcells)
{
	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells - 1);
    type_t dy = (frame->ymax - frame->ymin)/(nbgridcells - 1);
        //
    gradient2_grid_GPU(grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d, frame, lens, nhalos, dx, dy, nbgridcells, nbgridcells, 0, 0);
}
//
//
//
void gradient2_grid_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{

	int nBlocks_gpu = 0;
	// Define the number of threads per block the GPU will use
	cudaDeviceProp properties_gpu;

	cudaGetDeviceProperties(&properties_gpu, 0); // Get properties of 0th GPU in use

/*
	if (properties_gpu.maxThreadsDim[0]<threadsPerBlock)
	{
		fprintf(stderr, "ERROR: The GPU has to support at least %u threads per block.\n", threadsPerBlock);
		exit(-1);
	}
	else
	{
		nBlocks_gpu = properties_gpu.maxGridSize[0] / threadsPerBlock;  // Get the maximum number of blocks with the chosen number of threads
		// per Block that the GPU supports
	}
*/
	grid_param *frame_gpu;
	Potential_SOA *lens_gpu,*lens_kernel;
	int *type_gpu;
	type_t *lens_x_gpu, *lens_y_gpu, *b0_gpu, *angle_gpu, *epot_gpu, *rcore_gpu, *rcut_gpu, *anglecos_gpu, *anglesin_gpu;
	type_t *grid_grad2_a_gpu, *grid_grad2_b_gpu , *grid_grad2_c_gpu, *grid_grad2_d_gpu;

	lens_gpu = (Potential_SOA *) malloc(sizeof(Potential_SOA));
	lens_gpu->type = (int *) malloc(sizeof(int));

	// Allocate variables on the GPU
	cudasafe(cudaMalloc( (void**)&(lens_kernel), sizeof(Potential_SOA)),"Gradient2gpu.cu : Alloc Potential_SOA: " );
	cudasafe(cudaMalloc( (void**)&(type_gpu), nhalos*sizeof(int)),"Gradient2gpu.cu : Alloc type_gpu: " );
	cudasafe(cudaMalloc( (void**)&(lens_x_gpu), nhalos*sizeof(type_t)),"Gradient2gpu.cu : Alloc x_gpu: " );
	cudasafe(cudaMalloc( (void**)&(lens_y_gpu), nhalos*sizeof(type_t)),"Gradient2gpu.cu : Alloc y_gpu: " );
	cudasafe(cudaMalloc( (void**)&(b0_gpu), nhalos*sizeof(type_t)),"Gradient2gpu.cu : Alloc b0_gpu: " );
	cudasafe(cudaMalloc( (void**)&(angle_gpu), nhalos*sizeof(type_t)),"Gradient2gpu.cu : Alloc angle_gpu: " );
	cudasafe(cudaMalloc( (void**)&(epot_gpu), nhalos*sizeof(type_t)),"Gradient2gpu.cu : Alloc epot_gpu: " );
	cudasafe(cudaMalloc( (void**)&(rcore_gpu), nhalos*sizeof(type_t)),"Gradient2gpu.cu : Alloc rcore_gpu: " );
	cudasafe(cudaMalloc( (void**)&(rcut_gpu), nhalos*sizeof(type_t)),"Gradient2gpu.cu : Alloc rcut_gpu: " );
	cudasafe(cudaMalloc( (void**)&(anglecos_gpu), nhalos*sizeof(type_t)),"Gradient2gpu.cu : Alloc anglecos_gpu: " );
	cudasafe(cudaMalloc( (void**)&(anglesin_gpu), nhalos*sizeof(type_t)),"Gradient2gpu.cu : Alloc anglesin_gpu: " );
	cudasafe(cudaMalloc( (void**)&(frame_gpu), sizeof(grid_param)),"Gradient2gpu.cu : Alloc frame_gpu: " );
	cudasafe(cudaMalloc( (void**)&(grid_grad2_a_gpu), (nbgridcells_x) * (nbgridcells_y) *sizeof(type_t)),"Gradient2gpu.cu : Alloc source_a_gpu: " );
	cudasafe(cudaMalloc( (void**)&(grid_grad2_b_gpu), (nbgridcells_x) * (nbgridcells_y) *sizeof(type_t)),"Gradient2gpu.cu : Alloc source_b_gpu: " );
	cudasafe(cudaMalloc( (void**)&(grid_grad2_c_gpu), (nbgridcells_x) * (nbgridcells_y) *sizeof(type_t)),"Gradient2gpu.cu : Alloc source_c_gpu: " );
	cudasafe(cudaMalloc( (void**)&(grid_grad2_d_gpu), (nbgridcells_x) * (nbgridcells_y) *sizeof(type_t)),"Gradient2gpu.cu : Alloc source_d_gpu: " );
	// Copy values to the GPU
	//
	cudasafe(cudaMemcpy(type_gpu,lens->type , nhalos*sizeof(int),cudaMemcpyHostToDevice ),"Gradient2gpu.cu : Copy type_gpu: " );
	cudasafe(cudaMemcpy(lens_x_gpu,lens->position_x , nhalos*sizeof(type_t),cudaMemcpyHostToDevice ),"Gradient2gpu.cu : Copy x_gpu: " );
	cudasafe(cudaMemcpy(lens_y_gpu,lens->position_y , nhalos*sizeof(type_t), cudaMemcpyHostToDevice),"Gradient2gpu.cu : Copy y_gpu: " );
	cudasafe(cudaMemcpy(b0_gpu,lens->b0 , nhalos*sizeof(type_t), cudaMemcpyHostToDevice),"Gradient2pu.cu : Copy b0_gpu: " );
	cudasafe(cudaMemcpy(angle_gpu,lens->ellipticity_angle , nhalos*sizeof(type_t), cudaMemcpyHostToDevice),"Gradient2gpu.cu : Copy angle_gpu: " );
	cudasafe(cudaMemcpy(epot_gpu, lens->ellipticity_potential, nhalos*sizeof(type_t),cudaMemcpyHostToDevice ),"Gradient2gpu.cu : Copy epot_gpu: " );
	cudasafe(cudaMemcpy(rcore_gpu, lens->rcore, nhalos*sizeof(type_t),cudaMemcpyHostToDevice ),"Gradient2gpu.cu : Copy rcore_gpu: " );
	cudasafe(cudaMemcpy(rcut_gpu, lens->rcut, nhalos*sizeof(type_t), cudaMemcpyHostToDevice),"Gradient2gpu.cu : Copy rcut_gpu: " );
	cudasafe(cudaMemcpy(anglecos_gpu, lens->anglecos, nhalos*sizeof(type_t),cudaMemcpyHostToDevice ),"Gradient2gpu.cu : Copy anglecos: " );
	cudasafe(cudaMemcpy(anglesin_gpu, lens->anglesin, nhalos*sizeof(type_t), cudaMemcpyHostToDevice),"Gradient2gpu.cu : Copy anglesin: " );
	cudasafe(cudaMemcpy(frame_gpu, frame, sizeof(grid_param), cudaMemcpyHostToDevice),"Gradient2gpu.cu : Copy fame_gpu: " );
	//
	lens_gpu->type 			= type_gpu;
	lens_gpu->position_x 		= lens_x_gpu;
	lens_gpu->position_y 		= lens_y_gpu;
	lens_gpu->b0 			= b0_gpu;
	lens_gpu->ellipticity_angle 	= angle_gpu;
	lens_gpu->ellipticity_potential = epot_gpu;
	lens_gpu->rcore 		= rcore_gpu;
	lens_gpu->rcut 			= rcut_gpu;
	lens_gpu->anglecos 		= anglecos_gpu;
	lens_gpu->anglesin 		= anglesin_gpu;
	//
	cudaMemcpy(lens_kernel, lens_gpu, sizeof(Potential_SOA), cudaMemcpyHostToDevice);
	//
	type_t time = -myseconds();
	//module_potentialDerivatives_totalGradient_SOA_CPU_GPU(grid_grad_x_gpu, grid_grad_y_gpu, frame_gpu, lens_kernel, nbgridcells_x, nhalos);
	module_potentialDerivatives_totalGradient2_SOA_CPU_GPU(grid_grad2_a_gpu, grid_grad2_b_gpu, grid_grad2_c_gpu, grid_grad2_d_gpu, frame_gpu, lens_kernel, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
	//
	//cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU");
	cudaDeviceSynchronize();
	time += myseconds();
	//std::cout << "	kernel time = " << time << " s." << std::endl;
	//

	cudasafe(cudaMemcpy( grid_grad2_a, grid_grad2_a_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost )," --- Gradient2gpu.cu : Copy source_a_gpu: " );
	cudasafe(cudaMemcpy( grid_grad2_b, grid_grad2_b_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost)," --- Gradient2gpu.cu : Copy source_b_gpu: " );
	cudasafe(cudaMemcpy( grid_grad2_c, grid_grad2_c_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost )," --- Gradient2gpu.cu : Copy source_c_gpu: " );
	cudasafe(cudaMemcpy( grid_grad2_d, grid_grad2_d_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost)," --- Gradient2gpu.cu : Copy source_d_gpu: " );
	//
	//printf("-----> %f %f \n",grid_grad_x[Nx], grid_grad_y[Ny]);
	// Free GPU memory
	cudaFree(lens_kernel);
	cudaFree(type_gpu);
	cudaFree(lens_x_gpu);
	cudaFree(lens_y_gpu);
	cudaFree(b0_gpu);
	cudaFree(angle_gpu);
	cudaFree(epot_gpu);
	cudaFree(rcore_gpu);
	cudaFree(rcut_gpu);
	cudaFree(anglecos_gpu);
	cudaFree(anglesin_gpu);
	cudaFree(grid_grad2_a_gpu);
	cudaFree(grid_grad2_b_gpu);
	cudaFree(grid_grad2_c_gpu);
	cudaFree(grid_grad2_d_gpu);
}

void 
module_potentialDerivatives_totalGradient2_SOA_CPU_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct grid_param *frame, const struct Potential_SOA *lens_gpu, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
        int GRID_SIZE_X = (nbgridcells_x + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
        int GRID_SIZE_Y = (nbgridcells_y + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y;
        //
        printf("grid_size_x = %d, grid_size_y = %d, nbgridcells_x = %d, nbgridcells_y = %d, istart = %d, jstart = %d (split)\n", GRID_SIZE_X, GRID_SIZE_Y, nbgridcells_x, nbgridcells_y, istart, jstart);
        //
        dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
        dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);
        //printf("nhalos = %d, size of shared memory = %lf\n", nhalos, (double) (8*nhalos + BLOCK_SIZE_X*nbgridcells/BLOCK_SIZE_Y)*sizeof(double));
        printf("nhalos = %d, size of shared memory = %lf (split)\n", nhalos, (type_t) (8*nhalos + BLOCK_SIZE_X*BLOCK_SIZE_Y)*sizeof(type_t));
        //
        cudaMemset(grid_grad2_a, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        cudaMemset(grid_grad2_b, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        cudaMemset(grid_grad2_c, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        cudaMemset(grid_grad2_d, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        //module_potentialDerivatives_totalGradient_SOA_GPU<<<grid, threads>>> (grid_grad_x, grid_grad_y, lens, frame, nhalos, nbgridcells_x);
        module_potentialDerivatives_totalGradient2_SOA_GPU<<<grid, threads>>> (grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,  lens_gpu, frame, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
        cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU_8_SOA_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//
//
void
module_potentialDerivatives_Kmap_SOA_CPU_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct grid_param *frame, const struct Potential_SOA *lens_gpu, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
        int GRID_SIZE_X = (nbgridcells_x + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
        int GRID_SIZE_Y = (nbgridcells_y + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y;
        //
        printf("grid_size_x = %d, grid_size_y = %d, nbgridcells_x = %d, nbgridcells_y = %d, istart = %d, jstart = %d (split)\n", GRID_SIZE_X, GRID_SIZE_Y, nbgridcells_x, nbgridcells_y, istart, jstart);
        //
        dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
        dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);
        //printf("nhalos = %d, size of shared memory = %lf\n", nhalos, (double) (8*nhalos + BLOCK_SIZE_X*nbgridcells/BLOCK_SIZE_Y)*sizeof(double));
        printf("nhalos = %d, size of shared memory = %lf (split)\n", nhalos, (type_t) (8*nhalos + BLOCK_SIZE_X*BLOCK_SIZE_Y)*sizeof(type_t));
        //
        cudaMemset(grid_grad2_a, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        cudaMemset(grid_grad2_b, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        cudaMemset(grid_grad2_c, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        cudaMemset(grid_grad2_d, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        //module_potentialDerivatives_totalGradient_SOA_GPU<<<grid, threads>>> (grid_grad_x, grid_grad_y, lens, frame, nhalos, nbgridcells_x);
        module_potentialDerivatives_Kmap_SOA_GPU<<<grid, threads>>> (grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,  lens_gpu, frame, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
        cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU_8_SOA_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//


