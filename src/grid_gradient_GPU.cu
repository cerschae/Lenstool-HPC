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

@brief: Function for first order derivative computation over a grid

*/
#include <cuda.h>
#include <fstream>
#include "grid_gradient_GPU.cuh"
#include "gradient_GPU.cuh"

#define BLOCK_SIZE_X 16
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


//#define cudasafe 


void cudasafe( cudaError_t error, const char* message)
  {
     if(error!=cudaSuccess) { fprintf(stderr,"ERROR: %s : %i\n",message,error); exit(-1); }
  }


extern "C" 
{
	type_t myseconds();
}

__global__ 
void module_potentialDerivatives_totalGradient_SOA_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int nhalos);

//
void 
module_potentialDerivatives_totalGradient_SOA_CPU_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart);
//
void
module_potentialDerivatives_totalGradient_SOA_CPU_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos);
//
void gradient_grid_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart);

//
//void
//module_potentialDerivatives_totalGradient_SOA_CPU_GPU_v2(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_cpu, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos);
//
//
//
void gradient_grid_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos ,int nbgridcells)
{
	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells - 1);
        type_t dy = (frame->ymax - frame->ymin)/(nbgridcells - 1);
        //
        gradient_grid_GPU(grid_grad_x, grid_grad_y, frame, lens, nhalos, dx, dy, nbgridcells, nbgridcells, 0, 0);
}
//
//
//
void gradient_grid_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
	grid_param *frame_gpu;
	Potential_SOA *lens_gpu,*lens_kernel;
	int *type_gpu, *N_types_gpu;
	type_t *lens_x_gpu, *lens_y_gpu, *b0_gpu, *angle_gpu, *epot_gpu, *rcore_gpu, *rcut_gpu, *anglecos_gpu, *anglesin_gpu;

	lens_gpu = (Potential_SOA *) malloc(sizeof(Potential_SOA));
	lens_gpu->type = (int *) malloc(sizeof(int));
	// Allocate variables on the GPU
	cudasafe(cudaMalloc( (void**)&(lens_kernel), sizeof(Potential_SOA)),"Gradientgpu.cu : Alloc Potential_SOA: " );
	cudasafe(cudaMalloc( (void**)&(type_gpu), nhalos*sizeof(int)),"Gradientgpu.cu : Alloc type_gpu: " );
	cudasafe(cudaMalloc( (void**)&(lens_x_gpu), nhalos*sizeof(type_t)),"Gradientgpu.cu : Alloc x_gpu: " );
	cudasafe(cudaMalloc( (void**)&(lens_y_gpu), nhalos*sizeof(type_t)),"Gradientgpu.cu : Alloc y_gpu: " );
	cudasafe(cudaMalloc( (void**)&(b0_gpu), nhalos*sizeof(type_t)),"Gradientgpu.cu : Alloc b0_gpu: " );
	cudasafe(cudaMalloc( (void**)&(angle_gpu), nhalos*sizeof(type_t)),"Gradientgpu.cu : Alloc angle_gpu: " );
	cudasafe(cudaMalloc( (void**)&(epot_gpu), nhalos*sizeof(type_t)),"Gradientgpu.cu : Alloc epot_gpu: " );
	cudasafe(cudaMalloc( (void**)&(rcore_gpu), nhalos*sizeof(type_t)),"Gradientgpu.cu : Alloc rcore_gpu: " );
	cudasafe(cudaMalloc( (void**)&(rcut_gpu), nhalos*sizeof(type_t)),"Gradientgpu.cu : Alloc rcut_gpu: " );
	cudasafe(cudaMalloc( (void**)&(anglecos_gpu), nhalos*sizeof(type_t)),"Gradientgpu.cu : Alloc anglecos_gpu: " );
	cudasafe(cudaMalloc( (void**)&(anglesin_gpu), nhalos*sizeof(type_t)),"Gradientgpu.cu : Alloc anglesin_gpu: " );
	cudasafe(cudaMalloc( (void**)&(frame_gpu), sizeof(grid_param)),"Gradientgpu.cu : Alloc frame_gpu: " );
	cudasafe(cudaMalloc( (void**)&(N_types_gpu), 100*sizeof(int)),"Gradientgpu.cu : Alloc N_types: " );
	//
	// Copy values to the GPU
	//
	cudasafe(cudaMemcpy(type_gpu,lens->type , nhalos*sizeof(int),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy type_gpu: " );
	cudasafe(cudaMemcpy(lens_x_gpu,lens->position_x , nhalos*sizeof(type_t),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy x_gpu: " );
	cudasafe(cudaMemcpy(lens_y_gpu,lens->position_y , nhalos*sizeof(type_t), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy y_gpu: " );
	cudasafe(cudaMemcpy(b0_gpu,lens->b0 , nhalos*sizeof(type_t), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy b0_gpu: " );
	cudasafe(cudaMemcpy(angle_gpu,lens->ellipticity_angle , nhalos*sizeof(type_t), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy angle_gpu: " );
	cudasafe(cudaMemcpy(epot_gpu, lens->ellipticity_potential, nhalos*sizeof(type_t),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy epot_gpu: " );
	cudasafe(cudaMemcpy(rcore_gpu, lens->rcore, nhalos*sizeof(type_t),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy rcore_gpu: " );
	cudasafe(cudaMemcpy(rcut_gpu, lens->rcut, nhalos*sizeof(type_t), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy rcut_gpu: " );
	cudasafe(cudaMemcpy(anglecos_gpu, lens->anglecos, nhalos*sizeof(type_t),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy anglecos: " );
	cudasafe(cudaMemcpy(anglesin_gpu, lens->anglesin, nhalos*sizeof(type_t), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy anglesin: " );
	cudasafe(cudaMemcpy(N_types_gpu, lens->N_types, 100*sizeof(int), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy anglesin: " );
	cudasafe(cudaMemcpy(frame_gpu, frame, sizeof(grid_param), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy frame_gpu: " );
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
	lens_gpu->N_types 		= N_types_gpu;
	//
	//cudaMemcpy(lens_kernel, lens_gpu, sizeof(Potential_SOA), cudaMemcpyHostToDevice);
	//
        cudasafe(cudaMalloc( (void**)&(frame_gpu), sizeof(grid_param)),"Gradientgpu.cu : Alloc frame_gpu: " );
        cudasafe(cudaMemcpy(frame_gpu, frame, sizeof(grid_param), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy fame_gpu: " );
        cudasafe(cudaMemcpy(lens_kernel, lens_gpu, sizeof(Potential_SOA), cudaMemcpyHostToDevice), "Gradientgpu.cu: Copy lens_kernel");
	//
	type_t *grid_grad_x_gpu, *grid_grad_y_gpu ;
	cudasafe(cudaMalloc( (void**)&(grid_grad_x_gpu), (nbgridcells_x) * (nbgridcells_y) *sizeof(type_t)),"Gradientgpu.cu : Alloc source_x_gpu: " );
	cudasafe(cudaMalloc( (void**)&(grid_grad_y_gpu), (nbgridcells_x) * (nbgridcells_y) *sizeof(type_t)),"Gradientgpu.cu : Alloc source_y_gpu: " );
	//
	type_t time = -myseconds();
	//module_potentialDerivatives_totalGradient_SOA_CPU_GPU(grid_grad_x_gpu, grid_grad_y_gpu, frame_gpu, lens_kernel, nbgridcells_x, nhalos);
	//@@module_potentialDerivatives_totalGradient_SOA_CPU_GPU(grid_grad_x_gpu, grid_grad_y_gpu, frame_gpu, lens_kernel, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
	module_potentialDerivatives_totalGradient_SOA_CPU_GPU(grid_grad_x_gpu, grid_grad_y_gpu, frame_gpu, lens_kernel, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
	//
	cudaDeviceSynchronize();
	cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU");
	time += myseconds();
	//std::cout << "	kernel time = " << time << " s." << std::endl;
	//

	cudasafe(cudaMemcpy( grid_grad_x, grid_grad_x_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost),"Gradientgpu.cu : Copy source_x_gpu D2H" );
	cudasafe(cudaMemcpy( grid_grad_y, grid_grad_y_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost),"Gradientgpu.cu : Copy source_y_gpu D2H" );
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
	cudaFree(frame_gpu);
	cudaFree(N_types_gpu);
	cudaFree(grid_grad_x_gpu);
	cudaFree(grid_grad_y_gpu);
}


#if 0
void gradient_grid_GPU_UM(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
        type_t *grid_grad_x_gpu, *grid_grad_y_gpu ;
        //
        cudasafe(cudaMalloc( (void**)&(grid_grad_x_gpu), (nbgridcells_x) * (nbgridcells_y) *sizeof(type_t)),"Gradientgpu.cu : Alloc source_x_gpu: " );
        cudasafe(cudaMalloc( (void**)&(grid_grad_y_gpu), (nbgridcells_x) * (nbgridcells_y) *sizeof(type_t)),"Gradientgpu.cu : Alloc source_y_gpu: " );
        //
        grid_param *frame_gpu;
        cudasafe(cudaMalloc( (void**)&(frame_gpu), sizeof(grid_param)),"Gradientgpu.cu : Alloc frame_gpu: " );
        cudasafe(cudaMemcpy(frame_gpu, frame, sizeof(grid_param), cudaMemcpyHostToDevice),"UM Gradientgpu.cu : Copy frame_gpu: " );
        //
        type_t time = -myseconds();
        //module_potentialDerivatives_totalGradient_SOA_CPU_GPU(grid_grad_x_gpu, grid_grad_y_gpu, frame_gpu, lens_kernel, nbgridcells_x, nhalos);
        //@@module_potentialDerivatives_totalGradient_SOA_CPU_GPU(grid_grad_x_gpu, grid_grad_y_gpu, frame_gpu, lens_kernel, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
        module_potentialDerivatives_totalGradient_SOA_CPU_GPU(grid_grad_x_gpu, grid_grad_y_gpu, frame_gpu, lens, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);

        //
        //cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU");
        cudaDeviceSynchronize();
        time += myseconds();
        std::cout << "        kernel time (UM) = " << time << " s." << std::endl; fflush(stdout);
        //
        cudasafe(cudaMemcpy( grid_grad_x, grid_grad_x_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost )," --- Gradientgpu.cu : Copy source_x_gpu H2D " );
        cudasafe(cudaMemcpy( grid_grad_y, grid_grad_y_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost)," --- Gradientgpu.cu : Copy source_y_gpu H2D " );
        //
        cudaFree(grid_grad_x_gpu);
        cudaFree(grid_grad_y_gpu);
        cudaFree(frame_gpu);
}
#endif



void gradient_grid_GPU_UM(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
	//type_t *grid_grad_x_gpu, *grid_grad_y_gpu ;
	//
        //cudasafe(cudaMalloc( (void**)&(grid_grad_x_gpu), (nbgridcells_x) * (nbgridcells_y) *sizeof(type_t)),"Gradientgpu.cu : Alloc source_x_gpu: " );
        //cudasafe(cudaMalloc( (void**)&(grid_grad_y_gpu), (nbgridcells_x) * (nbgridcells_y) *sizeof(type_t)),"Gradientgpu.cu : Alloc source_y_gpu: " );
	//
	grid_param *frame_gpu;
	cudasafe(cudaMalloc( (void**)&(frame_gpu), sizeof(grid_param)),"Gradientgpu.cu : Alloc frame_gpu: " );
	cudasafe(cudaMemcpy(frame_gpu, frame, sizeof(grid_param), cudaMemcpyHostToDevice),"UM Gradientgpu.cu : Copy frame_gpu: " );
        //
        type_t time = -myseconds();
        //module_potentialDerivatives_totalGradient_SOA_CPU_GPU(grid_grad_x_gpu, grid_grad_y_gpu, frame_gpu, lens_kernel, nbgridcells_x, nhalos);
        //@@module_potentialDerivatives_totalGradient_SOA_CPU_GPU(grid_grad_x_gpu, grid_grad_y_gpu, frame_gpu, lens_kernel, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
        module_potentialDerivatives_totalGradient_SOA_CPU_GPU(grid_grad_x, grid_grad_y, frame_gpu, lens, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
	cudasafe(cudaGetLastError(), "module_potentialDerivatives_totalGradient_SOA_CPU_GPU");

        //
        //cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU");
        cudaDeviceSynchronize();
        time += myseconds();
        std::cout << "        kernel time (UM) = " << time << " s." << std::endl; fflush(stdout);
        //
        //cudasafe(cudaMemcpy( grid_grad_x, grid_grad_x_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost )," --- Gradientgpu.cu : Copy source_x_gpu H2D " );
        //cudasafe(cudaMemcpy( grid_grad_y, grid_grad_y_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost)," --- Gradientgpu.cu : Copy source_y_gpu H2D " );
	//
	//cudaFree(grid_grad_x_gpu);
	//cudaFree(grid_grad_y_gpu);
	cudaFree(frame_gpu);
}
//
#if 0
__global__
void
module_potentialDerivatives_totalGradient_8_SOA_GPU_loc(type_t *grid_grad_x, type_t *grid_grad_y, const struct Potential_SOA
 *lens, const struct grid_param *frame, int nbgridcells, int shalos, int nhalos)
{
        //
        //asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        //
        type_t grad_x, grad_y;
        type_t clumpgrad_x, clumpgrad_y;
        type_t image_point_x, image_point_y;
        //
        SHARED type_t cosi      [200];
        SHARED type_t sinu      [200];
        SHARED type_t rci       [200];
        SHARED type_t b0        [200];
        SHARED type_t epsi      [200];
        SHARED type_t position_x[200];
        SHARED type_t position_y[200];
        SHARED type_t rsqe      [200];
        //SHARED type_t sgrad_x   [(BLOCK_SIZE_X + 1)*BLOCK_SIZE_Y];
        //SHARED type_t sgrad_y   [(BLOCK_SIZE_X + 1)*BLOCK_SIZE_Y];

        //SHARED type_t sonepeps  [200];
        //SHARED type_t sonemeps  [200];
        //
        grad_x         = 0;
        grad_y         = 0;
        //
        int row        = blockIdx.y*blockDim.y + threadIdx.y;
        int col        = blockIdx.x*blockDim.x + threadIdx.x;
        //
        //int loc_row    = threadIdx.x;
        //int loc_col    = threadIdx.y*blockDim.x + threadIdx.x;
        //
        //int grid_size  = nbgridcells/blockDim.y;
        //
        //if (threadIdx.x == 0) printf("%d %d %d: row = %d, col = %d, grid_size = %d\n", blockIdx.y, gridDim.y, threadIdx.y, row, col, grid_size);
        //
        type_t dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
        type_t dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
        //if (threadIdx.x == 0) printf("dx = %f, dy = %f\n", dx, dy);
        //
        image_point_x = frame->xmin + col*dx;
        image_point_y = frame->ymin + row*dy;
        //
        //int iloc  = threadIdx.x*blockDim.y + threadIdx.y;
        int iglob = row*nbgridcells + col;
        int numThreads = blockDim.x*blockDim.y;
        //
        for (int i = 0; i < (nhalos + numThreads - 1)/numThreads; ++i)
        {
                int iloc  = threadIdx.y*blockDim.x + threadIdx.x + i*numThreads;
                if (iloc < nhalos)
                {
                        cosi[iloc]       = __ldg(&lens->anglecos             [shalos + iloc]);
                        sinu[iloc]       = __ldg(&lens->anglesin             [shalos + iloc]);
                        position_x[iloc] = __ldg(&lens->position_x           [shalos + iloc]);
                        position_y[iloc] = __ldg(&lens->position_y           [shalos + iloc]);
                        rci[iloc]        = __ldg(&lens->rcore                [shalos + iloc]);
                        b0[iloc]         = __ldg(&lens->b0                   [shalos + iloc]);
                        epsi[iloc]       = __ldg(&lens->ellipticity_potential[shalos + iloc]);
                        rsqe[iloc]       = sqrt(epsi[iloc]);
                }
        }
        __syncthreads();
        //
        if ((row < nbgridcells) && (col < nbgridcells))
        {
                //
                for(int i = 0; i < nhalos; i++)
                {
                        //int index  = iloc;
#if 1
                        type_t rc      = rci[i];
                        type_t eps     = epsi[i];
                        type_t onemeps = 1 - eps;
                        type_t onepeps = 1 + eps;
                        //
                        type_t sqe = rsqe[i];
                        type_t cx1 = onemeps/onepeps;
                        //
                        //
                        //type_t zci_im = 1;
                        type_t zci_im  = -0.5*(1. - eps*eps)/sqe;
                        type_t inv_onepeps = 1./(onepeps*onepeps);
                        type_t inv_onemeps = 1./(onemeps*onemeps);
#endif
                        //
                        {
                                //KERNEL_8;
                                type_t grad_x = grad_y = 0.;
                                type_t true_coord_y = image_point_y - position_y[i];
                                type_t true_coord_x = image_point_x - position_x[i];
                                KERNEL_8_reg(0);
                                grid_grad_x[iglob +  0] += grad_x;
                                grid_grad_y[iglob +  0] += grad_y;
                        }
                        /*
                        {
                                //KERNEL_8;
                                type_t grad_x = grad_y = 0.;
                                type_t true_coord_y = image_point_y - position_y[i];
                                type_t true_coord_x = image_point_x - position_x[i] + BLOCK_SIZE_X/2*dx;
                                KERNEL_8_reg(0);
                                grid_grad_x[iglob + BLOCK_SIZE_X/2] += grad_x;
                                grid_grad_y[iglob + BLOCK_SIZE_X/2] += grad_y;
                        }
                        */
                //
                }
        }
}


#endif
void 
module_potentialDerivatives_totalGradient_SOA_CPU_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_gpu, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
	//printf("module_potentialDerivatives_totalGradient_SOA_CPU_GPU-1\n");
        //
        int GRID_SIZE_X = (nbgridcells_x + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
        int GRID_SIZE_Y = (nbgridcells_y + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y;
	//
        printf("grid_size_x = %d, grid_size_y = %d, nbgridcells_x = %d, nbgridcells_y = %d, istart = %d, jstart = %d (split)\n", GRID_SIZE_X, GRID_SIZE_Y, nbgridcells_x, nbgridcells_y, istart, jstart);
        //
        dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
        dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);
        //
        //@@printf("nhalos = %d, size of shared memory = %lf (split)\n", nhalos, (type_t) (8*nhalos + BLOCK_SIZE_X*BLOCK_SIZE_Y)*sizeof(type_t));
        //
        cudaMemset(grid_grad_x, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        cudaMemset(grid_grad_y, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
	//
        //module_potentialDerivatives_totalGradient_SOA_GPU<<<grid, threads>>> (grid_grad_x, grid_grad_y, lens, frame, nhalos, nbgridcells_x);
        module_potentialDerivatives_totalGradient_SOA_GPU<<<grid, threads>>> (grid_grad_x, grid_grad_y, lens_gpu, frame, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
        cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU-1");
	//
        cudaDeviceSynchronize();
        //@@printf("GPU kernel done...\n");
}
//
//
//
void
module_potentialDerivatives_totalGradient_SOA_CPU_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos)
{
	//printf("module_potentialDerivatives_totalGradient_SOA_CPU_GPU-2\n");
	//
	int GRID_SIZE_X = (nbgridcells + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
	int GRID_SIZE_Y = (nbgridcells + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y; 
	//printf("grid_size_x = %d, grid_size_y = %d\n", GRID_SIZE_X, GRID_SIZE_Y);
	//
	type_t* timer = (type_t *) malloc((int) nbgridcells*nbgridcells*sizeof(type_t));
	//type_t* dtimer;
	//cudasafe(cudaMalloc( (void**)&(dtimer), (int) nbgridcells*nbgridcells*sizeof(type_t)),"Gradientgpu.cu : totalGradient_SOA_CPU_GPU: " );
	//
	//{
	//dim3 dimBlock(BLOCK_SIZE_X/1, BLOCK_SIZE_Y);
	//dim3 dimGrid (GRID_SIZE_X   , GRID_SIZE_Y );	
	dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
	dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);	
	//
	//int count = nhalos;
	//printf("nhalos = %d, size of shared memory = %lf\n", nhalos, (type_t) (8*nhalos + BLOCK_SIZE_X*nbgridcells/BLOCK_SIZE_Y)*sizeof(type_t));
	//@@printf("nhalos = %d, size of shared memory = %lf\n", nhalos, (type_t) (8*nhalos + BLOCK_SIZE_X*BLOCK_SIZE_Y)*sizeof(type_t));
	//
	cudaMemset(grid_grad_x, 0, nbgridcells*nbgridcells*sizeof(type_t));
	cudaMemset(grid_grad_y, 0, nbgridcells*nbgridcells*sizeof(type_t));
		//
	module_potentialDerivatives_totalGradient_SOA_GPU<<<grid, threads>>> (grid_grad_x, grid_grad_y, lens_gpu, frame, nbgridcells, nhalos);
	cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU-2");
		
	//}
	cudaDeviceSynchronize();
	//@@printf("GPU kernel done...\n");
}


