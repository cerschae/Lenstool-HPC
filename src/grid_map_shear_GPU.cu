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

@brief: Function for shear computation over a grid

*/


#include <fstream>
#include "grid_map_shear_GPU.cuh"
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

__device__ struct  ellipse formeli_shear(type_t a, type_t b, type_t c);


//GPU mapping function declaration to change when figured out linkage problems
__global__ void shear_1_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t dos, type_t z,int nbgridcells);
__global__ void shear_2_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t dos, type_t z,int nbgridcells);
__global__ void shear_3_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t dos, type_t z,int nbgridcells);
__global__ void shear_4_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t dos, type_t z,int nbgridcells);

////Map function selection
map_gpu_function_t select_map_shear_function(const struct runmode_param* runmode){

		if(runmode->shear == 1){
			return &shear_1_grid_CPU_GPU;
		}
		else if(runmode->shear == 2){
			return &shear_2_grid_CPU_GPU;
		}
		else if(runmode->shear == 3){
			return &shear_3_grid_CPU_GPU;
		}
		else if(runmode->shear == 4){
			return &shear_4_grid_CPU_GPU;
		}
		else{
			fprintf(stderr, "ERROR: Mass mode %d not supported yet \n",runmode->mass);
			exit(-1);
		}

	return 0;
}

////General Map calculation
void map_grid_shear_GPU(map_gpu_function_t mapfunction, type_t *map, const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos ,int nbgridcells,int mode_amp, type_t z )
{
	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells - 1);
    type_t dy = (frame->ymax - frame->ymin)/(nbgridcells - 1);
    //
    map_grid_shear_GPU(mapfunction,map,cosmo, frame, lens, nhalos,mode_amp,z, dx, dy, nbgridcells, nbgridcells, 0, 0);
}
//
void map_grid_shear_GPU(map_gpu_function_t mapfunction, type_t *map,const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int mode_amp, type_t z, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{

	grid_param *frame_gpu;
	Potential_SOA *lens_gpu,*lens_kernel;
	int *type_gpu;
	type_t *lens_x_gpu, *lens_y_gpu, *b0_gpu, *angle_gpu, *epot_gpu, *rcore_gpu, *rcut_gpu, *anglecos_gpu, *anglesin_gpu;
	type_t *grid_grad2_a_gpu, *grid_grad2_b_gpu , *grid_grad2_c_gpu, *grid_grad2_d_gpu, *map_gpu;

	type_t dl0s = module_cosmodistances_objectObject(lens->z[0], z, *cosmo);
	type_t dos = module_cosmodistances_observerObject(z, *cosmo);
	type_t dol = module_cosmodistances_observerObject(lens->z[0], *cosmo);
	//select_ratio_function(std::string mode, const struct runmode_param* runmode, type_t dls, type_t ds)

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
	cudasafe(cudaMalloc( (void**)&(map_gpu), (nbgridcells_x) * (nbgridcells_y) *sizeof(type_t)),"Gradient2gpu.cu : Alloc map: " );
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
	//
	module_potentialDerivatives_totalGradient2_SOA_CPU_GPU(grid_grad2_a_gpu, grid_grad2_b_gpu, grid_grad2_c_gpu, grid_grad2_d_gpu, frame_gpu, lens_kernel, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
	//
	mapfunction(map_gpu,grid_grad2_a_gpu, grid_grad2_b_gpu, grid_grad2_c_gpu, grid_grad2_d_gpu,dl0s,dos,dol,cosmo->h,z,nbgridcells_x,nbgridcells_y,frame);
	//cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU");
	cudaDeviceSynchronize();
	//
	cudasafe(cudaMemcpy( map, map_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost )," --- Gradient2gpu.cu : Copy source_a_gpu: " );
	//
	time += myseconds();
	std::cout << "	kernel time = " << time << " s." << std::endl;
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
	cudaFree(map_gpu);
}

//Shear NR 1
void shear_1_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame)
{
        int GRID_SIZE_X = (nbgridcells_x + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
        int GRID_SIZE_Y = (nbgridcells_y + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y;
        //
        //printf("grid_size_x = %d, grid_size_y = %d, nbgridcells_x = %d, nbgridcells_y = %d, istart = %d, jstart = %d (split)\n", GRID_SIZE_X, GRID_SIZE_Y, nbgridcells_x, nbgridcells_y, istart, jstart);
        //
        dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
        dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);
        //
        cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        shear_1_grid_GPU<<<grid, threads>>> (map,grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,dl0s,ds,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "amplif_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//
__global__ void shear_1_grid_GPU(type_t *shear,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t z,int nbgridcells)
{
	type_t A,B,C;
	ellipse amp;
	type_t dlsds= dl0s/ds;
	////
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    //
    if ((row < nbgridcells) && (col < nbgridcells))
    {

    	int index = row*nbgridcells + col;
        A = 1. - grid_grad2_a[index]*dlsds;   // 1 - DLS/DS * d2phixx
        B = - grid_grad2_b[index]*dlsds;   // - DLS/DS * d2phixy
        C = 1. - grid_grad2_c[index]*dlsds;   // 1 - DLS/DS * d2phiyy
        amp = formeli_shear(A, B, C);
        shear[index] = (amp.a - amp.b) /2. ;
    }
}
//Shear NR 2
void shear_2_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame)
{
        int GRID_SIZE_X = (nbgridcells_x + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
        int GRID_SIZE_Y = (nbgridcells_y + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y;
        //
        //printf("grid_size_x = %d, grid_size_y = %d, nbgridcells_x = %d, nbgridcells_y = %d, istart = %d, jstart = %d (split)\n", GRID_SIZE_X, GRID_SIZE_Y, nbgridcells_x, nbgridcells_y, istart, jstart);
        //
        dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
        dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);
        //
        cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        shear_2_grid_GPU<<<grid, threads>>> (map,grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,dl0s,ds,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "amplif_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//
__global__ void shear_2_grid_GPU(type_t *shear,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t z,int nbgridcells)
{
	type_t A,B,C;
	ellipse amp;
	type_t dlsds= dl0s/ds;
	type_t q;
	////
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    //
    if ((row < nbgridcells) && (col < nbgridcells))
    {

    	int index = row*nbgridcells + col;
        A = 1. - grid_grad2_a[index]*dlsds;   // 1 - DLS/DS * d2phixx
        B = - grid_grad2_b[index]*dlsds;   // - DLS/DS * d2phixy
        C = 1. - grid_grad2_c[index]*dlsds;   // 1 - DLS/DS * d2phiyy
        amp = formeli_shear(A, B, C);
        q = amp.a / amp.b;
        shear[index] = fabs((q * q - 1.) / (q * q + 1.));
    }
}
//Shear NR 3
void shear_3_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame)
{
        int GRID_SIZE_X = (nbgridcells_x + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
        int GRID_SIZE_Y = (nbgridcells_y + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y;
        //
        //printf("grid_size_x = %d, grid_size_y = %d, nbgridcells_x = %d, nbgridcells_y = %d, istart = %d, jstart = %d (split)\n", GRID_SIZE_X, GRID_SIZE_Y, nbgridcells_x, nbgridcells_y, istart, jstart);
        //
        dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
        dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);
        //
        cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        shear_3_grid_GPU<<<grid, threads>>> (map,grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,dl0s,ds,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "amplif_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//
__global__ void shear_3_grid_GPU(type_t *shear,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t z,int nbgridcells)
{
	type_t dlsds= dl0s/ds;
	////
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    //
    if ((row < nbgridcells) && (col < nbgridcells))
    {

    	int index = row*nbgridcells + col;
        shear[index] = (grid_grad2_a[index] - grid_grad2_c[index])*dlsds / 2. ;
        //if(col == 0 and row == 0)printf(" GPUUUU:  Grad  %f %f %f  ABC %f %f %f dlsds %f ab %f %f %f \n",grid_grad2_a[index]*dlsds,grid_grad2_b[index]*dlsds, grid_grad2_c[index]*dlsds,A,B,C,dlsds,amp.a,amp.b,(A - C)*(A - C) + 4*B*B);
    }
}
//Shear NR 4
void shear_4_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame)
{
        int GRID_SIZE_X = (nbgridcells_x + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
        int GRID_SIZE_Y = (nbgridcells_y + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y;
        //
        //printf("grid_size_x = %d, grid_size_y = %d, nbgridcells_x = %d, nbgridcells_y = %d, istart = %d, jstart = %d (split)\n", GRID_SIZE_X, GRID_SIZE_Y, nbgridcells_x, nbgridcells_y, istart, jstart);
        //
        dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
        dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);
        //
        //printf("nhalos = %d, size of shared memory = %lf (split)\n", nhalos, (type_t) (8*nhalos + BLOCK_SIZE_X*BLOCK_SIZE_Y)*sizeof(type_t));
        //
        //std::cerr << "Amplif 1 :"  <<dl0s<<" " <<ds<<" " << std::endl;
        //
        cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        shear_4_grid_GPU<<<grid, threads>>> (map,grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,dl0s,ds,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "amplif_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//
__global__ void shear_4_grid_GPU(type_t *shear,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t z,int nbgridcells)
{
	type_t dlsds= dl0s/ds;
	////
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    //
    if ((row < nbgridcells) && (col < nbgridcells))
    {

    	int index = row*nbgridcells + col;
    	shear[index] = grid_grad2_b[index] * dlsds ;
        //if(col == 0 and row == 0)printf(" GPUUUU:  Grad  %f %f %f  ABC %f %f %f dlsds %f ab %f %f %f \n",grid_grad2_a[index]*dlsds,grid_grad2_b[index]*dlsds, grid_grad2_c[index]*dlsds,A,B,C,dlsds,amp.a,amp.b,(A - C)*(A - C) + 4*B*B);
    }
}

__device__ struct  ellipse formeli_shear(type_t a, type_t b, type_t c)
{
    struct  ellipse eli;
    type_t  e, delta, lambda, mu;

    // eq carateristique : det(M-xI) = 0
    delta = (a - c)*(a - c) + 4*b*b;    //  4*gamma^2 (cf phd_JPK eq 2.56)
    e = sqrt(delta);  /*e is 2 * shear, ie 2*gamma*/
    lambda = .5*(a + c + e);  // 1 - k + gamma
    mu = .5*(a + c - e);      // 1 - k - gamma

    eli.a = lambda;
    eli.b = mu;
    if (lambda != mu && fabs(b) > 1e-5)
        eli.theta = atan2(lambda - a, b); // cf phd_JPK eq 2.58, and
    // tan(theta)= ( -cos(2theta) +- 1 ) / sin(2theta)
// ADDED by EJ 29/11/2007
    else if ( a >= c ) // ellipse aligned along the major axis of magnification
        eli.theta = 0.;
    else
        eli.theta = acos(-1.) / 2.;    // ellipse aligned along the minor axis of magnification

    return(eli);
}

