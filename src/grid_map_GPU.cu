//
//
//
#include <fstream>
#include "grid_map_GPU.cuh"
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
//
//
//void amplif_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t z,int mode_amp,int nhalos,int nbgridcells_x, int nbgridcells_y);
//map_gpu_function_t amplif_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t z,int mode_amp,int nhalos,int nbgridcells_x, int nbgridcells_y);

__global__ void amplif_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t z,int nbgridcells);
//
void map_grid_GPU(map_gpu_function_t mapfunction, type_t *map, const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos ,int nbgridcells,int mode_amp, type_t z )
{
	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells - 1);
    type_t dy = (frame->ymax - frame->ymin)/(nbgridcells - 1);
        //
    map_grid_GPU(mapfunction,map,cosmo, frame, lens, nhalos,mode_amp,z, dx, dy, nbgridcells, nbgridcells, 0, 0);
}
//
void map_grid_GPU(map_gpu_function_t mapfunction, type_t *map,const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int mode_amp, type_t z, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{

	int nBlocks_gpu = 0;
	// Define the number of threads per block the GPU will use
	cudaDeviceProp properties_gpu;

	cudaGetDeviceProperties(&properties_gpu, 0); // Get properties of 0th GPU in use

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

	grid_param *frame_gpu;
	Potential_SOA *lens_gpu,*lens_kernel;
	int *type_gpu;
	type_t *lens_x_gpu, *lens_y_gpu, *b0_gpu, *angle_gpu, *epot_gpu, *rcore_gpu, *rcut_gpu, *anglecos_gpu, *anglesin_gpu;
	type_t *grid_grad2_a_gpu, *grid_grad2_b_gpu , *grid_grad2_c_gpu, *grid_grad2_d_gpu, *map_gpu;

	type_t dl0s = module_cosmodistances_objectObject(lens->z[0], z, *cosmo);
	type_t dos = module_cosmodistances_observerObject(z, *cosmo);

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
	mapfunction(map_gpu,grid_grad2_a_gpu, grid_grad2_b_gpu, grid_grad2_c_gpu, grid_grad2_d_gpu,dl0s,z,mode_amp,nhalos,nbgridcells_x,nbgridcells_y);
	//amplif_grid_CPU_GPU(map_gpu,grid_grad2_a_gpu, grid_grad2_b_gpu, grid_grad2_c_gpu, grid_grad2_d_gpu,dl0s,z,mode_amp,nhalos,nbgridcells_x,nbgridcells_y);
	//cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU");
	cudaDeviceSynchronize();
	time += myseconds();
	//std::cout << "	kernel time = " << time << " s." << std::endl;
	//


	cudasafe(cudaMemcpy( map, map_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost )," --- Gradient2gpu.cu : Copy source_a_gpu: " );
	//cudasafe(cudaMemcpy( grid_grad2_b, grid_grad2_b_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost)," --- Gradient2gpu.cu : Copy source_b_gpu: " );
	//cudasafe(cudaMemcpy( grid_grad2_c, grid_grad2_c_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost )," --- Gradient2gpu.cu : Copy source_c_gpu: " );
	//cudasafe(cudaMemcpy( grid_grad2_d, grid_grad2_d_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost)," --- Gradient2gpu.cu : Copy source_d_gpu: " );
	//
	//printf("-----> %f %f \n",grid_grad_x[Nx], grid_grad_y[Ny]);
	// Free GPU memory
	cudaFree(lens_gpu);
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

void amplif5_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t z,int mode_amp, int nhalos,int nbgridcells_x, int nbgridcells_y)
{
        int GRID_SIZE_X = (nbgridcells_x + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
        int GRID_SIZE_Y = (nbgridcells_y + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y;
        //
        //printf("grid_size_x = %d, grid_size_y = %d, nbgridcells_x = %d, nbgridcells_y = %d, istart = %d, jstart = %d (split)\n", GRID_SIZE_X, GRID_SIZE_Y, nbgridcells_x, nbgridcells_y, istart, jstart);
        //
        dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
        dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);
        //
        int count = nhalos;
        //printf("nhalos = %d, size of shared memory = %lf\n", nhalos, (double) (8*nhalos + BLOCK_SIZE_X*nbgridcells/BLOCK_SIZE_Y)*sizeof(double));
        printf("nhalos = %d, size of shared memory = %lf (split)\n", nhalos, (type_t) (8*nhalos + BLOCK_SIZE_X*BLOCK_SIZE_Y)*sizeof(type_t));
        //
        cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        //module_potentialDerivatives_totalGradient_SOA_GPU<<<grid, threads>>> (grid_grad_x, grid_grad_y, lens, frame, nhalos, nbgridcells_x);
        amplif_grid_GPU<<<grid, threads>>> (map,grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,dl0s,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "amplif_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}

__global__
void
amplif_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t z,int nbgridcells)
{
    struct point image_point;
	struct matrix grad, clumpgrad;
	//
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    //
    if ((row < nbgridcells) && (col < nbgridcells))
    {
    	int index = row*nbgridcells + col;
        type_t kappa = (grid_grad2_a[index] + grid_grad2_c[index]) / 2.;
        ampli[index] = kappa;
    }

}

#if 0
void module_potentialDerivatives_totalGradient2_SOA_CPU_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct grid_param *frame, const struct Potential_SOA *lens_gpu, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
        int GRID_SIZE_X = (nbgridcells_x + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
        int GRID_SIZE_Y = (nbgridcells_y + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y;
        //
        printf("grid_size_x = %d, grid_size_y = %d, nbgridcells_x = %d, nbgridcells_y = %d, istart = %d, jstart = %d (split)\n", GRID_SIZE_X, GRID_SIZE_Y, nbgridcells_x, nbgridcells_y, istart, jstart);
        //
        dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
        dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);
        //
        int count = nhalos;
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
//
void
module_potentialDerivatives_totalGradient2_SOA_CPU_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct grid_param *frame, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos)
{
	//
	int GRID_SIZE_X = (nbgridcells + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
	int GRID_SIZE_Y = (nbgridcells + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y; 
	//
	type_t* timer = (type_t *) malloc((int) nbgridcells*nbgridcells*sizeof(type_t));
	type_t* dtimer;
	//
	dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
	dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);	
	//
	int count = nhalos;	
	//printf("nhalos = %d, size of shared memory = %lf\n", nhalos, (type_t) (8*nhalos + BLOCK_SIZE_X*nbgridcells/BLOCK_SIZE_Y)*sizeof(type_t));
	printf("nhalos = %d, size of shared memory = %lf\n", nhalos, (type_t) (8*nhalos + BLOCK_SIZE_X*BLOCK_SIZE_Y)*sizeof(type_t));
	//
	cudaMemset(grid_grad2_a, 0, nbgridcells*nbgridcells*sizeof(type_t));
	cudaMemset(grid_grad2_b, 0, nbgridcells*nbgridcells*sizeof(type_t));
	cudaMemset(grid_grad2_c, 0, nbgridcells*nbgridcells*sizeof(type_t));
	cudaMemset(grid_grad2_d, 0, nbgridcells*nbgridcells*sizeof(type_t));
	//
	//module_potentialDerivatives_totalGradient2_SOA_GPU<<<grid, threads>>> (grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d, lens_gpu, frame, nbgridcells, nhalos);
	cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU_8_SOA_GPU");
	//
	cudaDeviceSynchronize();
	printf("GPU kernel done...\n");
}
#endif

