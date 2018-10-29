/**
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)
* @date   July 2017
* @version 0,1
*
*/

#include <fstream>
#include "grid_map_dpl_GPU.cuh"
#include "gradient_GPU.cuh"
//#include "gradient_GPU.cu"
#include <structure_hpc.hpp>

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


#define cudasafe 

extern "C" 
{
	type_t myseconds();
}



//GPU mapping function declaration to change when figured out linkage problems
__global__ void dpl_grid_GPU(type_t *dpl, type_t dl0s, type_t ds, type_t z,int nbgridcells);
////Map function selection
#if 0
map_gpu_function_t select_map_dpl_function(const struct runmode_param* runmode){

		if(runmode->dpl == 1){
			return &dpl_1_grid_CPU_GPU;
		}
		else{
			fprintf(stderr, "ERROR: Mass mode %d not supported yet \n",runmode->mass);
			exit(-1);
		}

	return 0;
}
#endif

////General Map calculation
void map_grid_dpl_GPU(map_gpu_function_t mapfunction, type_t *grid_grad_x, type_t *grid_grad_y, const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos ,int nbgridcells,int mode_amp, type_t z )
{
	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells - 1);
    type_t dy = (frame->ymax - frame->ymin)/(nbgridcells - 1);
    //
    map_grid_dpl_GPU(mapfunction, grid_grad_x, grid_grad_y, cosmo, frame, lens, nhalos,mode_amp,z, dx, dy, nbgridcells, nbgridcells, 0, 0);
}
//
void map_grid_dpl_GPU(map_gpu_function_t mapfunction, type_t *grid_grad_x, type_t *grid_grad_y, const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int mode_amp, type_t z, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
	grid_param *frame_gpu;
	Potential_SOA *lens_gpu,*lens_kernel;
	int *type_gpu;
	type_t *lens_x_gpu, *lens_y_gpu, *b0_gpu, *angle_gpu, *epot_gpu, *rcore_gpu, *rcut_gpu, *anglecos_gpu, *anglesin_gpu;
	type_t *grid_grad_x_gpu, *grid_grad_y_gpu ;

	type_t dl0s = module_cosmodistances_objectObject(lens->z[0], z, *cosmo);
	type_t dos = module_cosmodistances_observerObject(z, *cosmo);
	type_t dol = module_cosmodistances_observerObject(lens->z[0], *cosmo);

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
	cudasafe(cudaMalloc( (void**)&(grid_grad_x_gpu), (nbgridcells_x) * (nbgridcells_y) *sizeof(type_t)),"Gradientgpu.cu : Alloc source_x_gpu: " );
	cudasafe(cudaMalloc( (void**)&(grid_grad_y_gpu), (nbgridcells_x) * (nbgridcells_y) *sizeof(type_t)),"Gradientgpu.cu : Alloc source_y_gpu: " );
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
	//
	cudasafe(cudaMemcpy(lens_kernel, lens_gpu, sizeof(Potential_SOA), cudaMemcpyHostToDevice), "Gradientgpu.cu: Copy lens_kernel");
	//
	type_t time = -myseconds();
	module_potentialDerivatives_totalGradient_SOA_CPU_GPU(grid_grad_x_gpu, grid_grad_y_gpu, frame_gpu, lens_kernel, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
	//
	dpl_grid_CPU_GPU(grid_grad_x_gpu,dl0s,dos,dol,cosmo->h,z,nbgridcells_x,nbgridcells_y,frame);
	dpl_grid_CPU_GPU(grid_grad_y_gpu,dl0s,dos,dol,cosmo->h,z,nbgridcells_x,nbgridcells_y,frame);
	//cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU");
	cudaDeviceSynchronize();
	time += myseconds();
	//std::cout << "	kernel time = " << time << " s." << std::endl;
	//

	cudasafe(cudaMemcpy( grid_grad_x, grid_grad_x_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost )," --- Gradientgpu.cu : Copy source_x_gpu: " );
	cudasafe(cudaMemcpy( grid_grad_y, grid_grad_y_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost)," --- Gradientgpu.cu : Copy source_y_gpu: " );
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
	cudaFree(grid_grad_x_gpu);
	cudaFree(grid_grad_y_gpu);
}

////Map functions
//DPL NR 1
void dpl_grid_CPU_GPU(type_t *map, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame)
{
        int GRID_SIZE_X = (nbgridcells_x + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
        int GRID_SIZE_Y = (nbgridcells_y + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y;
        //
        //printf("grid_size_x = %d, grid_size_y = %d, nbgridcells_x = %d, nbgridcells_y = %d, istart = %d, jstart = %d (split)\n", GRID_SIZE_X, GRID_SIZE_Y, nbgridcells_x, nbgridcells_y, istart, jstart);
        //
        dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
        dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);
        //
        //cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        dpl_grid_GPU<<<grid, threads>>> (map,dl0s,ds,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "dpl_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//
__global__ void dpl_grid_GPU(type_t *dpl, type_t dl0s, type_t ds, type_t z,int nbgridcells)
{
	//ellipse amp;
	type_t dlsds= dl0s/ds;
	////
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    //
    if ((row < nbgridcells) && (col < nbgridcells))
    {

    	int index = row*nbgridcells + col;
        dpl[index] = dpl[index] * dlsds;
        //if(col == 0 and row == 0)printf(" GPUUUU:  Grad  %f %f %f  ABC %f %f %f dlsds %f ab %f %f %f \n",grid_grad2_a[index]*dlsds,grid_grad2_b[index]*dlsds, grid_grad2_c[index]*dlsds,A,B,C,dlsds,amp.a,amp.b,(A - C)*(A - C) + 4*B*B);
    }
}


