/**
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)
* @date   July 2017
* @version 0,1
*
*/

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

__device__ struct  ellipse formeli_HPC(type_t a, type_t b, type_t c);


//GPU mapping function declaration to change when figured out linkage problems
__global__ void amplif_1_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t dos, type_t z,int nbgridcells);
__global__ void amplif_2_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t dos, type_t z,int nbgridcells);
__global__ void amplif_3_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t dos, type_t z,int nbgridcells);
__global__ void amplif_4_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t dos, type_t z,int nbgridcells);
__global__ void amplif_5_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t dos, type_t z,int nbgridcells);
__global__ void amplif_6_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t dos, type_t z,int nbgridcells);
__global__ void mass_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t mult, type_t ds, type_t dl, type_t h, type_t z,int nbgridcells);


////Map function selection
map_gpu_function_t select_map_function(std::string mode, const struct runmode_param* runmode){
	if (mode == "ampli"){
		if(runmode->amplif == 1){
			return &amplif_1_grid_CPU_GPU;
		}
		else if(runmode->amplif == 2){
			return &amplif_2_grid_CPU_GPU;
		}
		else if(runmode->amplif == 3){
			return &amplif_3_grid_CPU_GPU;
		}
		else if(runmode->amplif == 4){
			return &amplif_4_grid_CPU_GPU;
		}
		else if(runmode->amplif == 5){
			return &amplif_5_grid_CPU_GPU;
		}
		else if(runmode->amplif == 6){
			return &amplif_6_grid_CPU_GPU;
		}
		else{
			fprintf(stderr, "ERROR: Amplif mode %d not supported yet \n",runmode->amplif);
			exit(-1);
		}
	}
	else if(mode == "mass"){
		if(runmode->mass == 1){
			return &mass_1_grid_CPU_GPU;
		}
		else if(runmode->mass == 2){
			return &mass_2_grid_CPU_GPU;
		}
		else if(runmode->mass == 3){
			return &mass_3_grid_CPU_GPU;
		}
		else if(runmode->mass == 4){
			return &mass_4_grid_CPU_GPU;
		}
		else{
			fprintf(stderr, "ERROR: Mass mode %d not supported yet \n",runmode->mass);
			exit(-1);
		}
	}
	else{
		fprintf(stderr, "ERROR: No mode recognised \n");
		exit(-1);
	}
	return 0;
}

////Mass Map calculation, doesnt fit the bloody template...
void map_mass_grid_GPU(map_gpu_function_t mapfunction, type_t *map, const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos ,int nbgridcells,int mode_amp, type_t zl, type_t zs )
{
	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells - 1);
    type_t dy = (frame->ymax - frame->ymin)/(nbgridcells - 1);
    //
    map_mass_grid_GPU(mapfunction,map,cosmo, frame, lens, nhalos,mode_amp,zl, zs, dx, dy, nbgridcells, nbgridcells, 0, 0);
}
//
void map_mass_grid_GPU(map_gpu_function_t mapfunction, type_t *map,const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int mode_amp,  type_t zl, type_t zs, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
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
	type_t *lens_x_gpu, *lens_y_gpu, *b0_gpu, *angle_gpu, *epot_gpu, *rcore_gpu, *rcut_gpu, *anglecos_gpu, *anglesin_gpu, *dlsds_gpu;
	type_t *grid_grad2_a_gpu, *grid_grad2_b_gpu , *grid_grad2_c_gpu, *grid_grad2_d_gpu, *map_gpu;


	if(zl == 0) {
		zl = lens->z[0];
	}

	type_t dl0s = module_cosmodistances_objectObject(zl, zs, *cosmo);
	type_t dos = module_cosmodistances_observerObject(zs, *cosmo);
	type_t dol = module_cosmodistances_observerObject(zl, *cosmo);
	//select_ratio_function(std::string mode, const struct runmode_param* runmode, type_t dls, type_t ds);

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
	cudasafe(cudaMalloc( (void**)&(dlsds_gpu), nhalos*sizeof(type_t)),"Gradient2gpu.cu : Alloc dlsds_gpu: " );
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
	cudasafe(cudaMemcpy(dlsds_gpu, lens->dlsds, nhalos*sizeof(type_t), cudaMemcpyHostToDevice),"Gradient2gpu.cu : Copy dlsds: " );
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
	lens_gpu->dlsds 		= dlsds_gpu;
	//
	cudaMemcpy(lens_kernel, lens_gpu, sizeof(Potential_SOA), cudaMemcpyHostToDevice);
	//
	type_t time = -myseconds();
	//
#if 0
	int shalos = 0;
    while (shalos < nhalos)
    {
            int lens_type = lens->type[shalos];
            type_t z = lens->z[shalos];
            int count     = 1;
            if (shalos + count < nhalos){
            	std::cerr << shalos  << " " <<  count << " " << lens->type[shalos + count] << " " << lens->z[shalos + count] << " " << std::endl;
            	while (lens->type[shalos + count] == lens_type and lens->z[shalos + count] == z ){
            		count++;
            		if(shalos + count >= nhalos)
            			break;
            		//std::cerr << shalos  << " " <<  count << " " << lens->type[shalos + count] << " " << lens->z[shalos + count] << " " << std::endl;
            	}
            	std::cerr << shalos  << " " <<  count << " " << lens->type[shalos + count] << " " << lens->z[shalos + count] << " " << std::endl;
            }
            //if (shalos < nhalos) std::cerr << shalos  << " " <<  count << " " << lens->type[shalos + count] << " " << lens->z[shalos + count] << " " << std::endl;

			shalos += count;
		}
#endif
	//module_potentialDerivatives_totalGradient2_SOA_CPU_GPU(grid_grad2_a_gpu, grid_grad2_b_gpu, grid_grad2_c_gpu, grid_grad2_d_gpu, frame_gpu, lens_kernel, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
	module_potentialDerivatives_Kmap_SOA_CPU_GPU(grid_grad2_a_gpu, grid_grad2_b_gpu, grid_grad2_c_gpu, grid_grad2_d_gpu, frame_gpu, lens_kernel, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
	//
	cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU");
	//std::cerr << "ZS " << zs << " "<< dl0s << " "<< dos <<std::endl;
	mapfunction(map_gpu,grid_grad2_a_gpu, grid_grad2_b_gpu, grid_grad2_c_gpu, grid_grad2_d_gpu,dl0s,dos,dol,cosmo->h,zs,nbgridcells_x,nbgridcells_y,frame);
	cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU");
	cudaDeviceSynchronize();
	//
	cudasafe(cudaMemcpy( map, map_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost )," --- Gradient2gpu.cu : Copy source_a_gpu: " );
	//
	time += myseconds();
	std::cout << "	kernel time = " << time << " s." << std::endl;
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
	cudaFree(dlsds_gpu);
	cudaFree(grid_grad2_a_gpu);
	cudaFree(grid_grad2_b_gpu);
	cudaFree(grid_grad2_c_gpu);
	cudaFree(grid_grad2_d_gpu);
	cudaFree(map_gpu);
}

////General Map calculation
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


//allows for resizing of the source
void map_resizedgrid_GPU(map_gpu_function_t mapfunction, type_t *map,const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int mode_amp, type_t z, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
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

	//Create resized frame information for source amplification and the like
	grid_param resized_frame;
    type_t resized_dx = (frame->xmax - frame->xmin) / 6.;
    type_t resized_dy = (frame->ymax - frame->ymin) / 6.;
	resized_frame.xmax = frame->xmax - resized_dx;
	resized_frame.xmin = frame->xmin + resized_dx;
	resized_frame.ymax = frame->ymax - resized_dy;
	resized_frame.ymin = frame->xmin + resized_dy;


	grid_param *frame_gpu;
	Potential_SOA *lens_gpu,*lens_kernel;
	int *type_gpu;
	type_t *lens_x_gpu, *lens_y_gpu, *b0_gpu, *angle_gpu, *epot_gpu, *rcore_gpu, *rcut_gpu, *anglecos_gpu, *anglesin_gpu;
	type_t *grid_grad2_a_gpu, *grid_grad2_b_gpu , *grid_grad2_c_gpu, *grid_grad2_d_gpu, *map_gpu;

	type_t dl0s = module_cosmodistances_objectObject(lens->z[0], z, *cosmo);
	type_t dos = module_cosmodistances_observerObject(z, *cosmo);
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
	//mapfunction(map_gpu,grid_grad2_a_gpu, grid_grad2_b_gpu, grid_grad2_c_gpu, grid_grad2_d_gpu,dl0s,dos,z,mode_amp,nhalos,nbgridcells_x,nbgridcells_y);
	//amplif_grid_CPU_GPU(map_gpu,grid_grad2_a_gpu, grid_grad2_b_gpu, grid_grad2_c_gpu, grid_grad2_d_gpu,dl0s,z,mode_amp,nhalos,nbgridcells_x,nbgridcells_y);
	//cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU");
	cudaDeviceSynchronize();
	//
	cudasafe(cudaMemcpy( map, map_gpu, (nbgridcells_x)*(nbgridcells_y)*sizeof(type_t), cudaMemcpyDeviceToHost )," --- Gradient2gpu.cu : Copy source_a_gpu: " );
	//
	time += myseconds();
	std::cout << "	kernel time = " << time << " s." << std::endl;
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


////Map functions
//Amplification NR 1
void amplif_1_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame)
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
        cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        amplif_1_grid_GPU<<<grid, threads>>> (map,grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,dl0s,ds,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "amplif_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//
__global__ void amplif_1_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t z,int nbgridcells)
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
        amp = formeli_HPC(A, B, C);
        ampli[index] = 1. / (amp.a * amp.b);
    }
}
//Amplification NR 2
void amplif_2_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame)
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
        cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        amplif_2_grid_GPU<<<grid, threads>>> (map,grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,dl0s,ds,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "amplif_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//
__global__ void amplif_2_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t z,int nbgridcells)
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
        amp = formeli_HPC(A, B, C);
        ampli[index] = 1. / fabs(amp.a * amp.b);
    }
}
//Amplification NR 3
void amplif_3_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame)
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
        cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        amplif_3_grid_GPU<<<grid, threads>>> (map,grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,dl0s,ds,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "amplif_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//
__global__ void amplif_3_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t z,int nbgridcells)
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
        amp = formeli_HPC(A, B, C);
        ampli[index] = -2.5 * log10(fabs(amp.a * amp.b));
    }
}
//Amplification NR 4
void amplif_4_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame)
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
        cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        amplif_4_grid_GPU<<<grid, threads>>> (map,grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,dl0s,ds,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "amplif_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//
__global__ void amplif_4_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t z,int nbgridcells)
{
	//type_t kappa,ga1,ga2,gam,gp;
	////
	type_t dlsds = dl0s/ds;
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    //
    if ((row < nbgridcells) && (col < nbgridcells))
    {

    	int index = row*nbgridcells + col;
        type_t kappa = (grid_grad2_a[index] + grid_grad2_c[index])*dlsds / 2.;
        type_t ga1 = (grid_grad2_a[index] - grid_grad2_c[index])*dlsds / 2.;
        type_t ga2 = grid_grad2_b[index]*dlsds;
        type_t gam = sqrt(ga1 * ga1 + ga2 * ga2);
        type_t gp = gam / (1 - kappa);

        ampli[index] = (1 - kappa) * (1 + gp * gp) / (1 - gp * gp);
    }
}
//Amplification NR 5
void amplif_5_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame)
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
        cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        amplif_5_grid_GPU<<<grid, threads>>> (map,grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,dl0s,ds,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "amplif_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//
__global__ void amplif_5_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t z,int nbgridcells)
{
	////
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
//Amplification NR 6
void amplif_6_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame)
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
        cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        amplif_6_grid_GPU<<<grid, threads>>> (map,grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,dl0s, ds,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "amplif_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//
__global__ void amplif_6_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t z,int nbgridcells)
{
	////
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    //
    if ((row < nbgridcells) && (col < nbgridcells))
    {
    	int index = row*nbgridcells + col;
        type_t ga1 = (grid_grad2_a[index] - grid_grad2_c[index]) / 2.;
        type_t ga2 = grid_grad2_b[index];
        type_t gam = sqrt(ga1 * ga1 + ga2 * ga2);
        ampli[index] = gam;
    }
}

//Mass NR 1
void mass_1_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame)
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
        cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        mass_grid_GPU<<<grid, threads>>> (map,grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,1, dl0s,ds,h,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "amplif_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//Mass NR 2
void mass_2_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame)
{
        int GRID_SIZE_X = (nbgridcells_x + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
        int GRID_SIZE_Y = (nbgridcells_y + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y;
        //
        //printf("grid_size_x = %d, grid_size_y = %d, nbgridcells_x = %d, nbgridcells_y = %d, istart = %d, jstart = %d (split)\n", GRID_SIZE_X, GRID_SIZE_Y, nbgridcells_x, nbgridcells_y, istart, jstart);
        type_t dcrit = cH4piG * h / dl / dl0s * ds;  // in g/cm^2
        //
        dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
        dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);
        //
        //printf("nhalos = %d, size of shared memory = %lf (split)\n", nhalos, (type_t) (8*nhalos + BLOCK_SIZE_X*BLOCK_SIZE_Y)*sizeof(type_t));
        //
        cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        mass_grid_GPU<<<grid, threads>>> (map,grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,dcrit, dl0s,ds,h,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "amplif_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//Mass NR 3
void mass_3_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame)
{
        int GRID_SIZE_X = (nbgridcells_x + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
        int GRID_SIZE_Y = (nbgridcells_y + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y;
        //
    	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells_x - 1);
        type_t dy = (frame->ymax - frame->ymin)/(nbgridcells_y - 1);
        //printf("grid_size_x = %d, grid_size_y = %d, nbgridcells_x = %d, nbgridcells_y = %d, istart = %d, jstart = %d (split)\n", GRID_SIZE_X, GRID_SIZE_Y, nbgridcells_x, nbgridcells_y, istart, jstart);
        type_t conv = MCRIT12 / h * dx * dy * dl / dl0s * ds;  // in  10^12 M_sol/pixel
        //
        dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
        dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);
        //
        //printf("nhalos = %d, size of shared memory = %lf (split)\n", nhalos, (type_t) (8*nhalos + BLOCK_SIZE_X*BLOCK_SIZE_Y)*sizeof(type_t));
        //
        cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        mass_grid_GPU<<<grid, threads>>> (map,grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,conv, dl0s,ds,h,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "amplif_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//Mass NR 4
void mass_4_grid_CPU_GPU(type_t *map,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t dl0s, type_t ds, type_t dl, type_t h, type_t z, int nbgridcells_x, int nbgridcells_y, const struct grid_param *frame)
{
        int GRID_SIZE_X = (nbgridcells_x + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
        int GRID_SIZE_Y = (nbgridcells_y + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y;
        //
        type_t dcritA = cH0_4piG * h / dl / dl0s * ds;  // in 10^12 M_sol/kpc^2
        //
        dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
        dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);
        //
        //printf("nhalos = %d, size of shared memory = %lf (split)\n", nhalos, (type_t) (8*nhalos + BLOCK_SIZE_X*BLOCK_SIZE_Y)*sizeof(type_t));
        //
        cudaMemset(map, 0, nbgridcells_x*nbgridcells_y*sizeof(type_t));
        //
        mass_grid_GPU<<<grid, threads>>> (map,grid_grad2_a, grid_grad2_b,grid_grad2_c, grid_grad2_d,dcritA, dl0s,ds,h,z,nbgridcells_x);
        cudasafe(cudaGetLastError(), "amplif_grid_CPU_GPU");
        //
        cudaDeviceSynchronize();
        printf("GPU kernel done...\n");
}
//
__global__ void mass_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t mult, type_t dls, type_t ds,type_t h, type_t z,int nbgridcells)
{
	////
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    //type_t dlsds = dls / ds ;
    //
    if ((row < nbgridcells) && (col < nbgridcells))
    {
    	int index = row*nbgridcells + col;
        type_t ga1 = (grid_grad2_a[index] + grid_grad2_c[index]) * 0.5 ;//* dlsds ;
        ampli[index] = ga1*mult;
    }
}

__device__ struct  ellipse formeli_HPC(type_t a, type_t b, type_t c)
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

