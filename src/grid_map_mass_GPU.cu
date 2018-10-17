/**
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)
* @date   July 2017
* @version 0,1
*
*/

#include <fstream>
#include "grid_map_mass_GPU.cuh"
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



//GPU mapping function declaration to change when figured out linkage problems
__global__ void mass_grid_GPU(type_t *ampli,type_t *grid_grad2_a,type_t *grid_grad2_b,type_t *grid_grad2_c,type_t *grid_grad2_d, type_t mult, type_t ds, type_t dl, type_t h, type_t z,int nbgridcells);


////Map function selection
map_gpu_function_t select_map_mass_function(const struct runmode_param* runmode){

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
	return 0;
}

////Mass Map calculation, doesnt fit the bloody template...
void map_grid_mass_GPU(map_gpu_function_t mapfunction, type_t *map, const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos ,int nbgridcells,int mode_amp, type_t zl, type_t zs )
{
	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells - 1);
    type_t dy = (frame->ymax - frame->ymin)/(nbgridcells - 1);
    //
    map_grid_mass_GPU(mapfunction,map,cosmo, frame, lens, nhalos,mode_amp,zl, zs, dx, dy, nbgridcells, nbgridcells, 0, 0);
}
//
void map_grid_mass_GPU(map_gpu_function_t mapfunction, type_t *map,const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int mode_amp,  type_t zl, type_t zs, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{

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
    	//if(col == 0 and row == 0)printf(" MASS GPU %f %f\n",grid_grad2_a[index], grid_grad2_c[index]);
        type_t ga1 = (grid_grad2_a[index] + grid_grad2_c[index]) * 0.5 ;//* dlsds ;
        ampli[index] = ga1*mult;

    }
}
