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
extern __shared__ double shared[];
#else
#define SHARED 
#endif

#define Nx 1
#define Ny 0


#define cudasafe 

extern "C" 
{
	double myseconds();
}

__global__ void module_potentialDerivatives_totalGradient_SOA_GPU(double *grid_grad_x, double *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int nhalos);

void gradient_grid_GPU_sorted(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos ,int nbgridcells)
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
	double *lens_x_gpu, *lens_y_gpu, *b0_gpu, *angle_gpu, *epot_gpu, *rcore_gpu, *rcut_gpu, *anglecos_gpu, *anglesin_gpu;
	double *grid_grad_x_gpu, *grid_grad_y_gpu ;

	lens_gpu = (Potential_SOA *) malloc(sizeof(Potential_SOA));
	lens_gpu->type = (int *) malloc(sizeof(int));

	// Allocate variables on the GPU
	cudasafe(cudaMalloc( (void**)&(lens_kernel), sizeof(Potential_SOA)),"Gradientgpu.cu : Alloc Potential_SOA: " );
	cudasafe(cudaMalloc( (void**)&(type_gpu), nhalos*sizeof(int)),"Gradientgpu.cu : Alloc type_gpu: " );
	cudasafe(cudaMalloc( (void**)&(lens_x_gpu), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc x_gpu: " );
	cudasafe(cudaMalloc( (void**)&(lens_y_gpu), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc y_gpu: " );
	cudasafe(cudaMalloc( (void**)&(b0_gpu), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc b0_gpu: " );
	cudasafe(cudaMalloc( (void**)&(angle_gpu), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc angle_gpu: " );
	cudasafe(cudaMalloc( (void**)&(epot_gpu), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc epot_gpu: " );
	cudasafe(cudaMalloc( (void**)&(rcore_gpu), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc rcore_gpu: " );
	cudasafe(cudaMalloc( (void**)&(rcut_gpu), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc rcut_gpu: " );
	cudasafe(cudaMalloc( (void**)&(anglecos_gpu), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc anglecos_gpu: " );
	cudasafe(cudaMalloc( (void**)&(anglesin_gpu), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc anglesin_gpu: " );
	cudasafe(cudaMalloc( (void**)&(frame_gpu), sizeof(grid_param)),"Gradientgpu.cu : Alloc frame_gpu: " );
	cudasafe(cudaMalloc( (void**)&(grid_grad_x_gpu), (nbgridcells) * (nbgridcells) *sizeof(double)),"Gradientgpu.cu : Alloc source_x_gpu: " );
	cudasafe(cudaMalloc( (void**)&(grid_grad_y_gpu), (nbgridcells) * (nbgridcells) *sizeof(double)),"Gradientgpu.cu : Alloc source_y_gpu: " );

	// Copy values to the GPU
	cudasafe(cudaMemcpy(type_gpu,lens->type , nhalos*sizeof(int),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy type_gpu: " );
	cudasafe(cudaMemcpy(lens_x_gpu,lens->position_x , nhalos*sizeof(double),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy x_gpu: " );
	cudasafe(cudaMemcpy(lens_y_gpu,lens->position_y , nhalos*sizeof(double), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy y_gpu: " );
	cudasafe(cudaMemcpy(b0_gpu,lens->b0 , nhalos*sizeof(double), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy b0_gpu: " );
	cudasafe(cudaMemcpy(angle_gpu,lens->ellipticity_angle , nhalos*sizeof(double), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy angle_gpu: " );
	cudasafe(cudaMemcpy(epot_gpu, lens->ellipticity_potential, nhalos*sizeof(double),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy epot_gpu: " );
	cudasafe(cudaMemcpy(rcore_gpu, lens->rcore, nhalos*sizeof(double),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy rcore_gpu: " );
	cudasafe(cudaMemcpy(rcut_gpu, lens->rcut, nhalos*sizeof(double), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy rcut_gpu: " );
	cudasafe(cudaMemcpy(anglecos_gpu, lens->anglecos, nhalos*sizeof(double),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy anglecos: " );
	cudasafe(cudaMemcpy(anglesin_gpu, lens->anglesin, nhalos*sizeof(double), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy anglesin: " );
	cudasafe(cudaMemcpy(frame_gpu, frame, sizeof(grid_param), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy fame_gpu: " );
	//
	lens_gpu->type = type_gpu;
	lens_gpu->position_x = lens_x_gpu;
	lens_gpu->position_y = lens_y_gpu;
	lens_gpu->b0 = b0_gpu;
	lens_gpu->ellipticity_angle = angle_gpu;
	lens_gpu->ellipticity_potential = epot_gpu;
	lens_gpu->rcore = rcore_gpu;
	lens_gpu->rcut = rcut_gpu;
	lens_gpu->anglecos = anglecos_gpu;
	lens_gpu->anglesin = anglesin_gpu;
	//
	cudaMemcpy(lens_kernel, lens_gpu, sizeof(Potential_SOA), cudaMemcpyHostToDevice);
	//
	double time = -myseconds();
	module_potentialDerivatives_totalGradient_SOA_CPU_GPU(grid_grad_x_gpu, grid_grad_y_gpu, frame_gpu, lens, lens_kernel, nbgridcells, nhalos);
	//
	//cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU");
	cudaDeviceSynchronize();
	time += myseconds();
	std::cout << "	kernel time = " << time << " s." << std::endl;
	//
	cudasafe(cudaMemcpy( grid_grad_x, grid_grad_x_gpu, (nbgridcells)*(nbgridcells)*sizeof(double), cudaMemcpyDeviceToHost )," --- Gradientgpu.cu : Copy source_x_gpu: " );
	cudasafe(cudaMemcpy( grid_grad_y, grid_grad_y_gpu, (nbgridcells)*(nbgridcells)*sizeof(double), cudaMemcpyDeviceToHost)," --- Gradientgpu.cu : Copy source_y_gpu: " );
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
	cudaFree(grid_grad_x_gpu);
	cudaFree(grid_grad_y_gpu);
}


__global__
void
module_potentialDerivatives_totalGradient_8_SOA_GPU_loc(double *grid_grad_x, double *grid_grad_y, const struct Potential_SOA
 *lens, const struct grid_param *frame, int nbgridcells, int shalos, int nhalos)
{
        //
        //asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        //
        double grad_x, grad_y;
        double clumpgrad_x, clumpgrad_y;
        double image_point_x, image_point_y;
        //
        SHARED double cosi      [200];
        SHARED double sinu      [200];
        SHARED double rci       [200];
        SHARED double b0        [200];
        SHARED double epsi      [200];
        SHARED double position_x[200];
        SHARED double position_y[200];
        SHARED double rsqe      [200];
        //SHARED double sgrad_x   [(BLOCK_SIZE_X + 1)*BLOCK_SIZE_Y];
        //SHARED double sgrad_y   [(BLOCK_SIZE_X + 1)*BLOCK_SIZE_Y];

        //SHARED double sonepeps  [200];
        //SHARED double sonemeps  [200];
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
        double dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
        double dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
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
                        double rc      = rci[i];
                        double eps     = epsi[i];
                        double onemeps = 1 - eps;
                        double onepeps = 1 + eps;
                        //
                        double sqe = rsqe[i];
                        double cx1 = onemeps/onepeps;
                        //
                        //
                        //double zci_im = 1;
                        double zci_im  = -0.5*(1. - eps*eps)/sqe;
                        double inv_onepeps = 1./(onepeps*onepeps);
                        double inv_onemeps = 1./(onemeps*onemeps);
#endif
                        //
                        {
                                //KERNEL_8;
                                double grad_x = grad_y = 0.;
                                double true_coord_y = image_point_y - position_y[i];
                                double true_coord_x = image_point_x - position_x[i];
                                KERNEL_8_reg(0);
                                grid_grad_x[iglob +  0] += grad_x;
                                grid_grad_y[iglob +  0] += grad_y;
                        }
                        /*
                        {
                                //KERNEL_8;
                                double grad_x = grad_y = 0.;
                                double true_coord_y = image_point_y - position_y[i];
                                double true_coord_x = image_point_x - position_x[i] + BLOCK_SIZE_X/2*dx;
                                KERNEL_8_reg(0);
                                grid_grad_x[iglob + BLOCK_SIZE_X/2] += grad_x;
                                grid_grad_y[iglob + BLOCK_SIZE_X/2] += grad_y;
                        }
                        */
                //
                }
        }
}


void
module_potentialDerivatives_totalGradient_SOA_CPU_GPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_cpu, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos)
{
	struct point grad, clumpgrad;
	//
	grad.x = clumpgrad.x = 0.;
	grad.y = clumpgrad.y = 0.;
	//
	int shalos = 0;
	//
	int GRID_SIZE_X = (nbgridcells + BLOCK_SIZE_X - 1)/BLOCK_SIZE_X; // number of blocks
	int GRID_SIZE_Y = (nbgridcells + BLOCK_SIZE_Y - 1)/BLOCK_SIZE_Y; 
	printf("grid_size_x = %d, grid_size_y = %d\n", GRID_SIZE_X, GRID_SIZE_Y);
	//
	double* timer = (double *) malloc((int) nbgridcells*nbgridcells*sizeof(double)); 
	double* dtimer;
	//cudasafe(cudaMalloc( (void**)&(dtimer), (int) nbgridcells*nbgridcells*sizeof(double)),"Gradientgpu.cu : totalGradient_SOA_CPU_GPU: " );
	//
	{
		dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y/1);
		dim3 grid   (GRID_SIZE_X , GRID_SIZE_Y);	
		//
		int count = nhalos;	
		printf("nhalos = %d, size of shared memory = %lf\n", nhalos, (double) (8*nhalos + BLOCK_SIZE_X*BLOCK_SIZE_Y)*sizeof(double));
		//
		cudaMemset(grid_grad_x, 0, nbgridcells*nbgridcells*sizeof(double));
		cudaMemset(grid_grad_y, 0, nbgridcells*nbgridcells*sizeof(double));
		//
		//module_potentialDerivatives_totalGradient_8_SOA_GPU_cur<<<grid, threads>>> (grid_grad_x, grid_grad_y, lens_gpu, frame, nbgridcells, shalos, count);
		module_potentialDerivatives_totalGradient_SOA_GPU<<<grid, threads>>> (grid_grad_x, grid_grad_y, lens_gpu, frame, nbgridcells, nhalos);
		//module_potentialDerivatives_totalGradient_8_SOA_GPU_SM2<<<grid, threads>>> (grid_grad_x, grid_grad_y, lens_gpu, frame, nbgridcells, shalos, count);
		//module_potentialDerivatives_totalGradient_8_SOA_GPU_SM4<<<dimGrid, dimBlock, (8*nhalos + BLOCK_SIZE_X*nbgridcells/BLOCK_SIZE_Y)*sizeof(double)>>> (grid_grad_x, grid_grad_y, lens_gpu, frame, nbgridcells, shalos, count/*, dtimer*/);
		cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU_8_SOA_GPU_SM4");
		
	}
	cudaDeviceSynchronize();
	printf("GPU kernel done... ");fflush(stdout);
	//cudasafe(cudaMemcpy( timer, dtimer, nbgridcells*nbgridcells*sizeof(double), cudaMemcpyDeviceToHost ),"Gradientgpu.cu : dtimer memcpy " );
	
	//
	/*
	   {
	   dim3 dimBlock(BLOCK_SIZE/1, BLOCK_SIZE/2);
	   dim3 dimGrid(GRID_SIZE, GRID_SIZE);
	//
	int count = nhalos;
	cudaMemset(grid_grad_x, 0, nbgridcells*nbgridcells*sizeof(double));
	cudaMemset(grid_grad_y, 0, nbgridcells*nbgridcells*sizeof(double));
	//
	module_potentialDerivatives_totalGradient_8_SOA_GPU_SM<<<dimGrid, dimBlock>>> (grid_grad_x, grid_grad_y, lens_gpu, frame, nbgridcells, shalos, count);
	}
	 */
	//
	cudasafe(cudaGetLastError(), "module_potentialDerivative_totalGradient_SOA_CPU_GPU");
	//
	printf("\n");
	//module_potentialDerivatives_totalGradient_8_SOA_GPU_SM<<<dimGrid, dimBlock>>> (grid_grad_x, grid_grad_y, lens_gpu, frame, nbgridcells, shalos, count);
	//grad.x += clumpgrad.x;
	//grad.y += clumpgrad.y;

	//
	//	
	/*
	   while (shalos < nhalos)
	   {

	   int lens_type = lens_cpu->type[shalos];
	   int count     = 1;
	   while (lens_cpu->type[shalos + count] == lens_type) count++;
	//std::cerr << "type = " << lens_type << " " << count << " " << shalos << std::endl;
	//printf ("%d %d %d \n",lens_type,count,shalos);
	//
	clumpgrad = (*halo_func_GPU[lens_type]<<<dimGrid, dimBlock>>> )(lens_gpu, shalos, count);
	//
	grad.x += clumpgrad.x;
	grad.y += clumpgrad.y;
	shalos += count;
	}

	return(grad);
	 */
}


