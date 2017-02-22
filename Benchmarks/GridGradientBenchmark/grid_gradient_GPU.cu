#include <fstream>
#include "grid_gradient_GPU.cuh"

#define BLOCK_SIZE 16

void
module_potentialDerivatives_totalGradient_SOA_CPU_GPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_cpu, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos);

void
module_potentialDerivatives_totalGradient_SOA_CPU_GPU_v2(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_cpu, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos);

void calculate_cossin_values(double *theta_cos, double *theta_sin, double *angles, int nhalos ){
	for(int i = 0 ; i < nhalos; i++){
		theta_cos[i]=cos(angles[i]);
		theta_sin[i]=sin(angles[i]);
	}
}

void gradient_grid_GPU_sorted(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos ,int nbgridcells){


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
  Potential_SOA *lens_gpu,*lens_kernel ;
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


  //printf("%p \n", lens_gpu);
  //printf("%p \n", type_gpu);
  //printf("%p \n", lens_gpu->type);
  //fflush(stdout);
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

  cudaMemcpy(lens_kernel, lens_gpu, sizeof(Potential_SOA), cudaMemcpyHostToDevice);


  if (int((nbgridcells) * (nbgridcells)/threadsPerBlock) == 0){
    gradient_grid_kernel<<<1,threadsPerBlock>>>(grid_grad_x_gpu, grid_grad_y_gpu,frame_gpu,nhalos, nbgridcells, lens_kernel);
  }
  else{
    gradient_grid_kernel<<<(nbgridcells) * (nbgridcells)/threadsPerBlock,threadsPerBlock>>>(grid_grad_x_gpu, grid_grad_y_gpu,frame_gpu,nhalos, nbgridcells, lens_kernel);
  }
  cudasafe(cudaMemcpy( grid_grad_x, grid_grad_x_gpu, (nbgridcells) * (nbgridcells) *sizeof(double),cudaMemcpyDeviceToHost ),"Gradientgpu.cu : Copy source_x_gpu: " );
  cudasafe(cudaMemcpy( grid_grad_y, grid_grad_y_gpu, (nbgridcells) * (nbgridcells) *sizeof(double), cudaMemcpyDeviceToHost),"Gradientgpu.cu : Copy source_y_gpu: " );


  //printf("%f %f \n",grid_grad_x[0],grid_grad_y[0]);

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

/*
  for (int i = 0; i < nbgridcells; i++){
    for(int j = 0; j < nbgridcells; j++){
      printf(" %f",grid_grad_x[i*nbgridcells + j]);
    }
    printf("\n");
  }*/

}

void gradient_grid_GPU_multiple(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells){
  //get number of GPU devices
  int nDevices;
  cudaGetDeviceCount(&nDevices);

  // Initialise kernel variables, table for multiple devices
  grid_param *frame_gpu[nDevices];
  Potential_SOA *lens_gpu[nDevices],*lens_kernel[nDevices] ;
  int *type_gpu[nDevices];
  double *lens_x_gpu[nDevices], *lens_y_gpu[nDevices], *b0_gpu[nDevices], *angle_gpu[nDevices], *epot_gpu[nDevices], *rcore_gpu[nDevices], *rcut_gpu[nDevices],*anglecos_gpu[nDevices], *anglesin_gpu[nDevices];
  double *grid_grad_x_gpu[nDevices], *grid_grad_y_gpu[nDevices] ;



  // Initialise multiple device variables
  int Ndevice[nDevices], indexactual[nDevices];
  cudaStream_t stream[nDevices];

  indexactual[0] = 0 ;
  Ndevice[0] = (nbgridcells) * (nbgridcells)/nDevices;

  printf("Using %d Gpu's \n",nDevices );
  //printf("%d %d \n",indexactual[0], Ndevice[0]);

  for (int dev = 1; dev < nDevices; dev++) {

    Ndevice[dev] = (nbgridcells) * (nbgridcells)/nDevices;

    if(indexactual[dev]+Ndevice[dev] > (nbgridcells) * (nbgridcells)){
      Ndevice[dev] = (nbgridcells) * (nbgridcells) - indexactual[dev-1];
    }

    indexactual[dev] = indexactual[dev-1] + Ndevice[dev];
    //printf("%d %d \n",indexactual[dev], Ndevice[dev]);
  }

  for (int dev = 0; dev < nDevices; dev++) {


    cudaSetDevice(dev);

    lens_gpu[dev] = (Potential_SOA *) malloc(sizeof(Potential_SOA));
    lens_gpu[dev]->type = (int *) malloc(sizeof(int));

    // Allocate variables on the GPU
    cudasafe(cudaMalloc( (void**)&(lens_kernel[dev]), sizeof(Potential_SOA)),"Gradientgpu.cu : Alloc Potential_SOA: " );
    cudasafe(cudaMalloc( (void**)&(type_gpu[dev]), nhalos*sizeof(int)),"Gradientgpu.cu : Alloc type_gpu: " );
    cudasafe(cudaMalloc( (void**)&(lens_x_gpu[dev]), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc x_gpu: " );
    cudasafe(cudaMalloc( (void**)&(lens_y_gpu[dev]), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc y_gpu: " );
    cudasafe(cudaMalloc( (void**)&(b0_gpu[dev]), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc b0_gpu: " );
    cudasafe(cudaMalloc( (void**)&(angle_gpu[dev]), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc angle_gpu: " );
    cudasafe(cudaMalloc( (void**)&(epot_gpu[dev]), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc epot_gpu: " );
    cudasafe(cudaMalloc( (void**)&(rcore_gpu[dev]), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc rcore_gpu: " );
    cudasafe(cudaMalloc( (void**)&(rcut_gpu[dev]), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc rcut_gpu: " );
    cudasafe(cudaMalloc( (void**)&(frame_gpu[dev]), sizeof(grid_param)),"Gradientgpu.cu : Alloc frame_gpu: " );
    cudasafe(cudaMalloc( (void**)&(anglecos_gpu[dev]), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc anglecos_gpu: " );
    cudasafe(cudaMalloc( (void**)&(anglesin_gpu[dev]), nhalos*sizeof(double)),"Gradientgpu.cu : Alloc anglesin_gpu: " );
    cudasafe(cudaMalloc( (void**)&(grid_grad_x_gpu[dev]), Ndevice[dev] *sizeof(double)),"Gradientgpu.cu : Alloc source_x_gpu: " );
    cudasafe(cudaMalloc( (void**)&(grid_grad_y_gpu[dev]), Ndevice[dev] *sizeof(double)),"Gradientgpu.cu : Alloc source_y_gpu: " );

    cudaStreamCreate(&stream[dev]);

  }

  for (int dev = 0; dev < nDevices; dev++) {


    cudaSetDevice(dev);

    // Copy values to the GPU
    cudasafe(cudaMemcpyAsync(type_gpu[dev],lens->type , nhalos*sizeof(int),cudaMemcpyHostToDevice,stream[dev] ),"Gradientgpu.cu : Copy type_gpu: " );
    cudasafe(cudaMemcpyAsync(lens_x_gpu[dev],lens->position_x , nhalos*sizeof(double),cudaMemcpyHostToDevice,stream[dev] ),"Gradientgpu.cu : Copy x_gpu: " );
    cudasafe(cudaMemcpyAsync(lens_y_gpu[dev],lens->position_y , nhalos*sizeof(double), cudaMemcpyHostToDevice,stream[dev]),"Gradientgpu.cu : Copy y_gpu: " );
    cudasafe(cudaMemcpyAsync(b0_gpu[dev],lens->b0 , nhalos*sizeof(double), cudaMemcpyHostToDevice,stream[dev]),"Gradientgpu.cu : Copy b0_gpu: " );
    cudasafe(cudaMemcpyAsync(angle_gpu[dev],lens->ellipticity_angle , nhalos*sizeof(double), cudaMemcpyHostToDevice,stream[dev]),"Gradientgpu.cu : Copy angle_gpu: " );
    cudasafe(cudaMemcpyAsync(epot_gpu[dev], lens->ellipticity_potential, nhalos*sizeof(double),cudaMemcpyHostToDevice ,stream[dev]),"Gradientgpu.cu : Copy epot_gpu: " );
    cudasafe(cudaMemcpyAsync(rcore_gpu[dev], lens->rcore, nhalos*sizeof(double),cudaMemcpyHostToDevice ,stream[dev]),"Gradientgpu.cu : Copy rcore_gpu: " );
    cudasafe(cudaMemcpyAsync(rcut_gpu[dev], lens->rcut, nhalos*sizeof(double), cudaMemcpyHostToDevice,stream[dev]),"Gradientgpu.cu : Copy rcut_gpu: " );
    cudasafe(cudaMemcpyAsync(anglecos_gpu[dev], lens->anglecos, nhalos*sizeof(double),cudaMemcpyHostToDevice,stream[dev] ),"Gradientgpu.cu : Copy anglecos: " );
    cudasafe(cudaMemcpyAsync(anglesin_gpu[dev], lens->anglesin, nhalos*sizeof(double), cudaMemcpyHostToDevice,stream[dev]),"Gradientgpu.cu : Copy anglesin: " );
    cudasafe(cudaMemcpyAsync(frame_gpu[dev], frame, sizeof(grid_param), cudaMemcpyHostToDevice,stream[dev]),"Gradientgpu.cu : Copy fame_gpu: " );


    lens_gpu[dev]->type = type_gpu[dev];
    lens_gpu[dev]->position_x = lens_x_gpu[dev];
    lens_gpu[dev]->position_y = lens_y_gpu[dev];
    lens_gpu[dev]->b0 = b0_gpu[dev];
    lens_gpu[dev]->ellipticity_angle = angle_gpu[dev];
    lens_gpu[dev]->ellipticity_potential = epot_gpu[dev];
    lens_gpu[dev]->rcore = rcore_gpu[dev];
    lens_gpu[dev]->rcut = rcut_gpu[dev];
    lens_gpu[dev]->anglecos = anglecos_gpu[dev];
    lens_gpu[dev]->anglesin = anglesin_gpu[dev];

    cudaMemcpyAsync(lens_kernel[dev], lens_gpu[dev], sizeof(Potential_SOA), cudaMemcpyHostToDevice,stream[dev]);
  }

  for (int dev = 0; dev < nDevices; dev++) {


    cudaSetDevice(dev);
    int nBlocks_gpu = 0;
    cudaDeviceProp properties_gpu;
    cudaGetDeviceProperties(&properties_gpu, dev); // Get properties of 0th GPU in use

    if (properties_gpu.maxThreadsDim[0]<threadsPerBlock)
      {
      fprintf(stderr, "ERROR: The GPU has to support at least %u threads per block.\n", threadsPerBlock);
      exit(-1);
      }


    if (int((nbgridcells) * (nbgridcells)/threadsPerBlock) == 0){
    	gradient_grid_kernel_multiple<<<1,threadsPerBlock,0,stream[dev]>>>(grid_grad_x_gpu[dev], grid_grad_y_gpu[dev],frame_gpu[dev],nhalos, nbgridcells, lens_kernel[dev], indexactual[dev],Ndevice[dev]);
    }
    else{
    	gradient_grid_kernel_multiple<<<(nbgridcells) * (nbgridcells)/threadsPerBlock,threadsPerBlock,0,stream[dev]>>>(grid_grad_x_gpu[dev], grid_grad_y_gpu[dev],frame_gpu[dev],nhalos, nbgridcells, lens_kernel[dev], indexactual[dev],Ndevice[dev]);
    }
  }

  for (int dev = 0; dev < nDevices; dev++) {
    cudasafe(cudaMemcpyAsync( grid_grad_x + indexactual[dev], grid_grad_x_gpu[dev], Ndevice[dev] *sizeof(double),cudaMemcpyDeviceToHost ,stream[dev]),"Gradientgpu.cu : Copy source_x_gpu: " );
    cudasafe(cudaMemcpyAsync( grid_grad_y + indexactual[dev], grid_grad_y_gpu[dev], Ndevice[dev] *sizeof(double), cudaMemcpyDeviceToHost,stream[dev]),"Gradientgpu.cu : Copy source_y_gpu: " );
  }

  for (int dev = 0; dev < nDevices; dev++) {


    cudaSetDevice(dev);
      // Free GPU memory
    cudaFree(type_gpu[dev]);
    cudaFree(lens_x_gpu[dev]);
    cudaFree(lens_y_gpu[dev]);
    cudaFree(b0_gpu[dev]);
    cudaFree(angle_gpu[dev]);
    cudaFree(epot_gpu[dev]);
    cudaFree(rcore_gpu[dev]);
    cudaFree(rcut_gpu[dev]);
    cudaFree(anglecos_gpu[dev]);
    cudaFree(anglesin_gpu[dev]);
    cudaFree(grid_grad_x_gpu[dev]);
    cudaFree(grid_grad_y_gpu[dev]);
    cudaStreamDestroy(stream[dev]);

  }

}

#if 1
void gradient_grid_GPU_sub(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells, int indexactual, int Ncells ){


	// GPU Property query
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
	Potential_SOA *lens_gpu,*lens_kernel ;
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


	//printf("%p \n", lens_gpu);
	//printf("%p \n", type_gpu);
	//printf("%p \n", lens_gpu->type);
	//fflush(stdout);
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

	cudaMemcpy(lens_kernel, lens_gpu, sizeof(Potential_SOA), cudaMemcpyHostToDevice);
#if 0
	int BLOCK_SIZE = 16; // number of threads
	int GRID_SIZE = (nbgridcells + BLOCK_SIZE - 1)/BLOCK_SIZE; // number of blocks
	//
	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
	dim3 dimGrid(GRID_SIZE, GRID_SIZE);
	//
	if (int((nbgridcells) * (nbgridcells)/threadsPerBlock) == 0)
	{
		gradient_grid_kernel<<<1,threadsPerBlock>>>(grid_grad_x_gpu, grid_grad_y_gpu,frame_gpu,nhalos, nbgridcells, lens_kernel);
	}
	else
	{
		//gradient_grid_kernel<<<(nbgridcells) * (nbgridcells)/threadsPerBlock,threadsPerBlock>>>(grid_grad_x_gpu, grid_grad_y_gpu,frame_gpu,nhalos, nbgridcells, lens_kernel);
		//gradient_grid_kernel_v2<<<dimGrid, dimBlock>>>(grid_grad_x_gpu, grid_grad_y_gpu, frame_gpu, nhalos, nbgridcells, lens_kernel);
		//gradient_grid_kernel_v2<<<dimGrid, dimBlock>>>(grid_grad_x_gpu, grid_grad_y_gpu, frame_gpu, nhalos, nbgridcells, lens_kernel);


	if (int((Ncells)/threadsPerBlock) == 0){
	gradient_grid_kernel_multiple<<<1,threadsPerBlock>>>(grid_grad_x_gpu, grid_grad_y_gpu,frame_gpu,nhalos, nbgridcells, lens_kernel, indexactual, Ncells);
	}
	else{
	gradient_grid_kernel_multiple<<<(Ncells)/threadsPerBlock,threadsPerBlock>>>(grid_grad_x_gpu, grid_grad_y_gpu,frame_gpu,nhalos, nbgridcells, lens_kernel, indexactual, Ncells);
	}
#endif
		module_potentialDerivatives_totalGradient_SOA_CPU_GPU(grid_grad_x_gpu, grid_grad_y_gpu, frame_gpu, lens, lens_kernel, nbgridcells, nhalos);
	cudasafe(cudaMemcpy( grid_grad_x, grid_grad_x_gpu, (nbgridcells) * (nbgridcells) *sizeof(double),cudaMemcpyDeviceToHost ),"Gradientgpu.cu : Copy source_x_gpu: " );
	cudasafe(cudaMemcpy( grid_grad_y, grid_grad_y_gpu, (nbgridcells) * (nbgridcells) *sizeof(double), cudaMemcpyDeviceToHost),"Gradientgpu.cu : Copy source_y_gpu: " );


	//printf("%f %f \n",grid_grad_x[0],grid_grad_y[0]);

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

	/*
	   for (int i = 0; i < nbgridcells; i++){
	   for(int j = 0; j < nbgridcells; j++){
	   printf(" %f",grid_grad_x[i*nbgridcells + j]);
	   }
	   printf("\n");
	   }*/

}


__global__ void gradient_grid_kernel(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcells, const struct Potential_SOA *lens) {

	//*grad_x = *grad_y = 0.;

	int bid=blockIdx.x; // index of the block (and of the set of images)
	int tid=threadIdx.x; // index of the thread within the block

	double dx,dy;        //pixelsize
	int grid_dim, index;
	struct point image_point, Grad;
	dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
	dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
	grid_dim = (nbgridcells);

	index = bid * threadsPerBlock + tid ;

	while(index < grid_dim*grid_dim){

		grid_grad_x[index] = 0.;
		grid_grad_y[index] = 0.;

		image_point.x = frame->xmin + (index/grid_dim)*dx;
		image_point.y = frame->ymin + (index % grid_dim)*dy;

		Grad = module_potentialDerivatives_totalGradient_SOA_GPU(&image_point, lens, Nlens);

		grid_grad_x[index] = Grad.x;
		grid_grad_y[index] = Grad.y;

		bid += gridDim.x;
		index = bid * threadsPerBlock + tid;
	}
}

__global__ void gradient_grid_kernel_v2(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcells, const struct Potential_SOA *lens) {

	//*grad_x = *grad_y = 0.;

	int bid = blockIdx.x; // index of the block (and of the set of images)
	int tid = threadIdx.x; // index of the thread within the block

	double dx,dy;        //pixelsize
	int grid_dim, index;
	struct point image_point, Grad;
	//
	dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
	dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
	//
	grid_dim = (nbgridcells);
	//
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	//
	index = col*nbgridcells + row ;
	//
	//while(index < grid_dim*grid_dim){

	//grid_grad_x[index] = 0.;
	//grid_grad_y[index] = 0.;

	image_point.x = frame->xmin + col*dx;
	image_point.y = frame->ymin + row*dy;

	Grad = module_potentialDerivatives_totalGradient_SOA_GPU(&image_point, lens, Nlens);

	grid_grad_x[index] = Grad.x;
	grid_grad_y[index] = Grad.y;

	bid += gridDim.x;
	index = bid * threadsPerBlock + tid;
	//}
}

/*
void gradient_grid_general_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcells, const struct Potential_SOA *lens){
  int bid=0; // index of the block (and of the set of images)
  int tid=0; // index of the thread within the block

  double dx,dy,x_pos,y_pos;        //pixelsize
  int grid_dim, index;
  point Grad, image_point, true_coord_rotation;
  double      R;
  dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
  dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
  grid_dim = (nbgridcells);

  index = bid ;

  while(index < grid_dim*grid_dim){

    grid_grad_x[index] = 0.;
    grid_grad_y[index] = 0.;

    image_point.x = frame->xmin + (index/grid_dim)*dx;
    image_point.y = frame->ymin + (index % grid_dim)*dy;

    Grad = module_potentialDerivatives_totalGradient_SOA(&image_point, lens, Nlens);

    grid_grad_x[index] = Grad.x;
    grid_grad_y[index] = Grad.y;

    bid += 1;
    index = bid * 1 + tid;
  }
}

*/

__global__
void 
inline
module_potentialDerivatives_totalGradient_8_SOA_GPU(double *grid_grad_x, double *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int shalos, int nhalos)
{
        //asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        //
        struct point grad, clumpgrad, image_point;
        grad.x = 0;
        grad.y = 0;
	//
	int row = blockIdx.y * blockDim.y + threadIdx.y;
        int col = blockIdx.x * blockDim.x + threadIdx.x;
	//
	if ((row < nbgridcells) && (col < nbgridcells))
	{
		//
		int index = col*nbgridcells + row;
		//
		//grid_grad_x[index] = 0.;
		//grid_grad_y[index] = 0.;
		//
		double dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
		double dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
		//
#if 0
		__shared__ double img_pt[2];
		if ((row == 0) && (col == 0))
		{
			img_pt[0] = frame->xmin + col*dx;
			img_pt[1] = frame->ymin + row*dy;
		}
		__syncthreads();
#else
		image_point.x = frame->xmin + col*dx;
		image_point.y = frame->ymin + row*dy;	
#endif
		//
		//
		for(int i = shalos; i < shalos + nhalos; i++)
		{
			//IACA_START;
			//
			struct point true_coord, true_coord_rot; //, result;
			//double       R, angular_deviation;
			complex      zis;
			//
			//result.x = result.y = 0.;
			//
#if 0
			true_coord.x = img_pt[0] - __ldg(&lens->position_x[i]);
			true_coord.y = img_pt[1] - __ldg(&lens->position_y[i]);
#else
			true_coord.x = image_point.x - __ldg(&lens->position_x[i]);
			true_coord.y = image_point.y - __ldg(&lens->position_y[i]);
#endif
			double cosi = __ldg(&lens->anglecos[i]);
			double sinu = __ldg(&lens->anglesin[i]);
			// positionning at the potential center
			// Change the origin of the coordinate system to the center of the clump
			double x = true_coord.x*cosi + true_coord.y*sinu;
			double y = true_coord.y*cosi - true_coord.x*sinu;
			//
			double eps = __ldg(&lens->ellipticity_potential[i]);
			//
			double sqe  = sqrt(eps);
			//
			double rem2 = x*x/((1. + eps)*(1. + eps)) + y*y/((1. - eps)*(1. - eps));
			//
			complex zci;
			complex znum, zden, zres;
			double norm;
			//
			zci.im  = -0.5*(1. - eps*eps)/sqe;
			//
			double rc  = __ldg(&lens->rcore[i]);
			double cx1  = (1. - eps)/(1. + eps);
			znum.re = cx1*x;
			znum.im = 2.*sqe*sqrt(rc*rc + rem2) - y/cx1;
			//
			zden.re = x;
			zden.im = 2.*rc*sqe - y;
			norm    = (zden.re*zden.re + zden.im*zden.im);     // zis = znum/zden
			//
			zis.re  = (znum.re*zden.re + znum.im*zden.im)/norm;
			zis.im  = (znum.im*zden.re - znum.re*zden.im)/norm;
			//
			norm    = zis.re;
			//
			zis.re  = log(sqrt(norm*norm + zis.im*zis.im));  // ln(zis) = ln(|zis|)+i.Arg(zis)
			zis.im  = atan2(zis.im, norm);
			//
			zres.re = zci.im*zis.im;   // Re( zci*ln(zis) )
			zres.im = zci.im*zis.re;   // Im( zci*ln(zis) )
			//
			double b0  = __ldg(&lens->b0[i]);
			grad.x += b0*(zres.re*cosi - zres.im*sinu);
			grad.y += b0*(zres.im*cosi + zres.re*sinu);
		}
		//IACA_END;
		//
		grid_grad_x[index] += grad.x;
		grid_grad_y[index] += grad.y;
	}
}

__global__
void
module_potentialDerivatives_totalGradient_8_SOA_GPU_v2(double *grid_grad_x, double *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame,
 int nbgridcells, int i, int nhalos)
{
        //asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        //
        struct point grad, clumpgrad, image_point;
        grad.x = 0;
        grad.y = 0;
        //
        int row = blockIdx.y * blockDim.y + threadIdx.y;
        int col = blockIdx.x * blockDim.x + threadIdx.x;
        //
        if ((row < nbgridcells) && (col < nbgridcells))
        {
                //
                int index = col*nbgridcells + row;
                //
                //grid_grad_x[index] = 0.;
                //grid_grad_y[index] = 0.;
                //
                double dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
                double dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
                //
#if 0
                __shared__ double img_pt[2];
                if ((row == 0) && (col == 0))
                {
                        img_pt[0] = frame->xmin + col*dx;
                        img_pt[1] = frame->ymin + row*dy;
                }
                __syncthreads();
#else
                image_point.x = frame->xmin + col*dx;
                image_point.y = frame->ymin + row*dy;
#endif
                //
                //
                //for(int i = shalos; i < shalos + nhalos; i++)
                //{
                        //IACA_START;
                        //
                        struct point true_coord, true_coord_rot; //, result;
                        //double       R, angular_deviation;
                        complex      zis;
                        //
                        //result.x = result.y = 0.;
                        //
#if 0
                        true_coord.x = img_pt[0] - __ldg(&lens->position_x[i]);
                        true_coord.y = img_pt[1] - __ldg(&lens->position_y[i]);
#else
                        true_coord.x = image_point.x - __ldg(&lens->position_x[i]);
                        true_coord.y = image_point.y - __ldg(&lens->position_y[i]);
#endif
                        double cosi = __ldg(&lens->anglecos[i]);
                        double sinu = __ldg(&lens->anglesin[i]);
                        // positionning at the potential center
                        // Change the origin of the coordinate system to the center of the clump
                        double x = true_coord.x*cosi + true_coord.y*sinu;
                        double y = true_coord.y*cosi - true_coord.x*sinu;
                        //
                        double eps = __ldg(&lens->ellipticity_potential[i]);
                        //
                        double sqe  = sqrt(eps);
                        //
                        double rem2 = x*x/((1. + eps)*(1. + eps)) + y*y/((1. - eps)*(1. - eps));
                        //
                        complex zci;
                        complex znum, zden, zres;
                        double norm;
                        //
                        zci.im  = -0.5*(1. - eps*eps)/sqe;
                        //
                        double rc  = __ldg(&lens->rcore[i]);
                        double cx1  = (1. - eps)/(1. + eps);
                        znum.re = cx1*x;
                        znum.im = 2.*sqe*sqrt(rc*rc + rem2) - y/cx1;
                        //
                        zden.re = x;
                        zden.im = 2.*rc*sqe - y;
                        norm    = (zden.re*zden.re + zden.im*zden.im);     // zis = znum/zden
                        //
                        zis.re  = (znum.re*zden.re + znum.im*zden.im)/norm;
                        zis.im  = (znum.im*zden.re - znum.re*zden.im)/norm;
                        //
                        //
                        double b0  = __ldg(&lens->b0[i]);
                        grad.x += b0*(zres.re*cosi - zres.im*sinu);
                        grad.y += b0*(zres.im*cosi + zres.re*sinu);
                //}
                //IACA_END;
                //
                grid_grad_x[index] += grad.x;
                grid_grad_y[index] += grad.y;
        }
}





/*
   typedef struct point (*halo_func_GPU_t) (const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);

   __constant__ halo_func_GPU_t halo_func_GPU[100] =
   {
   0, 0, 0, 0, 0, module_potentialDerivatives_totalGradient_5_SOA_GPU, 0, 0, module_potentialDerivatives_totalGradient_8_SOA_GPU,  0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0,  module_potentialDerivatives_totalGradient_81_SOA_GPU, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0
   };
 */

	void
module_potentialDerivatives_totalGradient_SOA_CPU_GPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_cpu, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos)
{
	struct point grad, clumpgrad;
	//
	grad.x = clumpgrad.x = 0;
	grad.y = clumpgrad.y = 0;
	int shalos = 0;
	int GRID_SIZE = (nbgridcells + BLOCK_SIZE - 1)/BLOCK_SIZE; // number of blocks
	//
	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
	dim3 dimGrid(GRID_SIZE, GRID_SIZE);	
	int count = nhalos;	
	cudaMemset(grid_grad_x, 0, nbgridcells*nbgridcells*sizeof(double));
	cudaMemset(grid_grad_y, 0, nbgridcells*nbgridcells*sizeof(double));
	module_potentialDerivatives_totalGradient_8_SOA_GPU<<<dimGrid, dimBlock>>> (grid_grad_x, grid_grad_y, lens_gpu, frame, nbgridcells, shalos, count);
	//grid.x += clumpgrad.x;
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


        void
module_potentialDerivatives_totalGradient_SOA_CPU_GPU_v2(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_cpu, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos)
{
        struct point grad, clumpgrad;
        //
        grad.x = clumpgrad.x = 0;
        grad.y = clumpgrad.y = 0;
        int shalos = 0;
        int GRID_SIZE = (nbgridcells + BLOCK_SIZE - 1)/BLOCK_SIZE; // number of blocks
        //
        dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
        dim3 dimGrid(GRID_SIZE, GRID_SIZE);
        int count = nhalos;
	//
	cudaMemset(grid_grad_x, 0, nbgridcells*nbgridcells*sizeof(double));
	cudaMemset(grid_grad_y, 0, nbgridcells*nbgridcells*sizeof(double));
	//
	for (int ii = 0; ii < nhalos; ++ii)
		module_potentialDerivatives_totalGradient_8_SOA_GPU<<<dimGrid, dimBlock>>> (grid_grad_x, grid_grad_y, lens_gpu, frame, nbgridcells, ii, 1);
	
        //grid.x += clumpgrad.x;
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

