#include <fstream>
#include "grid_gradient_GPU.cuh"



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
  double *lens_x_gpu, *lens_y_gpu, *b0_gpu, *angle_gpu, *epot_gpu, *rcore_gpu, *rcut_gpu;
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


void gradient_grid_pinned(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int *Nlens,int nbgridcells){


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
  double *lens_x_gpu, *lens_y_gpu, *b0_gpu, *angle_gpu, *epot_gpu, *rcore_gpu, *rcut_gpu;
  double *grid_grad_x_gpu, *grid_grad_y_gpu ;

  //SIS//

  // Allocate variables on the GPU
  cudasafe(cudaMalloc( (void**)&(lens_x_gpu), Nlens[0]*sizeof(double)),"Gradientgpu.cu : Alloc x_gpu: " );
  cudasafe(cudaMalloc( (void**)&(lens_y_gpu), Nlens[0]*sizeof(double)),"Gradientgpu.cu : Alloc y_gpu: " );
  cudasafe(cudaMalloc( (void**)&(b0_gpu), Nlens[0]*sizeof(double)),"Gradientgpu.cu : Alloc b0_gpu: " );
  cudasafe(cudaMalloc( (void**)&(angle_gpu), Nlens[0]*sizeof(double)),"Gradientgpu.cu : Alloc angle_gpu: " );
  cudasafe(cudaMalloc( (void**)&(epot_gpu), Nlens[0]*sizeof(double)),"Gradientgpu.cu : Alloc epot_gpu: " );
  cudasafe(cudaMalloc( (void**)&(rcore_gpu), Nlens[0]*sizeof(double)),"Gradientgpu.cu : Alloc rcore_gpu: " );
  cudasafe(cudaMalloc( (void**)&(rcut_gpu), Nlens[0]*sizeof(double)),"Gradientgpu.cu : Alloc rcut_gpu: " );
  cudasafe(cudaMalloc( (void**)&(frame_gpu), sizeof(grid_param)),"Gradientgpu.cu : Alloc frame_gpu: " );
  cudasafe(cudaMalloc( (void**)&(grid_grad_x_gpu), (nbgridcells) * (nbgridcells) *sizeof(double)),"Gradientgpu.cu : Alloc source_x_gpu: " );
  cudasafe(cudaMalloc( (void**)&(grid_grad_y_gpu), (nbgridcells) * (nbgridcells) *sizeof(double)),"Gradientgpu.cu : Alloc source_y_gpu: " );

  // Copy values to the GPU
  cudasafe(cudaMemcpy(lens_x_gpu,lens[0].position_x , Nlens[0]*sizeof(double),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy x_gpu: " );
  cudasafe(cudaMemcpy(lens_y_gpu,lens[0].position_y , Nlens[0]*sizeof(double), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy y_gpu: " );
  cudasafe(cudaMemcpy(b0_gpu,lens[0].b0 , Nlens[0]*sizeof(double), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy b0_gpu: " );
  cudasafe(cudaMemcpy(angle_gpu,lens[0].ellipticity_angle , Nlens[0]*sizeof(double), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy angle_gpu: " );
  cudasafe(cudaMemcpy(epot_gpu, lens[0].ellipticity_potential, Nlens[0]*sizeof(double),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy epot_gpu: " );
  cudasafe(cudaMemcpy(rcore_gpu, lens[0].rcore, Nlens[0]*sizeof(double),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy rcore_gpu: " );
  cudasafe(cudaMemcpy(rcut_gpu, lens[0].rcut, Nlens[0]*sizeof(double), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy rcut_gpu: " );
  cudasafe(cudaMemcpy(frame_gpu, frame, sizeof(grid_param), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy fame_gpu: " );

  //printf( "%d %d",nBlocks_gpu,(nbgridcells) * (nbgridcells)/threadsPerBlock );

  if (int((nbgridcells) * (nbgridcells)/threadsPerBlock) == 0){
    gradient_grid_sis_GPU<<<1,threadsPerBlock>>>(grid_grad_x_gpu, grid_grad_y_gpu,frame_gpu,Nlens[0], nbgridcells, lens_x_gpu, lens_y_gpu, b0_gpu, angle_gpu, epot_gpu, rcore_gpu, rcut_gpu);
  }
  else{
  gradient_grid_sis_GPU<<<(nbgridcells) * (nbgridcells)/threadsPerBlock,threadsPerBlock>>>(grid_grad_x_gpu, grid_grad_y_gpu,frame_gpu,Nlens[0], nbgridcells, lens_x_gpu, lens_y_gpu, b0_gpu, angle_gpu, epot_gpu, rcore_gpu, rcut_gpu);
  }
  cudasafe(cudaMemcpy( grid_grad_x, grid_grad_x_gpu, (nbgridcells) * (nbgridcells) *sizeof(double),cudaMemcpyDeviceToHost ),"Gradientgpu.cu : Copy source_x_gpu: " );
  cudasafe(cudaMemcpy( grid_grad_y, grid_grad_y_gpu, (nbgridcells) * (nbgridcells) *sizeof(double), cudaMemcpyDeviceToHost),"Gradientgpu.cu : Copy source_y_gpu: " );
  //printf ("%f %f \n", grid_grad_x[0],grid_grad_y[0]);

    // Free GPU memory
  cudaFree(lens_x_gpu);
  cudaFree(lens_y_gpu);
  cudaFree(b0_gpu);
  cudaFree(angle_gpu);
  cudaFree(epot_gpu);
  cudaFree(rcore_gpu);
  cudaFree(rcut_gpu);
  cudaFree(grid_grad_x_gpu);
  cudaFree(grid_grad_y_gpu);


  //PIEMD8//

  // Allocate variables on the GPU
  cudasafe(cudaMalloc( (void**)&(lens_x_gpu), Nlens[1]*sizeof(double)),"Gradientgpu.cu : Alloc x_gpu: " );
  cudasafe(cudaMalloc( (void**)&(lens_y_gpu), Nlens[1]*sizeof(double)),"Gradientgpu.cu : Alloc y_gpu: " );
  cudasafe(cudaMalloc( (void**)&(b0_gpu), Nlens[1]*sizeof(double)),"Gradientgpu.cu : Alloc b0_gpu: " );
  cudasafe(cudaMalloc( (void**)&(angle_gpu), Nlens[1]*sizeof(double)),"Gradientgpu.cu : Alloc angle_gpu: " );
  cudasafe(cudaMalloc( (void**)&(epot_gpu), Nlens[1]*sizeof(double)),"Gradientgpu.cu : Alloc epot_gpu: " );
  cudasafe(cudaMalloc( (void**)&(rcore_gpu), Nlens[1]*sizeof(double)),"Gradientgpu.cu : Alloc rcore_gpu: " );
  cudasafe(cudaMalloc( (void**)&(rcut_gpu), Nlens[1]*sizeof(double)),"Gradientgpu.cu : Alloc rcut_gpu: " );
  cudasafe(cudaMalloc( (void**)&(frame_gpu), sizeof(grid_param)),"Gradientgpu.cu : Alloc frame_gpu: " );
  cudasafe(cudaMalloc( (void**)&(grid_grad_x_gpu), (nbgridcells) * (nbgridcells) *sizeof(double)),"Gradientgpu.cu : Alloc source_x_gpu: " );
  cudasafe(cudaMalloc( (void**)&(grid_grad_y_gpu), (nbgridcells) * (nbgridcells) *sizeof(double)),"Gradientgpu.cu : Alloc source_y_gpu: " );

  // Copy values to the GPU
  cudasafe(cudaMemcpy(lens_x_gpu,lens[1].position_x , Nlens[1]*sizeof(double),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy x_gpu: " );
  cudasafe(cudaMemcpy(lens_y_gpu,lens[1].position_y , Nlens[1]*sizeof(double), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy y_gpu: " );
  cudasafe(cudaMemcpy(b0_gpu,lens[1].b0 , Nlens[1]*sizeof(double), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy b0_gpu: " );
  cudasafe(cudaMemcpy(angle_gpu,lens[1].ellipticity_angle , Nlens[1]*sizeof(double), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy angle_gpu: " );
  cudasafe(cudaMemcpy(epot_gpu, lens[1].ellipticity_potential, Nlens[1]*sizeof(double),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy epot_gpu: " );
  cudasafe(cudaMemcpy(rcore_gpu, lens[1].rcore, Nlens[1]*sizeof(double),cudaMemcpyHostToDevice ),"Gradientgpu.cu : Copy rcore_gpu: " );
  cudasafe(cudaMemcpy(rcut_gpu, lens[1].rcut, Nlens[1]*sizeof(double), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy rcut_gpu: " );
  cudasafe(cudaMemcpy(frame_gpu, frame, sizeof(grid_param), cudaMemcpyHostToDevice),"Gradientgpu.cu : Copy fame_gpu: " );

  //printf( "%d %d",nBlocks_gpu,(nbgridcells) * (nbgridcells)/threadsPerBlock );

  if (int((nbgridcells) * (nbgridcells)/threadsPerBlock) == 0){
    gradient_grid_piemd_GPU<<<1,threadsPerBlock>>>(grid_grad_x_gpu, grid_grad_y_gpu,frame_gpu,Nlens[1], nbgridcells, lens_x_gpu, lens_y_gpu, b0_gpu, angle_gpu, epot_gpu, rcore_gpu, rcut_gpu);
  }
  else{
  gradient_grid_piemd_GPU<<<(nbgridcells) * (nbgridcells)/threadsPerBlock,threadsPerBlock>>>(grid_grad_x_gpu, grid_grad_y_gpu,frame_gpu,Nlens[1], nbgridcells, lens_x_gpu, lens_y_gpu, b0_gpu, angle_gpu, epot_gpu, rcore_gpu, rcut_gpu);
  }
  cudasafe(cudaMemcpy( grid_grad_x, grid_grad_x_gpu, (nbgridcells) * (nbgridcells) *sizeof(double),cudaMemcpyDeviceToHost ),"Gradientgpu.cu : Copy source_x_gpu: " );
  cudasafe(cudaMemcpy( grid_grad_y, grid_grad_y_gpu, (nbgridcells) * (nbgridcells) *sizeof(double), cudaMemcpyDeviceToHost),"Gradientgpu.cu : Copy source_y_gpu: " );
  //printf ("%f %f \n", grid_grad_x[0],grid_grad_y[0]);

    // Free GPU memory
  cudaFree(lens_x_gpu);
  cudaFree(lens_y_gpu);
  cudaFree(b0_gpu);
  cudaFree(angle_gpu);
  cudaFree(epot_gpu);
  cudaFree(rcore_gpu);
  cudaFree(rcut_gpu);
  cudaFree(grid_grad_x_gpu);
  cudaFree(grid_grad_y_gpu);

  for (int i = 0; i < nbgridcells; i++){
    for(int j = 0; j < nbgridcells; j++){
      printf(" %f",grid_grad_x[i*nbgridcells + j]);
    }
    printf("\n");
  }
/*
  int nDevices;

    cudaGetDeviceCount(&nDevices);
    for (int i = 0; i < nDevices; i++) {
      cudaDeviceProp prop;
      cudaGetDeviceProperties(&prop, i);
      printf("Device Number: %d\n", i);
      printf("  Device name: %s\n", prop.name);
      printf("  Memory Clock Rate (KHz): %d\n",
             prop.memoryClockRate);
      printf("  Memory Bus Width (bits): %d\n",
             prop.memoryBusWidth);
      printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
             2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
    }

*/


}

void gradient_grid_pinned_multiple(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int *Nlens,int nbgridcells){
  //get number of GPU devices
  int nDevices;
  cudaGetDeviceCount(&nDevices);

  // Initialise kernel variables, table for multiple devices
  grid_param *frame_gpu[nDevices];
  double *lens_x_gpu[nDevices], *lens_y_gpu[nDevices], *b0_gpu[nDevices], *angle_gpu[nDevices], *epot_gpu[nDevices], *rcore_gpu[nDevices], *rcut_gpu[nDevices];
  double *grid_grad_x_gpu[nDevices], *grid_grad_y_gpu[nDevices] ;

  // Initialise multiple device variables
  int Ndevice[nDevices], indexactual[nDevices];
  cudaStream_t stream[nDevices];

  indexactual[0] = 0 ;
  Ndevice[0] = (nbgridcells) * (nbgridcells)/nDevices;
  printf("%d %d \n",indexactual[0], Ndevice[0]);

  for (int dev = 1; dev < nDevices; dev++) {

    Ndevice[dev] = (nbgridcells) * (nbgridcells)/nDevices;

    if(indexactual[dev]+Ndevice[dev] > (nbgridcells) * (nbgridcells)){
      Ndevice[dev] = (nbgridcells) * (nbgridcells) - indexactual[dev-1];
    }

    indexactual[dev] = indexactual[dev-1] + Ndevice[dev];
    printf("%d %d \n",indexactual[dev], Ndevice[dev]);
  }

  for (int dev = 0; dev < nDevices; dev++) {


    cudaSetDevice(dev);

    //SIS//


    //PIEMD8//

    // Allocate variables on the GPU
    cudasafe(cudaMalloc( (void**)&(lens_x_gpu[dev]), Nlens[1]*sizeof(double)),"Gradientgpu.cu : Alloc x_gpu: " );
    cudasafe(cudaMalloc( (void**)&(lens_y_gpu[dev]), Nlens[1]*sizeof(double)),"Gradientgpu.cu : Alloc y_gpu: " );
    cudasafe(cudaMalloc( (void**)&(b0_gpu[dev]), Nlens[1]*sizeof(double)),"Gradientgpu.cu : Alloc b0_gpu: " );
    cudasafe(cudaMalloc( (void**)&(angle_gpu[dev]), Nlens[1]*sizeof(double)),"Gradientgpu.cu : Alloc angle_gpu: " );
    cudasafe(cudaMalloc( (void**)&(epot_gpu[dev]), Nlens[1]*sizeof(double)),"Gradientgpu.cu : Alloc epot_gpu: " );
    cudasafe(cudaMalloc( (void**)&(rcore_gpu[dev]), Nlens[1]*sizeof(double)),"Gradientgpu.cu : Alloc rcore_gpu: " );
    cudasafe(cudaMalloc( (void**)&(rcut_gpu[dev]), Nlens[1]*sizeof(double)),"Gradientgpu.cu : Alloc rcut_gpu: " );
    cudasafe(cudaMalloc( (void**)&(frame_gpu[dev]), sizeof(grid_param)),"Gradientgpu.cu : Alloc frame_gpu: " );
    cudasafe(cudaMalloc( (void**)&(grid_grad_x_gpu[dev]), Ndevice[dev] *sizeof(double)),"Gradientgpu.cu : Alloc source_x_gpu: " );
    cudasafe(cudaMalloc( (void**)&(grid_grad_y_gpu[dev]), Ndevice[dev] *sizeof(double)),"Gradientgpu.cu : Alloc source_y_gpu: " );


    cudaStreamCreate(&stream[dev]);

  }

  for (int dev = 0; dev < nDevices; dev++) {


    cudaSetDevice(dev);

    // Copy values to the GPU
    cudasafe(cudaMemcpyAsync(lens_x_gpu[dev],lens[1].position_x , Nlens[1]*sizeof(double),cudaMemcpyHostToDevice,stream[dev] ),"Gradientgpu.cu : Copy x_gpu: " );
    cudasafe(cudaMemcpyAsync(lens_y_gpu[dev],lens[1].position_y , Nlens[1]*sizeof(double), cudaMemcpyHostToDevice,stream[dev]),"Gradientgpu.cu : Copy y_gpu: " );
    cudasafe(cudaMemcpyAsync(b0_gpu[dev],lens[1].b0 , Nlens[1]*sizeof(double), cudaMemcpyHostToDevice,stream[dev]),"Gradientgpu.cu : Copy b0_gpu: " );
    cudasafe(cudaMemcpyAsync(angle_gpu[dev],lens[1].ellipticity_angle , Nlens[1]*sizeof(double), cudaMemcpyHostToDevice,stream[dev]),"Gradientgpu.cu : Copy angle_gpu: " );
    cudasafe(cudaMemcpyAsync(epot_gpu[dev], lens[1].ellipticity_potential, Nlens[1]*sizeof(double),cudaMemcpyHostToDevice ,stream[dev]),"Gradientgpu.cu : Copy epot_gpu: " );
    cudasafe(cudaMemcpyAsync(rcore_gpu[dev], lens[1].rcore, Nlens[1]*sizeof(double),cudaMemcpyHostToDevice ,stream[dev]),"Gradientgpu.cu : Copy rcore_gpu: " );
    cudasafe(cudaMemcpyAsync(rcut_gpu[dev], lens[1].rcut, Nlens[1]*sizeof(double), cudaMemcpyHostToDevice,stream[dev]),"Gradientgpu.cu : Copy rcut_gpu: " );
    cudasafe(cudaMemcpyAsync(frame_gpu[dev], frame, sizeof(grid_param), cudaMemcpyHostToDevice,stream[dev]),"Gradientgpu.cu : Copy fame_gpu: " );
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
      gradient_grid_piemd_GPU_multiple<<<1,threadsPerBlock,0,stream[dev]>>>(grid_grad_x_gpu[dev], grid_grad_y_gpu[dev],frame_gpu[dev],Nlens[1], nbgridcells,Ndevice[dev],indexactual[dev], lens_x_gpu[dev], lens_y_gpu[dev], b0_gpu[dev], angle_gpu[dev], epot_gpu[dev], rcore_gpu[dev], rcut_gpu[dev]);
    }
    else{
      gradient_grid_piemd_GPU_multiple<<<(nbgridcells) * (nbgridcells)/threadsPerBlock,threadsPerBlock,0,stream[dev]>>>(grid_grad_x_gpu[dev], grid_grad_y_gpu[dev],frame_gpu[dev],Nlens[1], nbgridcells, Ndevice[dev],indexactual[dev], lens_x_gpu[dev], lens_y_gpu[dev], b0_gpu[dev], angle_gpu[dev], epot_gpu[dev], rcore_gpu[dev], rcut_gpu[dev]);
    }
  }

  for (int dev = 0; dev < nDevices; dev++) {
    cudasafe(cudaMemcpyAsync( grid_grad_x + indexactual[dev], grid_grad_x_gpu[dev], Ndevice[dev] *sizeof(double),cudaMemcpyDeviceToHost ,stream[dev]),"Gradientgpu.cu : Copy source_x_gpu: " );
    cudasafe(cudaMemcpyAsync( grid_grad_y + indexactual[dev], grid_grad_y_gpu[dev], Ndevice[dev] *sizeof(double), cudaMemcpyDeviceToHost,stream[dev]),"Gradientgpu.cu : Copy source_y_gpu: " );
  }

  for (int dev = 0; dev < nDevices; dev++) {


    cudaSetDevice(dev);
      // Free GPU memory
    cudaFree(lens_x_gpu[dev]);
    cudaFree(lens_y_gpu[dev]);
    cudaFree(b0_gpu[dev]);
    cudaFree(angle_gpu[dev]);
    cudaFree(epot_gpu[dev]);
    cudaFree(rcore_gpu[dev]);
    cudaFree(rcut_gpu[dev]);
    cudaFree(grid_grad_x_gpu[dev]);
    cudaFree(grid_grad_y_gpu[dev]);
    cudaStreamDestroy(stream[dev]);

  }
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


__global__ void gradient_grid_sis_GPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcells, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut) {

  //*grad_x = *grad_y = 0.;

  int bid=blockIdx.x; // index of the block (and of the set of images)
  int tid=threadIdx.x; // index of the thread within the block

  double dx,dy,x_pos,y_pos;        //pixelsize
  int grid_dim, index;
  struct point true_coord, true_coord_rotation;
  double      R;
  dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
  dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
  grid_dim = (nbgridcells);

  index = bid * threadsPerBlock + tid ;



  while(index < grid_dim*grid_dim){

    grid_grad_x[index] = 0.;
    grid_grad_y[index] = 0.;

    x_pos= frame->xmin + (index/grid_dim)*dx;
    y_pos= frame->ymin + (index % grid_dim)*dy;
    //printf("%f %f \n", x_pos, y_pos);

    for(int iterator=0; iterator < Nlens; iterator++){

      true_coord.x = x_pos - lens_x[iterator];
      true_coord.y = y_pos - lens_y[iterator];

      // Change the origin of the coordinate system to the center of the clump
      true_coord_rotation = rotateCoordinateSystem_GPU(true_coord, lens_angle[iterator]);

      R=sqrt(true_coord_rotation.x*true_coord_rotation.x*(1-lens_epot[iterator])+true_coord_rotation.y*true_coord_rotation.y*(1+lens_epot[iterator]));  //ellippot = ellipmass/3
      grid_grad_x[index] += (1-lens_epot[iterator])*lens_b0[iterator]*true_coord_rotation.x/(R);
      grid_grad_y[index] += (1+lens_epot[iterator])*lens_b0[iterator]*true_coord_rotation.y/(R);



    }




    bid += gridDim.x;
    index = bid * threadsPerBlock + tid;
  }


}

__global__ void gradient_grid_piemd_GPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcells, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut) {

  //*grad_x = *grad_y = 0.;

  int bid=blockIdx.x; // index of the block (and of the set of images)
  int tid=threadIdx.x; // index of the thread within the block

  double dx,dy,x_pos,y_pos;        //pixelsize
  int grid_dim, index;
  struct point true_coord, true_coord_rotation;
  complex      zis;
  dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
  dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
  grid_dim = (nbgridcells);

  index = bid * threadsPerBlock + tid ;



  while(index < grid_dim*grid_dim){

    //grid_grad_x[index] = 0.;
    //grid_grad_y[index] = 0.;

    x_pos= frame->xmin + (index/grid_dim)*dx;
    y_pos= frame->ymin + (index % grid_dim)*dy;
    //printf("%f %f \n", x_pos, y_pos);

    for(int iterator=0; iterator < Nlens; iterator++){

      true_coord.x = x_pos - lens_x[iterator];
      true_coord.y = y_pos - lens_y[iterator];

      // Change the origin of the coordinate system to the center of the clump
      true_coord_rotation = rotateCoordinateSystem_GPU(true_coord, lens_angle[iterator]);

      double x   = true_coord_rotation.x;
      double y   = true_coord_rotation.y;
      double eps = lens_epot[iterator];
      double rc  = lens_rcore[iterator];

      double sqe  = sqrt(eps);
      //
      double cx1  = (1. - eps)/(1. + eps);
      double cxro = (1. + eps)*(1. + eps);
      double cyro = (1. - eps)*(1. - eps);
      //
      double rem2 = x*x/cxro + y*y/cyro;

      complex zci, znum, zden, zres;
      double norm;
      //
      zci.re  = 0;
      zci.im  = -0.5*(1. - eps*eps)/sqe;
      //
      znum.re = cx1*x;
      znum.im = 2.*sqe*sqrt(rc*rc + rem2) - y/cx1;
      //
      zden.re = x;
      zden.im = 2.*rc*sqe - y;
      norm    = (zden.re*zden.re + zden.im*zden.im);     // zis = znum/zden
      //
      zis.re  = (znum.re*zden.re + znum.im*zden.im)/norm;
      zis.im  = (znum.im*zden.re - znum.re*zden.im)/norm;
      norm    = zis.re;
      zis.re  = log(sqrt(norm*norm + zis.im*zis.im));  // ln(zis) = ln(|zis|)+i.Arg(zis)
      zis.im  = atan2(zis.im, norm);
      //  norm = zis.re;
      zres.re = zci.re*zis.re - zci.im*zis.im;   // Re( zci*ln(zis) )
      zres.im = zci.im*zis.re + zis.im*zci.re;   // Im( zci*ln(zis) )
      //
      zis.re  = zres.re;
      zis.im  = zres.im;

      grid_grad_x[index] += lens_b0[iterator]*zis.re;
      grid_grad_y[index] += lens_b0[iterator]*zis.im;



    }




    bid += gridDim.x;
    index = bid * threadsPerBlock + tid;
  }


}

__global__ void gradient_grid_piemd_GPU_multiple(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcells, int Ndevice, int indexactual, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut) {

  //*grad_x = *grad_y = 0.;

  int bid=blockIdx.x; // index of the block (and of the set of images)
  int tid=threadIdx.x; // index of the thread within the block

  double dx,dy,x_pos,y_pos;        //pixelsize
  int grid_dim, index;
  struct point true_coord, true_coord_rotation;
  complex      zis;
  dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
  dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
  grid_dim = (nbgridcells);

  index = bid * threadsPerBlock + tid ;



  while(index < Ndevice){
    //printf("%d \n", index);

    grid_grad_x[index] = 0.;
    grid_grad_y[index] = 0.;

    x_pos= frame->xmin + ((indexactual +index)/grid_dim)*dx;
    y_pos= frame->ymin + ((indexactual +index) % grid_dim)*dy;
    //printf("%f %f \n", x_pos, y_pos);

    for(int iterator=0; iterator < Nlens; iterator++){

      true_coord.x = x_pos - lens_x[iterator];
      true_coord.y = y_pos - lens_y[iterator];

      // Change the origin of the coordinate system to the center of the clump
      true_coord_rotation = rotateCoordinateSystem_GPU(true_coord, lens_angle[iterator]);

      double x   = true_coord_rotation.x;
      double y   = true_coord_rotation.y;
      double eps = lens_epot[iterator];
      double rc  = lens_rcore[iterator];

      double sqe  = sqrt(eps);
      //
      double cx1  = (1. - eps)/(1. + eps);
      double cxro = (1. + eps)*(1. + eps);
      double cyro = (1. - eps)*(1. - eps);
      //
      double rem2 = x*x/cxro + y*y/cyro;

      complex zci, znum, zden, zres;
      double norm;
      //
      zci.re  = 0;
      zci.im  = -0.5*(1. - eps*eps)/sqe;
      //
      znum.re = cx1*x;
      znum.im = 2.*sqe*sqrt(rc*rc + rem2) - y/cx1;
      //
      zden.re = x;
      zden.im = 2.*rc*sqe - y;
      norm    = (zden.re*zden.re + zden.im*zden.im);     // zis = znum/zden
      //
      zis.re  = (znum.re*zden.re + znum.im*zden.im)/norm;
      zis.im  = (znum.im*zden.re - znum.re*zden.im)/norm;
      norm    = zis.re;
      zis.re  = log(sqrt(norm*norm + zis.im*zis.im));  // ln(zis) = ln(|zis|)+i.Arg(zis)
      zis.im  = atan2(zis.im, norm);
      //  norm = zis.re;
      zres.re = zci.re*zis.re - zci.im*zis.im;   // Re( zci*ln(zis) )
      zres.im = zci.im*zis.re + zis.im*zci.re;   // Im( zci*ln(zis) )
      //
      zis.re  = zres.re;
      zis.im  = zres.im;

      grid_grad_x[index] += lens_b0[iterator]*zis.re;
      grid_grad_y[index] += lens_b0[iterator]*zis.im;



    }




    bid += gridDim.x;

    index = bid * threadsPerBlock + tid;
  }


}



__global__ void piemd_GPU(double *grad_x, double *grad_y, point *pImage, int Nlens, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut) {

  //*grad_x = *grad_y = 0.;

  int bid=blockIdx.x; // index of the block (and of the set of images)
  int tid=threadIdx.x; // index of the thread within the block

  struct point true_coord, true_coord_rotation;
  complex      zis;
  //gridDim.x , threadsPerBlock;

  int iterator = bid * threadsPerBlock + tid ;
  //for(int iterator = 0; iterator < Nlens; ++iterator){
  while(iterator < Nlens){
    //printf(" %d %d %d %d \n",iterator , bid , threadsPerBlock , tid );
    true_coord.x = pImage->x - lens_x[iterator];
    true_coord.y = pImage->y - lens_y[iterator];

    // Change the origin of the coordinate system to the center of the clump
    true_coord_rotation = rotateCoordinateSystem_GPU(true_coord, lens_angle[iterator]);

    double x   = true_coord_rotation.x;
    double y   = true_coord_rotation.y;
    double eps = lens_epot[iterator];
    double rc  = lens_rcore[iterator];

    double sqe  = sqrt(eps);
    //
    double cx1  = (1. - eps)/(1. + eps);
    double cxro = (1. + eps)*(1. + eps);
    double cyro = (1. - eps)*(1. - eps);
    //
    double rem2 = x*x/cxro + y*y/cyro;

    complex zci, znum, zden, zres;
    double norm;
    //
    zci.re  = 0;
    zci.im  = -0.5*(1. - eps*eps)/sqe;
    //
    znum.re = cx1*x;
    znum.im = 2.*sqe*sqrt(rc*rc + rem2) - y/cx1;
    //
    zden.re = x;
    zden.im = 2.*rc*sqe - y;
    norm    = (zden.re*zden.re + zden.im*zden.im);     // zis = znum/zden
    //
    zis.re  = (znum.re*zden.re + znum.im*zden.im)/norm;
    zis.im  = (znum.im*zden.re - znum.re*zden.im)/norm;
    norm    = zis.re;
    zis.re  = log(sqrt(norm*norm + zis.im*zis.im));  // ln(zis) = ln(|zis|)+i.Arg(zis)
    zis.im  = atan2(zis.im, norm);
    //  norm = zis.re;
    zres.re = zci.re*zis.re - zci.im*zis.im;   // Re( zci*ln(zis) )
    zres.im = zci.im*zis.re + zis.im*zci.re;   // Im( zci*ln(zis) )
    //
    zis.re  = zres.re;
    zis.im  = zres.im;
    //
    //zres.re = zis.re*b0;
    //zres.im = zis.im*b0;
    //
    //

    //grad->x += lens_b0[iterator]*zis.re;
    //grad->y += lens_b0[iterator]*zis.im;
    atomicAdd_double(grad_x,lens_b0[iterator]*zis.re);
    atomicAdd_double(grad_y,lens_b0[iterator]*zis.im);

    //printf("%f %f \n ", lens_b0[iterator]*zis.re, lens_b0[iterator]*zis.im);
    //printf("%f %f \n ", *grad_x, *grad_y);

    bid += gridDim.x;
    iterator = bid * threadsPerBlock + tid ;

  }


}

__device__ static double atomicAdd_double(double* address, double val)
{
    unsigned long long int* address_as_ull =
                                          (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                        __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

