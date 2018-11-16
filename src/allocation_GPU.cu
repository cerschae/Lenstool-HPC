//
// This file is part of lenstoolhpc
// authors: gilles.fourestey@epfl.ch
// 
#include <math.h>
#include "module_cosmodistances.hpp"
#include "module_readParameters.hpp"
#ifdef __WITH_GPU
#include <cuda_runtime.h>
#include <cuda.h>
#endif

#include "allocation.hpp"

void PotentialSOAAllocation_GPU(Potential_SOA **lens_SOA, const int nhalos)
{
	
	Potential_SOA* p;
#if (defined __WITH_UM) 
#warning "Allocation Using unified memory"
	printf("--- Unified memory Allocation FTW\n");
	cudaError error = cudaMallocManaged(lens_SOA, sizeof(Potential_SOA));
	//if (error == 0) printf("Allocation error\n");
	p = *lens_SOA;
        cudaMallocManaged(&(p->type), nhalos*sizeof(int));
        cudaMallocManaged(&p->position_x         , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->position_y           , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->b0                   , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->ellipticity_angle    , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->ellipticity          , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->ellipticity_potential, nhalos*sizeof(size_t));
        cudaMallocManaged(&p->rcore                , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->rcut                 , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->vdisp                 , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->z                    , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->anglecos             , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->anglesin             , nhalos*sizeof(size_t));
	cudaMallocManaged(&p->SOA_index            , nhalos*sizeof(int));
	cudaMallocManaged(&p->N_types              , 100*sizeof(int));
#else
	PotentialSOAAlocation(lens_SOA, nhalos);
#if 0
	*lens_SOA = (Potential_SOA*) malloc(sizeof(Potential_SOA));
	p = *lens_SOA;
        p->type 		 = (int*) malloc(sizeof(int)*nhalos);
	p->position_x 		 = (double*) malloc(sizeof(double)*nhalos);
	p->position_y 		 = (double*) malloc(sizeof(double)*nhalos);
	p->b0 			 = (double*) malloc(sizeof(double)*nhalos);
	p->ellipticity_angle 	 = (double*) malloc(sizeof(double)*nhalos);
	p->ellipticity 		 = (double*) malloc(sizeof(double)*nhalos);
	p->ellipticity_potential = (double*) malloc(sizeof(double)*nhalos);
	p->rcore 		 = (double*) malloc(sizeof(double)*nhalos);
	p->rcut 		 = (double*) malloc(sizeof(double)*nhalos);
	p->vdisp		 = (double*) malloc(sizeof(double)*nhalos);
	p->z 	                 = (double*) malloc(sizeof(double)*nhalos);
	p->anglecos 		 = (double*) malloc(sizeof(double)*nhalos);
	p->anglesin 		 = (double*) malloc(sizeof(double)*nhalos);
	p->SOA_index 		 = (int*) malloc(sizeof(int)*nhalos);
	p->N_types 		 = (int*) malloc(sizeof(int)*100);
#endif
#endif
}



void PotentialSOADeallocation_GPU(Potential_SOA *lens_SOA)
{
#if (defined __WITH_UM)
#warning "Using unified memory"
	printf("--- Unified memory Deallocation FTW\n");
	cudaFree(lens_SOA->type);
	cudaFree(lens_SOA->position_x);
	cudaFree(lens_SOA->position_y);
	cudaFree(lens_SOA->b0);
	cudaFree(lens_SOA->ellipticity_angle);
	cudaFree(lens_SOA->ellipticity);
	cudaFree(lens_SOA->ellipticity_potential);
	cudaFree(lens_SOA->rcore);
	cudaFree(lens_SOA->rcut);
	cudaFree(lens_SOA->vdisp);
	cudaFree(lens_SOA->z);
	cudaFree(lens_SOA->anglecos);
	cudaFree(lens_SOA->anglesin);
	cudaFree(lens_SOA->SOA_index);
#else
	PotentialSOADeallocation_GPU(lens_SOA);
#if 0
	free(lens_SOA->type);
	free(lens_SOA->position_x);
	free(lens_SOA->position_y);
	free(lens_SOA->b0);
	free(lens_SOA->ellipticity_angle);
	free(lens_SOA->ellipticity);
	free(lens_SOA->ellipticity_potential);
	free(lens_SOA->rcore);
	free(lens_SOA->rcut);             
	free(lens_SOA->z);             
	free(lens_SOA->anglecos);
	free(lens_SOA->anglesin);       
	free(lens_SOA);
#endif
#endif
}
