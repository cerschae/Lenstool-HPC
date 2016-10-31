/**
 * @file   main.cpp
 * @Author Christoph Schaaefer, EPFL (christophernstrerne.schaefer@epfl.ch)
 * @date   October 2016
 * @brief  Benchmark for gradhalo function
 */

#include <iostream>
#include <string.h>
#include <cuda_runtime.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
//
#include <mm_malloc.h>
//
#include "structure.h"
#include "timer.h"
#include "gradient.hpp"
//
int main()
{
	//Constant
	int small(10);
	int medium(100);
	int big(10000);
	//
	//Variable creation
	//
	runmode_param runmodesmall;
	runmode_param runmodemedium;
	runmode_param runmodebig;
	//
	point image;
	image.x = image.y = 2;

	//
	//
	//Initialisation
	//
	runmodesmall.nhalos  = small;
	runmodemedium.nhalos = medium;
	runmodebig.nhalos    = big;
	//
	// Setting up the AOS potential 
	//
	Potential   lens[big];
	Potential* ilens = &lens[0];
	for (int i = 0; i < big; ++i, ++ilens)
	{
		//ilens = &lens[i];
		//
		ilens->vdisp       = 1.;
		ilens->position.x  = ilens->position.y = 0.;
		ilens->type        = 8;
		ilens->ellipticity = 0.11;
		ilens->ellipticity_potential = 0.;
		ilens->ellipticity_angle = 0.;
		ilens->rcut        = 5.;
		ilens->rcore       = 1;
		ilens->weight      = 0;
		ilens->rscale      = 0;
		ilens->exponent    = 0;
		ilens->alpha       = 0.;
		ilens->einasto_kappacritic = 0;
		ilens->z           = 0.4;
		module_readParameters_calculatePotentialparameter(ilens);
	}	
	//
	// Setting up the SOA potential
	//
	Potential_SOA lens_soa;
	for (int i = 0; i < big; ++i)
	{
		lens_soa.type       = (int*) malloc(big*sizeof(int));
		lens_soa.vdisp      = (double*) malloc(big*sizeof(double));
		//
		lens_soa.position_x = (double*) malloc(big*sizeof(double));
		lens_soa.position_y = (double*) malloc(big*sizeof(double));
		lens_soa.weight = (double*) malloc(big*sizeof(double));
                lens_soa.b0 = (double*) malloc(big*sizeof(double));
                lens_soa.vdisp = (double*) malloc(big*sizeof(double));
                lens_soa.ellipticity_angle = (double*) malloc(big*sizeof(double));
                lens_soa.ellipticity = (double*) malloc(big*sizeof(double));
                lens_soa.ellipticity_potential = (double*) malloc(big*sizeof(double));
                lens_soa.rcore = (double*) malloc(big*sizeof(double));
                lens_soa.rcut = (double*) malloc(big*sizeof(double));
                lens_soa.rscale = (double*) malloc(big*sizeof(double));
                lens_soa.exponent = (double*) malloc(big*sizeof(double));
                lens_soa.alpha = (double*) malloc(big*sizeof(double));
                lens_soa.einasto_kappacritic = (double*) malloc(big*sizeof(double));
                lens_soa.z = (double*) malloc(big*sizeof(double));
                lens_soa.mag = (double*) malloc(big*sizeof(double));
                lens_soa.lum = (double*) malloc(big*sizeof(double));
                lens_soa.theta = (double*) malloc(big*sizeof(double));
                lens_soa.sigma = (double*) malloc(big*sizeof(double));
	}
	//
	for (int i = 0; i < big; ++i)
	{
		//ilens = &lens[i];
		//
		lens_soa.vdisp[i]         	  = 1.;
		//ilens->position.x  = ilens->position.y = 0.;
		lens_soa.position_x[i]            = 0.;
		lens_soa.position_y[i]  	  = 0.;
		lens_soa.type[i]        	  = 8;
		lens_soa.ellipticity[i] 	  = 0.11;
		lens_soa.ellipticity_potential[i] = 0.;
		lens_soa.ellipticity_angle[i]     = 0.;
		lens_soa.rcut[i] 	          = 5.;
		lens_soa.rcore[i]                 = 1;
		lens_soa.weight[i]                = 0;
		lens_soa.rscale[i]                = 0;
		lens_soa.exponent[i]              = 0;
		lens_soa.alpha[i]                 = 0.;
		lens_soa.einasto_kappacritic[i]   = 0;
		lens_soa.z[i]                     = 0.4;
		module_readParameters_calculatePotentialparameter_SOA(&lens_soa, i);
	}
	//
	std::cout << "Benchmark for Gradient Calculation "<< std::endl;
	point grad; // store the result
	//
	//gettimeofday(&t1, 0);
	//
	//module_potentialDerivatives_totalGradient(&runmodesmall,&image, lens);
	double t1;
	//
	t1 = -myseconds();
	grad = module_potentialDerivatives_totalGradient(&runmodebig, &image, lens);
	t1 += myseconds();
	std::cout << "Sample size " << big << ": " << t1 << "s., point = " << grad.x << " " << grad.y << std::endl;
	//
	t1 = -myseconds();
	grad = module_potentialDerivatives_totalGradient_SOA_AVX(&runmodebig, &image, &lens_soa, big);
	t1 += myseconds();
	std::cout << "Sample size " << big << ": " << t1 << "s., point = " << grad.x << " " << grad.y << std::endl;
	//
#if 0
	std::ofstream myfile;
	myfile.open ("BenchmarkGrad.txt");
	myfile << "Benchmark for Gradient Calculation "<< std::endl;
	myfile << "Sample size " << small << ": " << time1 << std::endl;
	myfile << "Sample size " << medium << ": " << time2 << std::endl;
	myfile << "Sample size " << big << ": " << time3 << std::endl;
	myfile.close();
#endif
}

