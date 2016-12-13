/**
 * @file   main.cpp
 * @Author Christoph Schaaefer, EPFL (christophernstrerne.schaefer@epfl.ch)
 * @date   October 2016
 * @brief  Benchmark for gradhalo function
 */

#include <iostream>
#include <iomanip>
#include <string.h>
#include <cuda_runtime.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
//
#include <mm_malloc.h>
//
//#include "../../../Projects/lenstool-6.8.1/include/structure.h"
//#include "structure.h"
#include "structure_hpc.h"
#include "timer.h"
#include "gradient.hpp"
#include "gradient_avx.hpp"
#include "gradient_avx512f.hpp"
#include "setup.hpp"
//
#define NN 1000
//
//struct pot      lens_ref[NLMAX];
//
int main()
{
	double t1, t2;
	//
	//Variable creation
	//
	point image;
	//
	//Initialisation
	//
        int nlenses;
        double x, y;
        double sol_grad_x, sol_grad_y;
	point grad; // store the result
        //
        Potential*   lens;
	//
	setup_jauzac(&lens, &nlenses, &image.x, &image.y, &sol_grad_x, &sol_grad_y);
	t1 = -myseconds();
	for (int ii = 0; ii < NN; ++ii)
		grad = module_potentialDerivatives_totalGradient(nlenses, &image, lens);
        t1 += myseconds();
	//printf("---> grad = %f %f\n", grad.x, grad.y);
	//
	// Setting up the AOS potential 
	//
	point image_soa;
	int nlenses_soa;
	double x_soa, y_soa;
	double sol_grad_x_soa, sol_grad_y_soa;
	//
	Potential_SOA   lens_soa;
	setup_jauzac_SOA(&lens_soa, &nlenses_soa, &image_soa.x, &image_soa.y, &sol_grad_x_soa, &sol_grad_y_soa);
	//
	std::cout << "Benchmark for Gradient Calculation using  " << nlenses_soa << " lenses, image: " << image_soa.x << " " << image_soa.y << std::endl;
	point grad_soa; // store the result
	//
	//gettimeofday(&t1, 0);
	//
	//module_potentialDerivatives_totalGradient(&runmodesmall,&image, lens);
	//
	t2 = -myseconds();
	for (int ii = 0; ii < NN; ++ii)
	{
		grad_soa = module_potentialDerivatives_totalGradient_81_SOA_AVX512(&image_soa, &lens_soa, nlenses_soa);	 
		//grad_soa = module_potentialDerivatives_totalGradient_81_SOA_AVX(&image_soa, &lens_soa, nlenses_soa);	 
		//grad_soa = module_potentialDerivatives_totalGradient_SOA_AVX(&image_soa, &lens_soa, nlenses_soa);	 
	}	
	t2 += myseconds();
	//
	std::cout << " ref sol  = " << std::setprecision(15) << sol_grad_x << " " << std::setprecision(15) << sol_grad_y << std::endl;
	std::cout << " grad     = " << std::setprecision(15) << grad.x << " " << std::setprecision(15) << grad.y << ", time = " << t1 << " s. " << std::endl;
	std::cout << " grad avx = " << std::setprecision(15) << grad_soa.x << " " << std::setprecision(15) << grad_soa.y << ", time = " << t2 << " s. " << std::endl;
	//
	//
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
	free(lens);
}

