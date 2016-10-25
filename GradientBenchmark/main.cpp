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
#include "structure.h"
#include "timer.h"
//
/** for both gradient and second derivatives **/
static struct point rotateCoordinateSystem(struct point P, double theta);

/** gradient **/
struct point module_potentialDerivatives_totalGradient(const runmode_param *runmode, const struct point *pImage, const struct Potential *lens);
struct point potentialDerivatives_totalGradient(const runmode_param *runmode, const struct point *pImage, const struct Potential *lens, int nhalos);
//
static struct point grad_halo(const struct point *pImage, const struct Potential *lens);

/** PIEMD **/
static complex piemd_1derivatives_ci05(double x, double y, double eps, double rc);

/** Potential **/
void module_readParameters_calculatePotentialparameter(Potential *lens);

int main()
{
	//Constant
	int small(10);
	int medium(100);
	int big(1000);
	//
	//Variable creation
	//
	runmode_param runmodesmall;
	runmode_param runmodemedium;
	runmode_param runmodebig;

	point image;

	Potential  *ilens;
	Potential lens[big];
	//
	//Initialisation
	//
	runmodesmall.nhalos  = small;
	runmodemedium.nhalos = medium;
	runmodebig.nhalos    = big;
	//
	image.x = image.y = 2;
	//
	for (int i = 0; i <big; ++i)
	{
		ilens = &lens[i];
		//
		ilens->vdisp       = 1.;
		ilens->position.x  = ilens->position.y = 0.;
		ilens->type        = 8;
		ilens->ellipticity = 0.11;
		ilens->ellipticity_potential = 0.;
		ilens->ellipticity_angle = 0.;
		ilens->rcut 	   = 5.;
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
	grad = potentialDerivatives_totalGradient(&runmodebig, &image, lens, big);
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

