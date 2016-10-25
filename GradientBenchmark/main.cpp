/**
* @file   main.cpp
* @Author Christoph Schaaefer, EPFL (christophernstrerne.schaefer@epfl.ch)
* @date   October 2016
* @brief  Benchmark for gradhalo function
*/

#include <iostream>
#include <string.h>
#include "structure.h"
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <Grad.h>
#include <GradTest.h>


int main()
{

//Constant
int small(10);
int medium(100);
int big(1000);

//Variable creation
struct timeval t1, t2, t3, t4;
runmode_param runmodesmall;
runmode_param runmodemedium;
runmode_param runmodebig;

point image;

Potential  *ilens;
Potential lens[big];

//Initialisation

runmodesmall.nhalos = small;
runmodemedium.nhalos = medium;
runmodebig.nhalos = big;
image.x = image.y = 2;

for (int i = 0; i <big; ++i){
	ilens = &lens[i];
	
    ilens->position.x = ilens->position.y = 0.;
    ilens->type = 8;
    ilens->ellipticity = 0.11;
    ilens->ellipticity_potential = 0.;
    ilens->ellipticity_angle = 0.;
    ilens->vdisp = 1.;
    ilens->rcut = 5.;
    ilens->rcore = 1;
    ilens->weight = 0;
    ilens->rscale = 0;
    ilens->exponent = 0;
    ilens->alpha = 0.;
    ilens->einasto_kappacritic = 0;
    ilens->z = 0.4;
    module_readParameters_calculatePotentialparameter(ilens);

}

/** SoA part  **/

//Init PotentialSet

PotentialSet lenses;
lenses.type = 	new int[big];
lenses.x  = 	new double[big];
lenses.y = 		new double[big];
lenses.b0 = 	new double[big];
lenses.ellipticity_angle = new double[big];
lenses.ellipticity = new double[big];
lenses.ellipticity_potential = new double[big];
lenses.rcore = 	new double[big];
lenses.rcut = 	new double[big];
lenses.z = 		new double[big];

for (int i = 0; i <big; ++i){
	lenses.type[i] = 	lens[i].type;
	lenses.x[i]  = 		lens[i].position.x;
	lenses.y[i] = 		lens[i].position.y;
	lenses.b0[i] = 		lens[i].b0;
	lenses.ellipticity_angle[i] = lens[i].ellipticity_angle;
	lenses.ellipticity[i] = lens[i].ellipticity;
	lenses.ellipticity_potential[i] = lens[i].ellipticity_potential;
	lenses.rcore[i] = 	lens[i].rcore;
	lenses.rcut[i] = 	lens[i].rcut;
	lenses.z[i] = 		lens[i].z;

}


gettimeofday(&t1, 0);
module_potentialDerivatives_totalGradient(&runmodesmall,&image, &lenses);
gettimeofday(&t2, 0);
module_potentialDerivatives_totalGradient(&runmodemedium,&image, &lenses);
gettimeofday(&t3, 0);
point grad = module_potentialDerivatives_totalGradient(&runmodebig,&image, &lenses);
gettimeofday(&t4, 0);

double time1 = (1000000.0*(t2.tv_sec-t1.tv_sec) + t2.tv_usec-t1.tv_usec)/1000000.0;
double time2 = (1000000.0*(t3.tv_sec-t2.tv_sec) + t3.tv_usec-t2.tv_usec)/1000000.0;
double time3 = (1000000.0*(t4.tv_sec-t3.tv_sec) + t4.tv_usec-t3.tv_usec)/1000000.0;



std::cout << "Benchmark for Gradient SOA Calculation "<< std::endl;
std::cout << "Sample size " << small << ": " << time1 << std::endl;
std::cout << "Sample size " << medium << ": " << time2 << std::endl;
std::cout << "Sample size " << big << ": " << time3 << std::endl;
std::cout << "Grad " << grad.x << " and " << grad.y << std::endl;

gradtest();

std::ofstream myfile;
myfile.open ("BenchmarkGradSoA.txt");
myfile << "Benchmark for Gradient SOA Calculation "<< std::endl;
myfile << "Sample size " << small << ": " << time1 << std::endl;
myfile << "Sample size " << medium << ": " << time2 << std::endl;
myfile << "Sample size " << big << ": " << time3 << std::endl;
myfile.close();



}
