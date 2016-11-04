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

int main()
{

//Constant
int small(10);
int medium(100);
int big(1000);

//Number of SIS and PIEMD
int Nsis(0);
int Npiemd(0);

//Variable creation
struct timeval t1, t2, t3, t4;

point image;

Potential  *ilens;
Potential lens[big];

//Initialisation

image.x = image.y = 2;

for (int i = 0; i <big; ++i){

	ilens = &lens[i];
	
    ilens->position.x = ilens->position.y = 0.;
    if ( i < 0 ){
    	ilens->type = 5;
    	Nsis +=1;
    }
    else{
    	ilens->type = 8;
    	Npiemd +=1;
    }
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

/** SoA part + Sorting by lens type:**/

/** Sorting by lens type will remove an if condition from our loop.
 * Putting "if" inside loops prevents certain types of optimization like e.g. loop unrolling to expose
 * Instruction Level Parallelism. This is particularly true for GPUs. **/

//Init PotentialSet
PotentialSet lenses[2];
int Nlensessmall[2];
int Nlensesmedium[2];
int Nlensesbig[2];
PotentialSet lensespiemd;
PotentialSet lensessis;
PotentialSet  *ilenses;

lensespiemd.type = 	new int[Npiemd];
lensespiemd.x  = 	new double[Npiemd];
lensespiemd.y = 		new double[Npiemd];
lensespiemd.b0 = 	new double[Npiemd];
lensespiemd.ellipticity_angle = new double[Npiemd];
lensespiemd.ellipticity = new double[Npiemd];
lensespiemd.ellipticity_potential = new double[Npiemd];
lensespiemd.rcore = 	new double[Npiemd];
lensespiemd.rcut = 	new double[Npiemd];
lensespiemd.z = 		new double[Npiemd];

lensessis.type = 	new int[Nsis];
lensessis.x  = 	new double[Nsis];
lensessis.y = 		new double[Nsis];
lensessis.b0 = 	new double[Nsis];
lensessis.ellipticity_angle = new double[Nsis];
lensessis.ellipticity = new double[Nsis];
lensessis.ellipticity_potential = new double[Nsis];
lensessis.rcore = 	new double[Nsis];
lensessis.rcut = 	new double[Nsis];
lensessis.z = 		new double[Nsis];

int i_sis(0), i_piemd(0);
int *i_point;

for (int i = 0; i <big; ++i){
	if (lens[i].type == 5){
		ilenses = &lensessis;
		i_point = &i_sis;
	}
	else{
		ilenses = &lensespiemd;
		i_point = &i_piemd;
	}
	//std::cout << *i_point << std::endl;
	ilenses->x[*i_point]  = 		lens[i].position.x;
	ilenses->y[*i_point] = 		lens[i].position.y;
	ilenses->b0[*i_point] = 		lens[i].b0;
	ilenses->ellipticity_angle[*i_point] = lens[i].ellipticity_angle;
	ilenses->ellipticity[*i_point] = lens[i].ellipticity;
	ilenses->ellipticity_potential[*i_point] = lens[i].ellipticity_potential;
	ilenses->rcore[*i_point] = 	lens[i].rcore;
	ilenses->rcut[*i_point] = 	lens[i].rcut;
	ilenses->z[*i_point] = 		lens[i].z;
	if (lens[i].type == 5){
		i_sis +=1;
	}
	else{
		i_piemd +=1;
	}

}

lenses[0] = lensessis;
lenses[1] = lensespiemd;

if(Nsis < small){
	Nlensessmall[0] = Nsis;
	Nlensessmall[1] = small-Nsis;
}
else{
	Nlensessmall[0] = small;
	Nlensessmall[1] = 0;
}

if(Nsis < medium){
	Nlensesmedium[0] = Nsis;
	Nlensesmedium[1] = medium-Nsis;
}
else{
	Nlensesmedium[0] = medium;
	Nlensesmedium[1] = 0;
}

if(Nsis < big){
	Nlensesbig[0] = Nsis;
	Nlensesbig[1] = big-Nsis;
}
else{
	Nlensesbig[0] = big;
	Nlensesbig[1] = 0;
}

std::cout  << Nlensessmall[0] << " " << Nlensessmall[1] << std::endl;
std::cout  << Nlensesmedium[0] << " " << Nlensesmedium[1] << std::endl;
std::cout  << Nlensesbig[0] << " " << Nlensesbig[1] << std::endl;

gettimeofday(&t1, 0);
module_potentialDerivatives_totalGradient(Nlensessmall,&image, lenses);
gettimeofday(&t2, 0);
module_potentialDerivatives_totalGradient(Nlensesmedium,&image, lenses);
gettimeofday(&t3, 0);
point grad = module_potentialDerivatives_totalGradient(Nlensesbig,&image, lenses);
gettimeofday(&t4, 0);

double time1 = (1000000.0*(t2.tv_sec-t1.tv_sec) + t2.tv_usec-t1.tv_usec)/1000000.0;
double time2 = (1000000.0*(t3.tv_sec-t2.tv_sec) + t3.tv_usec-t2.tv_usec)/1000000.0;
double time3 = (1000000.0*(t4.tv_sec-t3.tv_sec) + t4.tv_usec-t3.tv_usec)/1000000.0;

std::cout.precision(10);

std::cout << "Benchmark for Gradient SOA Calculation "<< std::endl;
std::cout << "Sample size " << small << ": " << time1 << std::endl;
std::cout << "Sample size " << medium << ": " << time2 << std::endl;
std::cout << "Sample size " << big << ": " << time3 << std::endl;
std::cout << "Sample Grad " << grad.x << " " << grad.y << std::endl;

std::ofstream myfile;
myfile.open ("BenchmarkGradSoA.txt");
myfile << "Benchmark for Gradient SOA Calculation "<< std::endl;
myfile << "Sample size " << small << ": " << time1 << std::endl;
myfile << "Sample size " << medium << ": " << time2 << std::endl;
myfile << "Sample size " << big << ": " << time3 << std::endl;
myfile.close();



}
