//
// This file is part of lenstoolhpc
// authors: gilles.fourestey@epfl.ch
// 
#include <math.h>
#include "module_cosmodistances.hpp"
#include "module_readParameters.hpp"

#include "allocation.hpp"

void PotentialSOAAllocation(Potential_SOA **lens_SOA, const int nhalos)
{
	
	Potential_SOA* p;
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
	p->lum 		 = (double*) malloc(sizeof(double)*nhalos);
	p->mag 		 = (double*) malloc(sizeof(double)*nhalos);
	p->dlsds 		 = (double*) malloc(sizeof(double)*nhalos);
	p->SOA_index 		 = (int*) malloc(sizeof(int)*nhalos);
	p->N_types 		 = (int*) malloc(sizeof(int)*100);
}



void PotentialSOADeallocation(Potential_SOA *lens_SOA)
{
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
	free(lens_SOA->lum);
	free(lens_SOA->mag);
	free(lens_SOA->dlsds);
	free(lens_SOA->SOA_index);
	free(lens_SOA->N_types);
	free(lens_SOA);
}
