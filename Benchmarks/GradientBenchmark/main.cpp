/**
 * @file   main.cpp
 * @Author Christoph Schaaefer, EPFL (christophernstrerne.schaefer@epfl.ch)
 * @date   October 2016
 * @brief  Benchmark for gradhalo function
 */

#include <iostream>
#include <iomanip>
#include <string.h>
//#include <cuda_runtime.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
//
//
#include <mm_malloc.h>

//#define __WITH_LENSTOOL 0
//
#ifdef __WITH_LENSTOOL
#warning "linking with libtool..."
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
//#include <liblenstool>
#endif
//#include "../../../Projects/lenstool-6.8.1/include/structure.h"
//#include "structure.h"
#include "structure_hpc.h"
#include "timer.h"
#include "gradient.hpp"
#include "gradient_avx.hpp"
#include "gradient_avx512f.hpp"
#include "setup.hpp"
#include <type.h>
//
#define NN 10000
//
//
//
#ifdef __WITH_LENSTOOL
struct g_mode   M;
struct g_pot    P[NPOTFILE];
struct g_pixel  imFrame, wFrame, ps, PSF;
struct g_cube   cubeFrame;
struct g_dyn    Dy;      //   //TV


struct g_source S;
struct g_image  I;
struct g_grille G;
struct g_msgrid H;  // multi-scale grid
struct g_frame  F;
struct g_large  L;
struct g_cosmo  C;
struct g_cline  CL;
struct g_observ O;
struct pot      lens[NLMAX];
struct pot      lmin[NLMAX], lmax[NLMAX], prec[NLMAX];
struct g_cosmo  clmin, clmax;       /*cosmological limits*/
struct galaxie  smin[NFMAX], smax[NFMAX];       // limits on source parameters
struct ipot     ip;
struct MCarlo   mc;
struct vfield   vf;
struct vfield   vfmin,vfmax; // limits on velocity field parameters
struct cline    cl[NIMAX];
lensdata *lens_table;

int  block[NLMAX][NPAMAX];      /*switch for the lens optimisation*/
int  cblock[NPAMAX];                /*switch for the cosmological optimisation*/
int  sblock[NFMAX][NPAMAX];                /*switch for the source parameters*/
int  vfblock[NPAMAX];                /*switch for the velocity field parameters*/
double excu[NLMAX][NPAMAX];
double excd[NLMAX][NPAMAX];
/* supplments tableaux de valeurs pour fonctions g pour Einasto
 *  * Ce sont trois variables globales qu'on pourra utiliser dans toutes les fonctions du projet
 *  */

#define CMAX 20
#define LMAX 80

float Tab1[LMAX][CMAX];
float Tab2[LMAX][CMAX];
float Tab3[LMAX][CMAX];


int      nrline, ntline, flagr, flagt;
long int  narclet;

struct point    gimage[NGGMAX][NGGMAX], gsource_global[NGGMAX][NGGMAX];
struct biline   radial[NMAX], tangent[NMAX];
struct galaxie  arclet[NAMAX], source[NFMAX], image[NFMAX][NIMAX];
struct galaxie  cimage[NFMAX];
struct pointgal     gianti[NPMAX][NIMAX];

struct point    SC;
double elix;
double alpha_e;

double *v_xx;
double *v_yy;
double **map_p;
double **tmp_p;
double **map_axx;
double **map_ayy;
#endif
//struct pot      lens_ref[NLMAX];
//
//#include <iitnotify.h>
int main()
{
	double t0, t1, t2, t3;
	//
	//Variable creation
	//
	point image;
	//
	//Initialisation
	//
        int nlenses;
        type_t x, y;
        type_t sol_grad_x, sol_grad_y;
	point grad; // store the result
	//pot* lens;
#ifdef __WITH_LENSTOOL
	setup_jauzac_LT(&nlenses, &image.x, &image.y, &sol_grad_x, &sol_grad_y);
	//
	struct point grad_lt, Grad;
	//printf("Number lenses = %d, b0 = %f\n", nlenses, lens[0].b0);
	//
	t0 = -myseconds();	
	for (int ii = 0; ii < NN; ++ii)
	{
		grad_lt.x = grad_lt.y = 0.;
		for (long int jj = 0; jj < nlenses; ++jj)
		{

			//printf("%f\n", lens[ii].b0);
			Grad = e_grad_pot(&image, jj);
			//
			grad_lt.x += Grad.x;
			grad_lt.y += Grad.y;
			//printf("%f %f\n", grad_lt.x, grad_lt.y);
		}
	}
	t0 += myseconds();
#endif
	//
	//
	//setup_jauzac(Potential** lens, int* nlenses, double* x, double* y, double* sol_grad_x, double* sol_grad_y)
	//
	//
	Potential*   lens_aos;
	setup_jauzac(&lens_aos, &nlenses, &image.x, &image.y, &sol_grad_x, &sol_grad_y);
	t1 = -myseconds();
	for (int ii = 0; ii < NN; ++ii)
		grad = module_potentialDerivatives_totalGradient(nlenses, &image, lens_aos);
	t1 += myseconds();
	//printf("---> grad = %f %f\n", grad.x, grad.y);
	//
	// Setting up the AOS potential 
	//
	point image_soa;
	int nlenses_soa;
	double x_soa, y_soa;
	type_t sol_grad_x_soa, sol_grad_y_soa;
	//
	Potential_SOA   lens_soa;
	setup_jauzac_SOA(&lens_soa, &nlenses_soa, &image_soa.x, &image_soa.y, &sol_grad_x_soa, &sol_grad_y_soa);
	//
	std::cout << "Benchmark for Gradient Calculation using  " << nlenses_soa << " lenses with type " << lens_soa.type[0] << ", image: " << image_soa.x << " " << image_soa.y << std::endl;
	// AVX version
	//
	point grad_soa_novec;
	double t21 = -myseconds();
        for (int ii = 0; ii < NN; ++ii)
        {
                //grad_soa = module_potentialDerivatives_totalGradient_SOA_AVX512(&image_soa, &lens_soa, nlenses_soa);
                //grad_soa = module_potentialDerivatives_totalGradient_8_SOA_AVX(&image_soa, &lens_soa, nlenses_soa);
                //
                //grad_soa_avx = module_potentialDerivatives_totalGradient_81_SOA_AVX(&image_soa, &lens_soa, 0, nlenses_soa)
;
                grad_soa_novec = module_potentialDerivatives_totalGradient_SOA_novec(&image_soa, &lens_soa, nlenses_soa);
                //grad_soa = module_potentialDerivatives_totalGradient_SOA(&image_soa, &lens_soa, 0, nlenses_soa);
        }
        t21 += myseconds();	
	//
	point grad_soa; // store the result
	t2 = -myseconds();
	for (int ii = 0; ii < NN; ++ii)
	{
		//grad_soa = module_potentialDerivatives_totalGradient_SOA_AVX512(&image_soa, &lens_soa, nlenses_soa);
		//grad_soa = module_potentialDerivatives_totalGradient_8_SOA_AVX(&image_soa, &lens_soa, nlenses_soa);
		//
		//grad_soa_avx = module_potentialDerivatives_totalGradient_81_SOA_AVX(&image_soa, &lens_soa, 0, nlenses_soa);
		grad_soa = module_potentialDerivatives_totalGradient_SOA(&image_soa, &lens_soa, nlenses_soa);
		//grad_soa = module_potentialDerivatives_totalGradient_SOA(&image_soa, &lens_soa, 0, nlenses_soa);
	}	
	t2 += myseconds();
	//__SSC_MARK(0x222);
	//
	// autovectorized version
	//
	point grad_soa_avx;
	t3 = -myseconds();
	for (int ii = 0; ii < NN; ++ii)
	{
		//grad_soa_avx = module_potentialDerivatives_totalGradient_SOA_AVX512(&image_soa, &lens_soa, nlenses_soa);
		//                grad_soa = module_potentialDerivatives_totalGradient_81_SOA_AVX(&image_soa, &lens_soa, nlenses_soa);
		//                                //grad_soa = module_potentialDerivatives_totalGradient_81_SOA_AVX(&image_soa, &lens_soa, nlenses_soa);
		//grad_soa = module_potentialDerivatives_totalGradient_81_SOA(&image_soa, &lens_soa, 0, nlenses_soa);
		//grad_soa = module_potentialDerivatives_totalGradient_8_SOA(&image_soa, &lens_soa, 0, nlenses_soa);
#ifdef 0 //_double
		grad_soa_avx = module_potentialDerivatives_totalGradient_SOA_AVX(&image_soa, &lens_soa, nlenses_soa);
#endif
	}
	t3 += myseconds();
	//
#ifdef _double
	std::cout << " Double calculation   = "  << std::endl;
#else
	std::cout << " Float  calculation   = "  << std::endl;
#endif

	std::cout << " ref sol   = " << std::setprecision(15) << sol_grad_x << " " << std::setprecision(15) << sol_grad_y << std::endl;
	//
#ifdef __WITH_LENSTOOL
	std::cout << " Lenstool sol   = " << std::setprecision(15) << grad_lt.x << " " << std::setprecision(15) << grad_lt.y << ", time = " << t0 << " s." << std::endl;
#endif
	//
	std::cout << " grad           = " << std::setprecision(15) << grad.x << " " << std::setprecision(15) << grad.y << ", time = " << t1 << " s., speedup = " << (double) t0/t1 << std::endl;
	//
	std::cout << " grad novec     = " << std::setprecision(15) << grad.x << " " << std::setprecision(15) << grad.y << ", time = " << t21 << " s., speedup = " << (double) t0/t21 << std::endl;
	//
	std::cout << " grad SIMD      = " << std::setprecision(15) << grad_soa.x << " " << std::setprecision(15) << grad_soa.y << ", time = " << t2 << " s., speedup = " << (double) t0/t2 << std::endl;
	//
	std::cout << " grad handcoded = " << std::setprecision(15) << grad_soa_avx.x << " " << std::setprecision(15) << grad_soa_avx.y << ", time = " << t3 << " s. speedup = " << (double) t0/t3 << std::endl;
	//
}

