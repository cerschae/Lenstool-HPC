/**
* @file   main.cpp
* @Author Christoph Schaaefer, EPFL (christophernstrerne.schaefer@epfl.ch)
* @date   October 2016
* @brief  Benchmark for gradhalo function
*/

#include <iostream>
#include <iomanip>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>
//
//#include <mm_malloc.h>
#include <omp.h>
//
//#include <cuda_runtime.h>
#include <structure_hpc.hpp>
#include "timer.h"
#include "gradient.hpp"
#include "chi_CPU.hpp"
#include "module_cosmodistances.hpp"
#include "module_readParameters.hpp"
#ifdef __WITH_GPU
#warning "GPU support enabled"
#include "grid_gradient_GPU.cuh"
#endif

#ifdef __WITH_LENSTOOL
#include "setup.hpp"
#warning "linking with lenstool..."
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
//
//
struct g_mode   M;
struct g_pot    P[NPOTFILE];
struct g_pixel  imFrame, wFrame, ps, PSF;
struct g_cube   cubeFrame;
struct g_dyn    Dy;      //   //TV
//
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
//
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

void
gradient_grid_GPU_sorted(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, int nbgridcells);
//
//
int module_readCheckInput_readInput(int argc, char *argv[])
{
	/// check if there is a correct number of arguments, and store the name of the input file in infile

	char* infile;
	struct stat  file_stat;

	// If we do not have 3 arguments, stop
	if ( argc != 3 )
	{
		fprintf(stderr, "\nUnexpected number of arguments\n");
		fprintf(stderr, "\nUSAGE:\n");
		fprintf(stderr, "lenstool  input_file  output_directorypath  [-n]\n\n");
		exit(-1);
	}
	else if ( argc == 3 )
		infile=argv[1];
	std::ifstream ifile(infile,std::ifstream::in); // Open the file


	int ts = (int) time (NULL);
	char buffer[10];
	std::stringstream ss;
	ss << ts;
	std::string trimstamp = ss.str();
	//
	std::string outdir = argv[2];
	outdir += "-";
	outdir += trimstamp;
	std::cout << outdir << std::endl;

	// check whether the output directory already exists
	if (stat(outdir.c_str(), &file_stat) < 0){
		mkdir(outdir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH );
	}
	else {
		printf("Error : Directory %s already exists. Specify a non existing directory.\n",argv[2]);
		exit(-1);
	}

	// check whether the input file exists. If it could not be opened (ifile = 0), it does not exist
	if(ifile){
		ifile.close();
	}
	else{
		printf("The file %s does not exist, please specify a valid file name\n",infile);
		exit(-1);
	}

	return 0;
}
//
//
//
int main(int argc, char *argv[])
{
	//
	// Setting Up the problem
	//

	// This module function reads the terminal input when calling LENSTOOL and checks that it is correct
	// Otherwise it exits LENSTOOL
	// 
	char cwd[1024];
	if (getcwd(cwd, sizeof(cwd)) != NULL)
		fprintf(stdout, "Current working dir: %s\n", cwd);
	//
	module_readCheckInput_readInput(argc, argv);
	//
#if 1
	// This module function reads the cosmology parameters from the parameter file
	// Input: struct cosmologicalparameters cosmology, parameter file
	// Output: Initialized cosmology struct
	cosmo_param cosmology;  // Cosmology struct to store the cosmology data from the file
	std::string inputFile = argv[1];   // Input file
	module_readParameters_readCosmology(inputFile, cosmology);
	//
	// This module function reads the runmode paragraph and the number of sources, arclets, etc. in the parameter file.
	// The runmode_param stores the information of what exactly the user wants to do with lenstool.
	struct runmode_param runmode;
	module_readParameters_readRunmode(inputFile, &runmode);
	module_readParameters_debug_cosmology(runmode.debug, cosmology);
	module_readParameters_debug_runmode(runmode.debug, runmode);
	//
	//=== Declaring variables
	//
	struct grid_param frame;
	struct galaxy images[runmode.nimagestot];
	struct galaxy sources[runmode.nsets];
	//struct Potential lenses[runmode.nhalos + runmode.npotfile-1];
	struct Potential_SOA lenses_SOA_table[NTYPES];
	struct Potential_SOA lenses_SOA;
	struct cline_param cline;
	struct potfile_param potfile;
	//struct Potential potfilepotentials[runmode.npotfile];
	struct potentialoptimization host_potentialoptimization[runmode.nhalos];
	int nImagesSet[runmode.nsets]; // Contains the number of images in each set of images

	// This module function reads in the potential form and its parameters (e.g. NFW)
	// Input: input file
	// Output: Potentials and its parameters


	module_readParameters_PotentialSOA_direct(inputFile, &lenses_SOA, runmode.nhalos, runmode.n_tot_halos, cosmology);
	module_readParameters_debug_potential_SOA(1, lenses_SOA, runmode.nhalos);

	//module_readParameters_Potential(inputFile, lenses, runmode.nhalos);
	//Converts to SOA
	//module_readParameters_PotentialSOA(inputFile, lenses, &lenses_SOA, runmode.nhalos);
	//module_readParameters_debug_potential(runmode.debug, lenses, runmode.nhalos);
	// This module function reads in the potfiles parameters
	// Input: input file
	// Output: Potentials from potfiles and its parameters

	if (runmode.potfile == 1 )
	{
		module_readParameters_readpotfiles_param(inputFile, &potfile, cosmology);
		module_readParameters_debug_potfileparam(1, &potfile);
		module_readParameters_readpotfiles_SOA(&runmode, &cosmology,&potfile,&lenses_SOA);
		module_readParameters_debug_potential_SOA(1, lenses_SOA, runmode.n_tot_halos);

	}
	//
	// This module function reads in the grid form and its parameters
	// Input: input file
	// Output: grid and its parameters
	//
	module_readParameters_Grid(inputFile, &frame);
	//

	//
	//
	//
	std::cout << "--------------------------" << std::endl << std::endl; fflush(stdout);

	double t_1,t_2,t_3,t_4;
	//
	//
	//
#ifdef __WITH_LENSTOOL
	printf("Setting up lenstool using %d lenses...", runmode.n_tot_halos); fflush(stdout);
	convert_to_LT(&lenses_SOA, runmode.n_tot_halos);
	printf("ok\n");
#endif
	//
	// Lenstool-CPU Grid-Gradient
	//


	//Setting Test:
	type_t dx, dy;
	int grid_dim = runmode.nbgridcells;
	//
	dx = (frame.xmax - frame.xmin)/(runmode.nbgridcells-1);
	dy = (frame.ymax - frame.ymin)/(runmode.nbgridcells-1);
	//
	//

#ifdef __WITH_LENSTOOL
        std::cout << " CPU Test Lenstool    ... ";
        struct point Grad;
        double *grid_grad_x, *grid_grad_y;
        grid_grad_x = (double *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
        grid_grad_y = (double *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
        double t_lt = -myseconds();
#pragma omp parallel for
        for (int jj = 0; jj < runmode.nbgridcells; ++jj)
                for (int ii = 0; ii < runmode.nbgridcells; ++ii)
                {
                        //  (index < grid_dim*grid_dim)

                        int index = jj*runmode.nbgridcells + ii;

                        struct point image_point;
                        image_point.x = frame.xmin + ii*dx;
                        image_point.y = frame.ymin + jj*dy;
#if 1
                        G.nlens = runmode.n_tot_halos;
                        Grad = e_grad(&image_point);
                        grid_grad_x[index] = Grad.x;
                        grid_grad_y[index] = Grad.y;

#else
                        for (int lens = 0; lens < runmode.n_tot_halos; ++lens)
                        {

                                struct point Grad = e_grad_pot(&image_point, lens);
                                //printf("%f %f\n", Grad.x, Grad.y);
                                //
                                grid_grad_x[index] += Grad.x;
                                grid_grad_y[index] += Grad.y;
                        }
#endif
		}
	t_lt += myseconds();
	std::cout << " Time = " << t_lt << " s." << std::endl;
#endif

	//

	type_t* grid_gradient_x_cpu = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
    type_t* grid_gradient_y_cpu = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));

	memset(grid_gradient_x_cpu, 0, (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	memset(grid_gradient_y_cpu, 0, (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	
	std::cout << " CPU Test lenstool_hpc... "; 
	//
	t_1 = -myseconds();
	int Nstat = 1;
	//for(int ii = 0; ii < Nstat; ++ii) {
		gradient_grid_CPU(grid_gradient_x_cpu, grid_gradient_y_cpu, &frame, &lenses_SOA, runmode.n_tot_halos, runmode.nbgridcells);
		//gradient_grid_CPU_print(grid_gradient_x_cpu, grid_gradient_y_cpu, &frame, &lenses_SOA, runmode.nhalos, grid_dim);
	//	}
	t_1 += myseconds();
	//
	std::cout << " Time = " << std::setprecision(15) << t_1 << std::endl;
	//
	//
	//
	type_t *grid_gradient_x, *grid_gradient_y;

#ifdef __WITH_GPU
#warning "using GPUs..."
	// GPU test

	std::cout << " GPU Test... "; 

	type_t* grid_gradient_x_gpu = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
       type_t* grid_gradient_y_gpu = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	//
	memset(grid_gradient_x_gpu, 0, (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	memset(grid_gradient_y_gpu, 0, (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	//
	grid_gradient_x_gpu = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	grid_gradient_y_gpu = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));

	//t_2 = -myseconds();
	//Packaging the image to sourceplane conversion
	//gradient_grid_CPU(grid_gradient_x,grid_gradient_y, &frame, &lenses_SOA, runmode.nhalos, grid_dim);
	//t_2 += myseconds();
	//Some Sort of cache or initial overhead problem... alway takes 0.2 sec the first time
	gradient_grid_GPU(grid_gradient_x_gpu, grid_gradient_y_gpu, &frame, &lenses_SOA, runmode.nhalos, runmode.nbgridcells);

	t_2 = -myseconds();
	//test();
	//test2();
	//for(int ii = 0; ii < Nstat; ++ii) {
	gradient_grid_GPU(grid_gradient_x_gpu, grid_gradient_y_gpu, &frame, &lenses_SOA, runmode.n_tot_halos, runmode.nbgridcells);
	//}
	//module_potentialDerivatives_totalGradient_SOA_CPU_GPU(grid_gradient_x_gpu, grid_gradient_y_gpu, &frame, &lenses_SOA, runmode.nhalos, grid_dim);
	//gradient_gid_CPU(grid_gradient_x, grid_gradient_y, &frame, &lenses_SOA, runmode.nhalos, grid_dim);
	t_2 += myseconds();

	std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;

/*
	std::cout << " gradient_grid_CPU Brute Force Benchmark " << std::endl;
	std::cout << " Test 1: " << std::endl;
	std::cout << " Point 1 : " << std::setprecision(5) << test_point1_1.x << " "<< test_point1_1.y <<  std::endl;
	std::cout << " Gradient  " << std::setprecision(5) << grid_gradient_x_gpu[0] << " "<< grid_gradient_y[0] <<  std::endl;
	std::cout << " Test 2: " << std::endl;
	std::cout << " Point 2 : " << std::setprecision(5) << test_point2_2.x << " "<< test_point2_2.y <<  std::endl;
	std::cout << " Gradient  " << std::setprecision(5) << grid_gradient_x_gpu[runmode.nbgridcells+1] << " "<< grid_gradient_y[runmode.nbgridcells+1] <<  std::endl;
	std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;
*/
#endif
	std::ofstream myfile;
#ifdef __WITH_LENSTOOL 
	{
		type_t norm_x = 0.;
		type_t norm_y = 0.;

		type_t sum_x = 0.;
		type_t sum_y = 0.;
		//
		for (int ii = 0; ii < grid_dim*grid_dim; ++ii)
		{
			type_t g_x = grid_grad_x[ii];
			type_t g_y = grid_grad_y[ii];
			sum_x += grid_grad_x[ii]*grid_grad_x[ii];
			sum_y += grid_grad_y[ii]*grid_grad_y[ii];

			//
			type_t c_x = grid_gradient_x_cpu[ii];
			type_t c_y = grid_gradient_y_cpu[ii];
			//
			norm_x += (grid_grad_x[ii] - grid_gradient_x_cpu[ii])*(grid_grad_x[ii] - grid_gradient_x_cpu[ii]);
			norm_y += (grid_grad_y[ii] - grid_gradient_y_cpu[ii])*(grid_grad_y[ii] - grid_gradient_y_cpu[ii]);
		}
		//
		std::cout << "  l2 difference norm cpu = " << std::setprecision(15) << norm_x << " " << std::setprecision(15) << norm_y << std::endl;
		//std::cout << sum_x << " " << std::setprecision(15) << sum_y << std::setprecision(15) << std::endl;
#if 0

		myfile.open ("lenstool_grid_x.txt");

		for (int ii = 0; ii < grid_dim*grid_dim; ++ii)
		{
			myfile << ii << " " << grid_grad_x[ii]<< std::setprecision(15)  << " " << std::endl;
		}
		myfile.close();

		myfile.open ("lenstool_grid_y.txt");
		for (int ii = 0; ii < grid_dim*grid_dim; ++ii)
		{
			myfile << ii << " " << grid_grad_y[ii]<< std::setprecision(15)  << " " << std::endl;
		}
		myfile.close();
#endif

	}
	//
#ifdef __WITH_GPU
	{
		type_t norm_x = 0.;
		type_t norm_y = 0.;

		//
		for (int ii = 0; ii < grid_dim*grid_dim; ++ii)
		{
			type_t g_x = grid_grad_x[ii];
			type_t g_y = grid_grad_y[ii];
			//
			type_t c_x = grid_gradient_x_gpu[ii];
			type_t c_y = grid_gradient_y_gpu[ii];
			//
			norm_x += (grid_grad_x[ii] - grid_gradient_x_gpu[ii])*(grid_grad_x[ii] - grid_gradient_x_gpu[ii]);
			norm_y += (grid_grad_y[ii] - grid_gradient_y_gpu[ii])*(grid_grad_y[ii] - grid_gradient_y_gpu[ii]);

		}
		//
		std::cout << "  l2 difference norm gpu = " << std::setprecision(15) << norm_x << " " << std::setprecision(15) << norm_y << std::endl;
	}
#endif
#endif

#ifdef __WITH_GPU
	{
		type_t norm_x = 0.;
		type_t norm_y = 0.;

		type_t sum_x_cpu = 0.;
		type_t sum_y_cpu = 0.;
		type_t sum_x_gpu = 0.;
		type_t sum_y_gpu = 0.;
		//
		for (int ii = 0; ii < grid_dim*grid_dim; ++ii)
		{
			//
			sum_x_cpu += grid_gradient_x_cpu[ii]*grid_gradient_x_cpu[ii];
			sum_y_cpu += grid_gradient_y_cpu[ii]*grid_gradient_y_cpu[ii];
			sum_x_gpu += grid_gradient_x_gpu[ii]*grid_gradient_x_gpu[ii];
			sum_y_gpu += grid_gradient_y_gpu[ii]*grid_gradient_y_gpu[ii];

			norm_x += (grid_gradient_x_cpu[ii] - grid_gradient_x_gpu[ii])*(grid_gradient_x_cpu[ii] - grid_gradient_x_gpu[ii]);
			norm_y += (grid_gradient_y_cpu[ii] - grid_gradient_y_gpu[ii])*(grid_gradient_y_cpu[ii] - grid_gradient_y_gpu[ii]);
		}

		sum_x_cpu -= 4761763143.24101;
		sum_y_cpu -= 5412618205.81843;
		sum_x_gpu -= 4761763143.24101;
		sum_y_gpu -= 5412618205.81843;

		std::cout << "  l2 difference norm cpu-gpu = " << std::setprecision(15) << norm_x << " " << std::setprecision(15) << norm_y << std::endl;
		std::cout << "  sum x cpu = " << std::setprecision(15) << sum_x_cpu << " sum_y_cpu  " << std::setprecision(15) << sum_y_cpu << std::endl;
		std::cout << "  sum x gpu = " << std::setprecision(15) << sum_x_gpu << " sum_y_gpu  " << std::setprecision(15) << sum_y_gpu << std::endl;

	}
#endif
	std::cout << "Exiting..." << std::endl;
#endif

}
