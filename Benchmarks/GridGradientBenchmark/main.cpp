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
#include <structure_hpc.h>
#include "timer.h"
#include "gradient.hpp"
#include "chi_CPU.hpp"
#include "module_cosmodistances.h"
#include "module_readParameters.hpp"
#ifdef __USE_GPU
#include "grid_gradient_GPU.cuh"
#endif

#ifdef __WITH_LENSTOOL
#include "setup.hpp"
#warning "linking with libtool..."
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
gradient_grid_GPU_sorted(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, int nbgridcells);

//
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
	struct Potential lenses[runmode.nhalos + runmode.npotfile-1];
	struct Potential_SOA lenses_SOA_table[NTYPES];
	struct Potential_SOA lenses_SOA;
	struct cline_param cline;
	struct potfile_param potfile;
	struct Potential potfilepotentials[runmode.npotfile];
	struct potentialoptimization host_potentialoptimization[runmode.nhalos];
	int nImagesSet[runmode.nsets]; // Contains the number of images in each set of images

	// This module function reads in the potential form and its parameters (e.g. NFW)
	// Input: input file
	// Output: Potentials and its parameters


	module_readParameters_Potential(inputFile, lenses, runmode.nhalos);
	//Converts to SOA
	//module_readParameters_PotentialSOA(inputFile, lenses, lenses_SOA, runmode.Nlens);
	module_readParameters_PotentialSOA(inputFile, lenses, &lenses_SOA, runmode.nhalos);
	//module_readParameters_PotentialSOA_nonsorted(inputFile, lenses, &lenses_SOA_nonsorted, runmode.nhalos);
	module_readParameters_debug_potential(runmode.debug, lenses, runmode.nhalos);
	//std::cerr << lenses_SOA[1].b0[0] << " " << lenses[0].b0  << std::endl;
	// This module function reads in the potfiles parameters
	// Input: input file
	// Output: Potentials from potfiles and its parameters

	if (runmode.potfile == 1 )
	{
		module_readParameters_readpotfiles_param(inputFile, &potfile);
		module_readParameters_debug_potfileparam(runmode.debug, &potfile);
		module_readParameters_readpotfiles(&runmode,&potfile,lenses);
		module_readParameters_debug_potential(runmode.debug, lenses, runmode.nhalos + runmode.npotfile);

	}
	//
	// This module function reads in the grid form and its parameters
	// Input: input file
	// Output: grid and its parameters
	//
	module_readParameters_Grid(inputFile, &frame);
	//
	if (runmode.image == 1 or runmode.inverse == 1 or runmode.time > 0)
	{

		// This module function reads in the strong lensing images
		module_readParameters_readImages(&runmode, images, nImagesSet);
		//runmode.nsets = runmode.nimagestot;
		for(int i = 0; i < runmode.nimagestot; ++i)
		{
			images[i].dls = module_cosmodistances_objectObject(lenses[0].z, images[i].redshift, cosmology);
			images[i].dos = module_cosmodistances_observerObject(images[i].redshift, cosmology);
			images[i].dr  = module_cosmodistances_lensSourceToObserverSource(lenses[0].z, images[i].redshift, cosmology);
		}
		module_readParameters_debug_image(runmode.debug, images, nImagesSet, runmode.nsets);
	}
	//
	if (runmode.inverse == 1)
	{
		// This module function reads in the potential optimisation limits
		module_readParameters_limit(inputFile, host_potentialoptimization, runmode.nhalos);
		module_readParameters_debug_limit(runmode.debug, host_potentialoptimization[0]);
	}
	//
	if (runmode.source == 1)
	{
		//Initialisation to default values.(Setting sources to z = 1.5 default value)
		for(int i = 0; i < runmode.nsets; ++i)
		{
			sources[i].redshift = 1.5;
		}
		// This module function reads in the strong lensing sources
		module_readParameters_readSources(&runmode, sources);
		//Calculating cosmoratios
		for(int i = 0; i < runmode.nsets; ++i)
		{
			sources[i].dls = module_cosmodistances_objectObject(lenses[0].z, sources[i].redshift, cosmology);
			sources[i].dos = module_cosmodistances_observerObject(sources[i].redshift, cosmology);
			sources[i].dr  = module_cosmodistances_lensSourceToObserverSource(lenses[0].z, sources[i].redshift, cosmology);
		}
		module_readParameters_debug_source(runmode.debug, sources, runmode.nsets);
	}
	//
	//
	//
	std::cout << "--------------------------" << std::endl << std::endl; fflush(stdout);

	double t_1,t_2,t_3,t_4;
	//
	//
	//
#ifdef __WITH_LENSTOOL
	printf("Setting up lenstool using %d lenses...", runmode.nhalos); fflush(stdout);
	convert_to_LT(&lenses_SOA, runmode.nhalos);	
	printf("ok\n");
#endif
	//
	// Lenstool-CPU Grid-Gradient
	//
	//Setting Test:
	double dx, dy;
	int grid_dim = runmode.nbgridcells;
	//
	dx = (frame.xmax - frame.xmin)/(runmode.nbgridcells-1);
	dy = (frame.ymax - frame.ymin)/(runmode.nbgridcells-1);
	//
	point test_point1_1, test_point2_2, test_result1_1, test_result2_2, test_pointN_N, test_resultN_N;
	double dlsds = images[0].dr;
	//
	test_point1_1.x = frame.xmin;
	test_point1_1.y = frame.ymin;
	test_point2_2.x = frame.xmin + dx;
	test_point2_2.y = frame.ymin + dy;
	test_pointN_N.x = frame.xmin + ((runmode.nbgridcells*runmode.nbgridcells-1)/runmode.nbgridcells)*dx;
	test_pointN_N.y = frame.ymin + ((runmode.nbgridcells*runmode.nbgridcells-1) % runmode.nbgridcells)*dy;

	//
	//
	//

#ifdef __WITH_LENSTOOL
        std::cout << " CPU Test Lenstool    ... ";
        struct point Grad;
        double *grid_grad_x, *grid_grad_y;
        grid_grad_x = (double *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(double));
        grid_grad_y = (double *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(double));
        double t_lt = -myseconds();
#pragma omp parallel for
//#pragma omp parallel for if (omp_get_num_threads() > 1) schedule(guided, 100)
        for (int jj = 0; jj < runmode.nbgridcells; ++jj)
                for (int ii = 0; ii < runmode.nbgridcells; ++ii)
                {
                        //  (index < grid_dim*grid_dim)

                        int index = jj*runmode.nbgridcells + ii;
                        //grid_grad_x[index] = 0.;
                        //grid_grad_y[index] = 0.;

                        struct point image_point;
                        image_point.x = frame.xmin + ii*dx;
                        image_point.y = frame.ymin + jj*dy;
#if 1
                        G.nlens = runmode.nhalos;
                        Grad = e_grad(&image_point);
                        grid_grad_x[index] = Grad.x;
                        grid_grad_y[index] = Grad.y;

#else
                        for (int lens = 0; lens < runmode.nhalos; ++lens)
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

	double* grid_gradient_x_cpu = (double *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(double));
        double* grid_gradient_y_cpu = (double *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(double));

	memset(grid_gradient_x_cpu, 0, (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(double));
	memset(grid_gradient_y_cpu, 0, (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(double));
	
	std::cout << " CPU Test lenstool_hpc... "; 
	//
	t_1 = -myseconds();
	//test_result1_1 = module_potentialDerivatives_totalGradient_SOA(&test_point1_1, &lenses_SOA, runmode.nhalos);
	//test_result2_2 = module_potentialDerivatives_totalGradient_SOA(&test_point2_2, &lenses_SOA, runmode.nhalos);
	//test_resultN_N = module_potentialDerivatives_totalGradient_SOA(&test_pointN_N, &lenses_SOA, runmode.nhalos);
	gradient_grid_CPU(grid_gradient_x_cpu, grid_gradient_y_cpu, &frame, &lenses_SOA, runmode.nhalos, grid_dim);	
	t_1 += myseconds();
	//
	std::cout << " Time = " << std::setprecision(15) << t_1 << std::endl;
	//
	//
	//
	double *grid_gradient_x, *grid_gradient_y;

#ifdef __USE_GPU
	// GPU test

	std::cout << " GPU Test... "; 

	double* grid_gradient_x_gpu = (double *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(double));
        double* grid_gradient_y_gpu = (double *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(double));
	//
	memset(grid_gradient_x_gpu, 0, (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(double));
	memset(grid_gradient_y_gpu, 0, (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(double));
	//
	grid_gradient_x_gpu = (double *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(double));
	grid_gradient_y_gpu = (double *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(double));

	//t_2 = -myseconds();
	//Packaging the image to sourceplane conversion
	//gradient_grid_CPU(grid_gradient_x,grid_gradient_y, &frame, &lenses_SOA, runmode.nhalos, grid_dim);
	//t_2 += myseconds();

	t_2 = -myseconds();
	//test();
	//test2();
	gradient_grid_GPU_sorted(grid_gradient_x_gpu, grid_gradient_y_gpu, &frame, &lenses_SOA, runmode.nhalos, grid_dim);
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

#ifdef __WITH_LENSTOOL 
	{
		double norm_x = 0.;
		double norm_y = 0.;
		//
		for (int ii = 0; ii < grid_dim*grid_dim; ++ii)
		{
			double g_x = grid_grad_x[ii];
			double g_y = grid_grad_y[ii];
			//
			double c_x = grid_gradient_x_cpu[ii];
			double c_y = grid_gradient_y_cpu[ii];
			//
			norm_x += (grid_grad_x[ii] - grid_gradient_x_cpu[ii])*(grid_grad_x[ii] - grid_gradient_x_cpu[ii]);
			norm_y += (grid_grad_y[ii] - grid_gradient_y_cpu[ii])*(grid_grad_y[ii] - grid_gradient_y_cpu[ii]);
		}
		//
		std::cout << "  l2 difference norm cpu = " << std::setprecision(15) << norm_x << " " << std::setprecision(15) << norm_y << std::endl;
	}
	//
#ifdef __USE_GPU
	{
		double norm_x = 0.;
		double norm_y = 0.;
		//
		for (int ii = 0; ii < grid_dim*grid_dim; ++ii)
		{
			double g_x = grid_grad_x[ii];
			double g_y = grid_grad_y[ii];
			//
			double c_x = grid_gradient_x_gpu[ii];
			double c_y = grid_gradient_y_gpu[ii];
			//
			norm_x += (grid_grad_x[ii] - grid_gradient_x_gpu[ii])*(grid_grad_x[ii] - grid_gradient_x_gpu[ii]);
			norm_y += (grid_grad_y[ii] - grid_gradient_y_gpu[ii])*(grid_grad_y[ii] - grid_gradient_y_gpu[ii]);
		}
		//
		std::cout << "  l2 difference norm gpu = " << std::setprecision(15) << norm_x << " " << std::setprecision(15) << norm_y << std::endl;
	}
#endif
#endif

#if 0
	t_2 = -myseconds();
	gradient_grid_GPU_sorted(grid_gradient_x,grid_gradient_y,&frame,&lenses_SOA,runmode.nhalos,grid_dim,0,grid_dim);
	t_2 += myseconds();

	std::cout << " gradient_grid_GPU_sorted Brute Force Benchmark " << std::endl;
	std::cout << " Test 1: " << std::endl;
	std::cout << " Point 1 : " << std::setprecision(5) << test_point1_1.x << " "<< test_point1_1.y <<  std::endl;
	std::cout << " Gradient  " << std::setprecision(5) << grid_gradient_x[0] << " "<< grid_gradient_y[0] <<  std::endl;
	std::cout << " Test 2: " << std::endl;
	std::cout << " Point 2 : " << std::setprecision(5) << test_point2_2.x << " "<< test_point2_2.y <<  std::endl;
	std::cout << " Gradient  " << std::setprecision(5) << grid_gradient_x[runmode.nbgridcells+1] << " "<< grid_gradient_y[runmode.nbgridcells+1] <<  std::endl;
	std::cout << " Test 3: " << std::endl;
	std::cout << " Point 3 : " << std::setprecision(5) << test_pointN_N.x << " "<< test_pointN_N.y <<  std::endl;
	std::cout << " Gradient  " << std::setprecision(5) << grid_gradient_x[runmode.nbgridcells*runmode.nbgridcells-1] << " "<< grid_gradient_y[runmode.nbgridcells*runmode.nbgridcells-1] <<  std::endl;
	std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;
#endif



}
