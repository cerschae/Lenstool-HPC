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
//#include <cuda.h>
#include "timer.h"
#include "gradient.hpp"
#include "chi_CPU.hpp"
#include "module_cosmodistances.hpp"
#include "module_readParameters.hpp"
#include "grid_gradient2_CPU.hpp"
#ifdef __WITH_GPU
#include "grid_gradient_GPU.cuh"
#include "grid_gradient2_GPU.cuh"
//#include "gradient2_GPU.cuh"
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
	struct Potential_SOA lenses_SOA_table[NTYPES];
	struct Potential_SOA lenses_SOA;
	struct cline_param cline;
	struct potfile_param potfile;
	struct potentialoptimization host_potentialoptimization[runmode.nhalos];
	int nImagesSet[runmode.nsets]; // Contains the number of images in each set of images

	// This module function reads in the potential form and its parameters (e.g. NFW)
	// Input: input file
	// Output: Potentials and its parameters


	module_readParameters_PotentialSOA_direct(inputFile, &lenses_SOA, runmode.nhalos, runmode.n_tot_halos, cosmology);
	module_readParameters_debug_potential_SOA(1, lenses_SOA, runmode.nhalos);

	//std::cerr <<"b0: "<< lenses_SOA.b0[0] << std::endl;

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

#include "gradient2.hpp"

	//Setting Test:
	type_t dx, dy;
	int grid_dim = runmode.nbgridcells;
	//
	dx = (frame.xmax - frame.xmin)/(runmode.nbgridcells-1);
	dy = (frame.ymax - frame.ymin)/(runmode.nbgridcells-1);
	//
	//
	//
#ifdef __WITH_LENSTOOL
        std::cout << " CPU Test Lenstool    ... ";
        matrix *grid_grad2;
        grid_grad2 = (matrix *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(matrix));
        double t_lt = -myseconds();
#pragma omp parallel for
//#pragma omp parallel for if (omp_get_num_threads() > 1) schedule(guided, 100)
        for (int jj = 0; jj < runmode.nbgridcells; ++jj){
                for (int ii = 0; ii < runmode.nbgridcells; ++ii)
                {
                        //  (index < grid_dim*grid_dim)
                    struct matrix Grad2;

                        int index = jj*runmode.nbgridcells + ii;

                        grid_grad2[index].a = 0.;
                        grid_grad2[index].b = 0.;
                        grid_grad2[index].c = 0.;
                        grid_grad2[index].d = 0.;

                        struct point image_point;
                        image_point.x = frame.xmin + ii*dx;
                        image_point.y = frame.ymin + jj*dy;
                        //std::cerr << "*" <<  image_point.x <<  " "<<image_point.y << std::endl;

                        for (int lens = 0; lens < runmode.n_tot_halos; ++lens)
                        {

                        		Grad2 = e_grad2_pot(&image_point, lens);

                                //

                                grid_grad2[index].a += Grad2.a;
                                grid_grad2[index].b += Grad2.b;
                                grid_grad2[index].c += Grad2.c;
                                grid_grad2[index].d += Grad2.d;
                                //std::cerr << "**" << Grad2.a << Grad2.b << Grad2.c << Grad2.d << std::endl;

                        }
                }
		}

	t_lt += myseconds();
	std::cout << " Time = " << t_lt << " s." << std::endl;
#endif

	//
	matrix* grid_gradient2_cpu = (matrix *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(matrix));
	memset(grid_gradient2_cpu, 0, (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(matrix));
	
	std::cout << " CPU Test lenstool_hpc... "; 
	//
	t_1 = -myseconds();

	int Nstat = 1;
	for(int ii = 0; ii < Nstat; ++ii) {
		gradient2_grid_CPU(grid_gradient2_cpu, &frame, &lenses_SOA, runmode.n_tot_halos, grid_dim);
		}
	t_1 += myseconds();
	//
	std::cout << " Time = " << std::setprecision(15) << t_1 << std::endl;
#if 1
#ifdef __WITH_GPU
	// GPU test

	std::cout << " GPU Test... "; 

	type_t* grid_gradient_a_gpu = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
    type_t* grid_gradient_b_gpu = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	type_t* grid_gradient_c_gpu = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
    type_t* grid_gradient_d_gpu = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	//
	memset(grid_gradient_a_gpu, 0, (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	memset(grid_gradient_b_gpu, 0, (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	memset(grid_gradient_c_gpu, 0, (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	memset(grid_gradient_d_gpu, 0, (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	//
	grid_gradient_a_gpu = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	grid_gradient_b_gpu = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	grid_gradient_c_gpu = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	grid_gradient_d_gpu = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));

	gradient2_grid_GPU(grid_gradient_a_gpu, grid_gradient_b_gpu,grid_gradient_c_gpu, grid_gradient_d_gpu, &frame, &lenses_SOA, runmode.n_tot_halos, grid_dim);

	t_2 = -myseconds();
	for(int ii = 0; ii < Nstat; ++ii) {
		gradient2_grid_GPU(grid_gradient_a_gpu, grid_gradient_b_gpu,grid_gradient_c_gpu, grid_gradient_d_gpu, &frame, &lenses_SOA, runmode.n_tot_halos, grid_dim);
	}
	t_2 += myseconds();
	std::cerr << "**" << grid_gradient_a_gpu[0] << grid_gradient_b_gpu[0] << grid_gradient_c_gpu[0] << grid_gradient_d_gpu[0] << std::endl;
	std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;

#endif
	std::ofstream myfile;
#ifdef __WITH_LENSTOOL 
	{
		type_t norm_a = 0.;
		type_t norm_b = 0.;
		type_t norm_c = 0.;
		type_t norm_d = 0.;
		//
		for (int ii = 0; ii < grid_dim*grid_dim; ++ii)
		{
			//
			norm_a += (grid_grad2[ii].a - grid_gradient2_cpu[ii].a)*(grid_grad2[ii].a - grid_gradient2_cpu[ii].a);
			norm_b += (grid_grad2[ii].b - grid_gradient2_cpu[ii].b)*(grid_grad2[ii].b - grid_gradient2_cpu[ii].b);
			norm_c += (grid_grad2[ii].c - grid_gradient2_cpu[ii].c)*(grid_grad2[ii].c - grid_gradient2_cpu[ii].c);
			norm_d += (grid_grad2[ii].d - grid_gradient2_cpu[ii].d)*(grid_grad2[ii].d - grid_gradient2_cpu[ii].d);
		}
		//
		std::cout << "  l2 difference norm cpu = " << std::setprecision(15) << norm_a << " " << std::setprecision(15) << norm_b << " " << std::setprecision(15) << norm_c << " " << std::setprecision(15) << norm_d << std::endl;

	}
	//
#if 1
#ifdef __WITH_GPU
	{
		type_t norm_a = 0.;
		type_t norm_b = 0.;
		type_t norm_c = 0.;
		type_t norm_d = 0.;
		//
		for (int ii = 0; ii < grid_dim*grid_dim; ++ii)
		{
			//
			norm_a += (grid_grad2[ii].a - grid_gradient_a_gpu[ii])*(grid_grad2[ii].a - grid_gradient_a_gpu[ii]);
			norm_b += (grid_grad2[ii].b - grid_gradient_b_gpu[ii])*(grid_grad2[ii].b - grid_gradient_b_gpu[ii]);
			norm_c += (grid_grad2[ii].c - grid_gradient_c_gpu[ii])*(grid_grad2[ii].c - grid_gradient_c_gpu[ii]);
			norm_d += (grid_grad2[ii].d - grid_gradient_d_gpu[ii])*(grid_grad2[ii].d - grid_gradient_d_gpu[ii]);
		}
		//
		std::cout << "  l2 difference norm gpu = " << std::setprecision(15) << norm_a << " " << std::setprecision(15) << norm_b << " " << std::setprecision(15) << norm_c << " " << std::setprecision(15) << norm_d << std::endl;
	}
#endif
#endif
#endif

#endif
}
