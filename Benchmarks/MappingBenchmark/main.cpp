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
#include "grid_amplif_CPU.hpp"
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


	module_readParameters_PotentialSOA_direct(inputFile, &lenses_SOA, runmode.nhalos, runmode.npotfile, cosmology);
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
		module_readParameters_readpotfiles_SOA(&runmode,&potfile,&lenses_SOA);
		module_readParameters_debug_potential_SOA(1, lenses_SOA, runmode.nhalos + runmode.npotfile);

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

#include "gradient2.hpp"

	//Setting Test:
	type_t dx, dy;
	int grid_dim = runmode.nbgridcells;
	//
	dx = (frame.xmax - frame.xmin)/(runmode.nbgridcells-1);
	dy = (frame.ymax - frame.ymin)/(runmode.nbgridcells-1);
#if 0
	//
	//mdci05(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0, &g05c);
	//
	point test_point1_1, test_point2_2, test_result1_1, test_result2_2, test_pointN_N, test_resultN_N;
	type_t dlsds = images[0].dr;
	//
	test_point1_1.x = frame.xmin;
	test_point1_1.y = frame.ymin;
	test_point2_2.x = frame.xmin + runmode.nbgridcells/2*dx;
	test_point2_2.y = frame.ymin + runmode.nbgridcells/2*dy;
	test_pointN_N.x = frame.xmin + ((runmode.nbgridcells*runmode.nbgridcells-1)/runmode.nbgridcells)*dx;
	test_pointN_N.y = frame.ymin + ((runmode.nbgridcells*runmode.nbgridcells-1) % runmode.nbgridcells)*dy;


	matrix test1, test_int, test2;
	test1.a  = 0;
	test1.b  = 0;
	test1.c  = 0;
	test1.d  = 0;
    for (int lens = 0; lens < runmode.nhalos; ++lens)
    {
    	//std::cerr <<"Point: " << test_point2_2.x << " " << test_point2_2.y << std::endl;
    	test_int = e_grad2_pot(&test_point2_2,  lens);
    	test1.a += test_int.a;
    	test1.b += test_int.b;
    	test1.c += test_int.c;
    	test1.d += test_int.d;
    	//std::cerr <<"**" << test1.a << " " << test_int.a << std::endl;
    }
	test2 = module_potentialDerivatives_totalGradient2_SOA(&test_point2_2, &lenses_SOA, runmode.nhalos);
	std::cerr << test1.a << " " << test2.a << std::endl;
	std::cerr << test1.b << " " << test2.b << std::endl;
	std::cerr << test1.c << " " << test2.c << std::endl;
	std::cerr << test1.d << " " << test2.d << std::endl;

#endif

	//
	//
	//

#ifdef __WITH_LENSTOOL
        std::cout << " CPU Test Lenstool    ... ";
        type_t *ampli;
        ampli = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
        double t_lt = -myseconds();
        type_t iamp = 5;
//#pragma omp parallel for if (omp_get_num_threads() > 1) schedule(guided, 100)
//#pragma omp parallel for
        for (int jj = 0; jj < runmode.nbgridcells; ++jj){
                for (int ii = 0; ii < runmode.nbgridcells; ++ii)
                {
                        //  (index < grid_dim*grid_dim)
                    struct matrix Grad2, Grad2_temp; type_t kappa, ga1, ga2, gam;

                        int index = jj*runmode.nbgridcells + ii;

                        ampli[index] = 0.;
                        Grad2.a = Grad2.b = Grad2.c = Grad2.d = 0 ;
                        struct point image_point;
                        image_point.x = frame.xmin + ii*dx;
                        image_point.y = frame.ymin + jj*dy;
                        //std::cerr << "*" <<  image_point.x <<  " "<<image_point.y << std::endl;

                        for (int lens = 0; lens < runmode.nhalos; ++lens)
                        {

                        		Grad2_temp = e_grad2_pot(&image_point, lens);
                        		Grad2.a += Grad2_temp.a;
                        		Grad2.b += Grad2_temp.b;
                        		Grad2.c += Grad2_temp.c;
                        		Grad2.d += Grad2_temp.d;
                        }
                		//
                		Grad2.a /= 1.0001;
                		Grad2.b /= 1.0001;
                		Grad2.c /= 1.0001;

                        kappa = (Grad2.a + Grad2.c) / 2.;
                        ga1 = (Grad2.a - Grad2.c) / 2.;
                        ga2 = Grad2.b;
                        gam = sqrt(ga1 * ga1 + ga2 * ga2);
                        if (iamp == 5)
                            ampli[index] = kappa;
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
    type_t *ampli_CPU;
    ampli_CPU = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	int Nstat = 1;
	for(int ii = 0; ii < Nstat; ++ii) {
		amplif_grid_CPU(ampli_CPU, &frame, &lenses_SOA, runmode.nhalos, grid_dim, 5);
		}
	t_1 += myseconds();
	//
	std::cout << " Time = " << std::setprecision(15) << t_1 << std::endl;
#if 0
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

	t_2 = -myseconds();
	for(int ii = 0; ii < Nstat; ++ii) {
		gradient2_grid_GPU(grid_gradient_a_gpu, grid_gradient_b_gpu,grid_gradient_c_gpu, grid_gradient_d_gpu, &frame, &lenses_SOA, runmode.nhalos, grid_dim);
	}
	t_2 += myseconds();
	std::cerr << "**" << grid_gradient_a_gpu[0] << grid_gradient_b_gpu[0] << grid_gradient_c_gpu[0] << grid_gradient_d_gpu[0] << std::endl;
	std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;

#endif
#endif
	std::ofstream myfile;
#ifdef __WITH_LENSTOOL 
	{
		type_t norm_a = 0.;
		//
		for (int ii = 0; ii < grid_dim*grid_dim; ++ii)
		{
			//
			norm_a += (ampli[ii] - ampli_CPU[ii])*(ampli[ii] - ampli_CPU[ii]);
		}
		//
		std::cout << "  l2 difference norm cpu = " << std::setprecision(15) << norm_a << std::endl;
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
#endif


	//
#if 0
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
		std::cout << "  l2 difference norm cpu = " << std::setprecision(15) << norm_a << " " << std::setprecision(15) << norm_b << " " << std::setprecision(15) << norm_c << " " << std::setprecision(15) << norm_d << std::endl;
	}
#endif
#endif


#if 0
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

#if 0
		std::cout << " CPUTEST " << std::endl;
		myfile.open ("Float_x_R_XMIN0.txt");
		for (int ii = 0; ii < grid_dim*grid_dim; ++ii)
		{
			myfile << ii << " " << grid_gradient_x_cpu[ii]<< std::setprecision(15)  << " " << std::endl;
		}
		myfile.close();
		myfile.open ("Float_y_R_XMIN0.txt");
		for (int ii = 0; ii < grid_dim*grid_dim; ++ii)
		{
			myfile << ii << " " << grid_gradient_y_cpu[ii]<< std::setprecision(15)  << " " << std::endl;
		}
		myfile.close();
#endif

#if 0
		std::cout << " GPUTEST " << std::endl;
		myfile.open ("Float_true_coord.x_2.txt");
		for (int ii = 0; ii < grid_dim*grid_dim; ++ii)
		{
			myfile << ii << " " << grid_gradient_x_gpu[ii]<< std::setprecision(15)  << " " << std::endl;
		}
		myfile.close();
		myfile.open ("Float_x_4.txt");
		for (int ii = 0; ii < grid_dim*grid_dim; ++ii)
		{
			myfile << ii << " " << grid_gradient_y_gpu[ii]<< std::setprecision(15)  << " " << std::endl;
		}
		myfile.close();
#endif
	}
#endif

#endif
}
