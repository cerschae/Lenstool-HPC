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
#include "module_writeFits.hpp"
#ifdef __WITH_GPU
#include "grid_gradient_GPU.cuh"
#include "grid_map_GPU.cuh"
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
#include<lt.h>
#include <stdlib.h>

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

double  **alloc_square_double_test(int nbr_lin,int nbr_col)
{
        auto     double   **square;
        register int      i, j;

        square = (double **) malloc((unsigned) nbr_lin*sizeof(double *));
        if (square != 0)
        {
                for (i=0; i<nbr_lin; i++)
            {
                square[i] = (double *)malloc((unsigned) nbr_col*sizeof(double));
                if (square[i] == 0) square = 0;
            }
        }

        for (i=0; i<nbr_lin; i++)
                for (j=0; j<nbr_col; j++)
                        square[i][j]=0.0;

        return(square);
}

int     **alloc_square_int_test(int nbr_lin,int nbr_col)
{
auto     int   **square;
register int      i, j;

square = (int **) malloc((unsigned) nbr_lin*sizeof(int *));
if (square != 0)
  {
  for (i=0; i<nbr_lin; i++)
    {
    square[i] = (int *)malloc((unsigned) nbr_col*sizeof(int));
    if (square[i] == 0) square = 0;
    }
  }

for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
square[i][j]=0;

return(square);
}

void    free_square_double_test(double **square,int nbr_lin)
{
        register  int    i;
        if(square!=NULL)
        {
            for (i=0; i<nbr_lin; i++)
                free(square[i]);

                free((double *) square);
        }
}

void    free_square_int_test(int **square,int nbr_lin)
{
        register  int    i;

        for (i=0; i<nbr_lin; i++)
                free((int *) square[i]);

        free((int *) square);
}



void
gradient_grid_GPU_sorted(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, int nbgridcells);
//
//
int module_readCheckInput_readInput(int argc, char *argv[], std::string *outdir)
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
	//std::string outdir = argv[2];
	*outdir = argv[2];
	*outdir += "-";
	*outdir += trimstamp;
	std::cout << *outdir << std::endl;

	// check whether the output directory already exists
	if (stat(outdir->c_str(), &file_stat) < 0){
		mkdir(outdir->c_str(), S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH );
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
	std::string path;
	module_readCheckInput_readInput(argc, argv, &path);
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
		module_readParameters_readpotfiles_SOA(&runmode, &cosmology,&potfile,&lenses_SOA);
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
	//
	type_t iamp = 5;

#ifdef __WITH_LENSTOOL
        std::cout << " CPU Test Lenstool    ... ";
        //type_t *ampli;
        //ampli = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));

        F.xmin =  F.ymin = frame.xmin;
        F.xmax = F.ymax = frame.xmax;
        G.nlens = runmode.nhalos;
        double **ampli;
        int **namp;
        double z = runmode.z_amplif;
        int np = runmode.nbgridcells;
        double dl0s = module_cosmodistances_objectObject(lens[0].z, z, cosmology);
        double dos = module_cosmodistances_observerObject(z, cosmology);
        double dlsds = dl0s / dos;
        point pi;
        int i,j;

        double t_lt = -myseconds();
//#pragma omp parallel for if (omp_get_num_threads() > 1) schedule(guided, 100)
//#pragma omp parallel for
#if 1
        ampli = (double **) alloc_square_double_test(np, np);
        namp = (int **) alloc_square_int_test(np, np);

        /* Make sure we have empty arrays */
        for (j = 0; j < np; j++)
            for (i = 0; i < np; i++)
            {
                ampli[i][j] = 0.;
                namp[i][j] = 0;
            }

        if (iamp > 0)
        {
//#pragma omp parallel for
            for (j = 0; j < np; j++)
            {
                struct matrix MA;
                struct ellipse amp;
                double kappa, ga1,ga2,gam,gp;
                pi.y = j * (F.ymax - F.ymin) / (np - 1) + F.ymin;
                for (i = 0; i < np; i++)
                {
                    pi.x = i * (F.xmax - F.xmin) / (np - 1) + F.xmin;
                    amp = e_unmag(&pi, dl0s, dos, z);
                    /*amplification*/
                    if (iamp == 1)
                        ampli[j][i] = 1. / (amp.a * amp.b);
                    /*absolute value of amplification*/
                    else if (iamp == 2)
                        ampli[j][i] = 1. / fabs(amp.a * amp.b);
                    /*amplification in magnitudes*/
                    else if (iamp == 3)
                        ampli[j][i] = -2.5 * log10(fabs(amp.a * amp.b));
                    /**/
                    else if (iamp == 4)
                    {
                        MA = e_grad2(&pi, dl0s, z);
                        MA.a /= dos;
                        MA.b /= dos;
                        MA.c /= dos;

                        kappa = (MA.a + MA.c) / 2.;
                        ga1 = (MA.a - MA.c) / 2.;
                        ga2 = MA.b;
                        gam = sqrt(ga1 * ga1 + ga2 * ga2); /*gamma*/
                        gp = gam / (1 - kappa);
                        ampli[j][i] = (1 - kappa) * (1 + gp * gp) / (1 - gp * gp);
                    }
                    else if (iamp == 5 || iamp == 6)
                    {
                        MA = e_grad2(&pi, dl0s, z);
                        MA.a /= dl0s;
                        MA.b /= dl0s;
                        MA.c /= dl0s;

                        kappa = (MA.a + MA.c) / 2.;
                        ga1 = (MA.a - MA.c) / 2.;
                        ga2 = MA.b;
                        gam = sqrt(ga1 * ga1 + ga2 * ga2);
                        if (iamp == 5)
                            ampli[j][i] = kappa;
                        else if (iamp == 6)
                            ampli[j][i] = gam;
                    }
                    /*amplification^-1*/
                    else
                        ampli[j][i] = (amp.a * amp.b);
                };
            };
        }
#endif
	t_lt += myseconds();
	std::cout << " Time = " << t_lt << " s." << std::endl;
#endif

	//
	std::cout << " CPU Test lenstool_hpc... "; 
	//
	 type_t *ampli_CPU;
	ampli_CPU = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	memset(ampli_CPU, 0, (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	int Nstat = 1;
	t_1 = -myseconds();
	for(int ii = 0; ii < Nstat; ++ii) {
		amplif_grid_CPU(ampli_CPU, &cosmology, &frame, &lenses_SOA, runmode.nhalos, grid_dim, runmode.amplif, runmode.z_amplif);
		}
	t_1 += myseconds();
	//
	std::cout << " Time = " << std::setprecision(15) << t_1 << std::endl;
#if 1
#ifdef __WITH_GPU
	// GPU test
	std::cout << " GPU Test... "; 
	//
	type_t* ampli_GPU = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	//
	memset(ampli_GPU, 0, (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	//
	ampli_GPU = (type_t *) malloc((int) (runmode.nbgridcells) * (runmode.nbgridcells) * sizeof(type_t));
	//
	t_2 = -myseconds();
	for(int ii = 0; ii < Nstat; ++ii) {
		map_gpu_function_t map_gpu_func = &amplif_grid_CPU_GPU;
		map_grid_GPU(map_gpu_func,ampli_GPU,&cosmology, &frame, &lenses_SOA, runmode.nhalos, grid_dim,runmode.amplif, runmode.z_amplif);
	}
	std::string file;
	file = path;
	file.append("/amplif");
	file.append(".fits");

	char file_char[file.length()+1];
	strcpy(file_char,file.c_str());

	module_writeFits_Image(file_char,ampli_GPU,grid_dim,grid_dim,frame.xmin,frame.xmax,frame.ymin,frame.ymax);

	//free(amplif);

	t_2 += myseconds();
	std::cerr << "**" << ampli_GPU[0] << std::endl;
	std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;

#endif
#endif
	std::ofstream myfile;
#ifdef __WITH_LENSTOOL 
	{
		type_t norm_a = 0.;
		//
		for (int ii = 0; ii < grid_dim; ++ii)
		{
			for (int jj = 0; jj < grid_dim; ++jj)
			{
				//std::cerr<< ii << " "<<  jj << " "  << ii*grid_dim +jj << std::endl;
				//std::cerr << ampli[ii][jj] << " ";
				norm_a += (ampli[ii][jj] - ampli_CPU[ii*grid_dim+jj])*(ampli[ii][jj] - ampli_CPU[ii*grid_dim+jj]);
			}
			//std::cerr << std::endl;
		}
		//
		std::cout << "  l2 difference norm cpu = " << std::setprecision(15) << norm_a << std::endl;
	}
#endif


	//
#if 1
#ifdef __WITH_GPU
	{
		type_t norm_a = 0.;
		//
		for (int ii = 0; ii < grid_dim; ++ii)
		{
			for (int jj = 0; jj < grid_dim; ++jj)
			{
				//std::cerr<< ii << " "<<  jj << " "  << ii*grid_dim +jj << std::endl;
				//std::cerr << ampli_GPU[ii*grid_dim+jj] << " " << ampli[ii][jj] << std::endl;
			norm_a += (ampli[ii][jj]  - ampli_GPU[ii*grid_dim+jj])*(ampli[ii][jj] - ampli_GPU[ii*grid_dim+jj]);
			}
			//std::cerr << std::endl;
		}
		//
		std::cout << "  l2 difference norm gpu = " << std::setprecision(15) << norm_a << std::endl;
	}
    free_square_double_test(ampli, np);
    free_square_int_test(namp, np);
#endif
#endif

}
