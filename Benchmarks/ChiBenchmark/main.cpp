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
#include <mm_malloc.h>
//
#include <structure_hpc.hpp>
#include "timer.h"
#include "gradient.hpp"
#include "chi_CPU.hpp"
#include "module_cosmodistances.hpp"
#include "module_readParameters.hpp"


#ifdef __WITH_MPI
#include<mpi.h>
#endif

#ifdef _OPENMP
#include<omp.h>
#endif


//#define __WITH_LENSTOOL 0
#ifdef __WITH_LENSTOOL
#warning "linking with libtool..."
#include <fonction.h>
#include <constant.h>
#include <dimension.h>
#include <structure.h>
#include <setup.hpp>
#endif

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


void chi_bruteforce_SOA_CPU_grid_gradient(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images);


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
	std::cout << "output directory: " << outdir << std::endl;

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

int main(int argc, char *argv[])
{
	int world_size = 1;
	int world_rank = 0;
#ifdef __WITH_MPI	
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);
	MPI_Barrier(MPI_COMM_WORLD);

	// Print off a hello world message
#endif
	int numthreads = 1;
#ifdef _OPENMP
#warning "using openmp"
#pragma omp parallel
	numthreads = omp_get_num_threads();
#endif
	//
	printf("Hello world from processor %s, rank %d out of %d processors and %d threads per rank\n", processor_name, world_rank, world_size, numthreads); fflush(stdout);
#ifdef __WITH_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	int verbose = (world_rank == 0);
	//
	if (verbose) printf("Lenstool-HPC\n\n");
	//
	double wallclock = myseconds();
	if (world_rank == 0) printf("Reading parameter file at time %f s...\n", myseconds() - wallclock);
	// Setting Up the problem
	//===========================================================================================================

	// This module function reads the terminal input when calling LENSTOOL and checks that it is correct
	// Otherwise it exits LENSTOOL
	if (world_rank == 0) module_readCheckInput_readInput(argc, argv);

	// This module function reads the cosmology parameters from the parameter file
	// Input: struct cosmologicalparameters cosmology, parameter file
	// Output: Initialized cosmology struct
	cosmo_param cosmology;  // Cosmology struct to store the cosmology data from the file
	std::string inputFile = argv[1];   // Input file
	module_readParameters_readCosmology(inputFile, cosmology);

	// This module function reads the runmode paragraph and the number of sources, arclets, etc. in the parameter file.
	// The runmode_param stores the information of what exactly the user wants to do with lenstool.
	struct runmode_param runmode;
	module_readParameters_readRunmode(inputFile, &runmode);

	if (world_rank == 0)
	{
		module_readParameters_debug_cosmology(runmode.debug, cosmology);
		module_readParameters_debug_runmode(runmode.debug, runmode);
	}

	//=== Declaring variables
	struct grid_param frame;
	struct galaxy images[runmode.nimagestot];
	struct galaxy sources[runmode.nsets];
	//struct Potential lenses[runmode.nhalos+runmode.npotfile-1];
	struct Potential_SOA lenses_SOA_table[NTYPES];
	struct Potential_SOA lenses_SOA;
	struct cline_param cline;
	struct potfile_param potfile;
	//struct Potential potfilepotentials[runmode.npotfile];
	struct potentialoptimization host_potentialoptimization[runmode.nhalos];
	int nImagesSet[runmode.nsets]; // Contains the number of images in each set of imagesnano

	// This module function reads in the potential form and its parameters (e.g. NFW)
	// Input: input file
	// Output: Potentials and its parameters

	module_readParameters_PotentialSOA_direct(inputFile, &lenses_SOA, runmode.nhalos, runmode.n_tot_halos, cosmology);
	module_readParameters_debug_potential_SOA(0, lenses_SOA, runmode.nhalos);

	// This module function reads in the potfiles parameters
	// Input: input file
	// Output: Potentials from potfiles and its parameters

	if (runmode.potfile == 1 )
	{
		module_readParameters_readpotfiles_param(inputFile, &potfile, cosmology);
		module_readParameters_debug_potfileparam(1, &potfile);
		module_readParameters_readpotfiles_SOA(&runmode, &cosmology,&potfile,&lenses_SOA);
		module_readParameters_debug_potential_SOA(0, lenses_SOA, runmode.n_tot_halos);

	}

	// This module function reads in the grid form and its parameters
	// Input: input file
	// Output: grid and its parameters

	module_readParameters_Grid(inputFile, &frame);

	if (runmode.image == 1 or runmode.inverse == 1 or runmode.time > 0){

		// This module function reads in the strong lensing images
		module_readParameters_readImages(&runmode, images, nImagesSet);
		//runmode.nsets = runmode.nimagestot;
		for(int i = 0; i < runmode.nimagestot; ++i){

			images[i].dls = module_cosmodistances_objectObject(lenses_SOA.z[0], images[i].redshift, cosmology);
			images[i].dos = module_cosmodistances_observerObject(images[i].redshift, cosmology);
			images[i].dr = module_cosmodistances_lensSourceToObserverSource(lenses_SOA.z[0], images[i].redshift, cosmology);

		}
		//module_readParameters_debug_image(runmode.debug, images, nImagesSet,runmode.nsets);

	}

	/*
	if (runmode.inverse == 1){

		// This module function reads in the potential optimisation limits
		module_readParameters_limit(inputFile,host_potentialoptimization,runmode.nhalos);
		module_readParameters_debug_limit(runmode.debug, host_potentialoptimization[0]);
	}
	*/


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
			sources[i].dls = module_cosmodistances_objectObject(lenses_SOA.z[0], sources[i].redshift, cosmology);
			sources[i].dos = module_cosmodistances_observerObject(sources[i].redshift, cosmology);
			sources[i].dr = module_cosmodistances_lensSourceToObserverSource(lenses_SOA.z[0], sources[i].redshift, cosmology);
		}
		module_readParameters_debug_source(runmode.debug, sources, runmode.nsets);
	}
#ifdef __WITH_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	//
	if (verbose) std::cout << "--------------------------" << std::endl << std::endl;
	//
	// Lenstool Bruteforce
	//===========================================================================================================
	if (verbose)
	{
#ifdef __WITH_LENSTOOL
		{
			printf("Calling lenstool at time %f s\n", myseconds() - wallclock);
			setup_lenstool();
			//
			double chi2;
			double lhood0(0);
			int error(0);
			double time;

			if ( M.ichi2 != 0 )
			{
				int error;
				//NPRINTF(stdout, "INFO: compute chires.dat\n");
				readConstraints();
				o_chires("chires.dat");
				time = -myseconds();
				error = o_chi_lhood0(&chi2, &lhood0, NULL);
				time += myseconds();
				printf("INFO: chi2 %lf  Lhood %lf\n", chi2, -0.5 * ( chi2 + lhood0 )  );
				o_global_free();
			}

			std::cout << "Lenstool 6.8.1 chi Benchmark: ";
			std::cout << " Chi : " << std::setprecision(15) << chi2 ;
			std::cout << " Time  " << std::setprecision(15) << time << std::endl;
			o_global_free();
		}
#endif
	}

	// Lenstool-GPU Bruteforce
	//===========================================================================================================
#if 0
	{
		printf("Calling lenstoolhpc orig at time %f s\n", myseconds() - wallclock);
		std::cout << "LenstoolHPC dist chi Benchmark:\n ";
		double chi2;
		double time;
		int error;
		time = -myseconds();
		chi_bruteforce_SOA_CPU_grid_gradient(&chi2, &error, &runmode, &lenses_SOA, &frame, nImagesSet, images);
		time += myseconds();

		std::cout << " Chi : " << std::setprecision(15) << chi2;
		std::cout << " Time  " << std::setprecision(15) << time << std::endl;
	}
#endif



#if 1
	{
		//std::cout << "MylenstoolHPC chi Benchmark:\n "; 
		if (verbose) printf("Calling lenstoolhpc at time %f s\n", myseconds() - wallclock);
		double chi2;
		double time;
		int error;
		time = -myseconds();
		mychi_bruteforce_SOA_CPU_grid_gradient(&chi2, &error, &runmode, &lenses_SOA, &frame, nImagesSet, images);
		time += myseconds();
		if (verbose)
		{
			std::cout << " Chi : " << std::setprecision(15) << chi2;
			std::cout << " Time  " << std::setprecision(15) << time << std::endl;
		}
	}
#endif

#if 1
	{
		//std::cout << "MylenstoolHPC chi Benchmark:\n ";
		if (verbose) printf("Calling lenstoolhpc with barycenter method at time %f s\n", myseconds() - wallclock);
		double chi2;
		double time;
		int error;
		time = -myseconds();
		mychi_bruteforce_SOA_CPU_grid_gradient_barycentersource(&chi2, &error, &runmode, &lenses_SOA, &frame, nImagesSet, images);
		time += myseconds();
		if (verbose)
		{
			std::cout << " Chi : " << std::setprecision(15) << chi2;
			std::cout << " Time  " << std::setprecision(15) << time << std::endl;
		}
	}
#endif

	if (verbose) printf("Ending execution at time %f s\n", myseconds() - wallclock);
#ifdef __WITH_MPI
	MPI_Finalize();
#endif

}
