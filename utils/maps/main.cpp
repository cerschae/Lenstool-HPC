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
#include "gradient2.hpp"
#include "chi_CPU.hpp"
#include "module_cosmodistances.hpp"
#include "module_readParameters.hpp"
#include "grid_gradient2_CPU.hpp"
#include "grid_amplif_CPU.hpp"
#include "module_writeFits.hpp"
#ifdef __WITH_GPU
#include "grid_gradient_GPU.cuh"
#include "grid_map_ampli_GPU.cuh"
#include "grid_map_pot_GPU.cuh"
#include "grid_map_shear_GPU.cuh"
#include "grid_map_dpl_GPU.cuh"
#include "grid_map_mass_GPU.cuh"
#include "grid_gradient2_GPU.cuh"
//#include "gradient_GPU.cuh"
#endif

#ifdef __WITH_LENSTOOL
//#include "setup.hpp"
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
	//*outdir += "-";
	//*outdir += trimstamp;
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
	module_readParameters_debug_runmode(1, runmode);
	//
	//=== Declaring variables
	//
	struct grid_param frame;
	struct galaxy images[runmode.nimagestot];
	struct galaxy sources[runmode.nsets];
	//struct Potential lenses[runmode.n_tot_halos];
	struct Potential_SOA lenses_SOA_table[NTYPES];
	struct Potential_SOA lenses_SOA;
	struct cline_param cline;
	struct potfile_param potfile[runmode.Nb_potfile];
	//struct Potential potfilepotentials[runmode.npotfile];
	struct potentialoptimization host_potentialoptimization[runmode.nhalos];
	int nImagesSet[runmode.nsets]; // Contains the number of images in each set of images
	//Bayesmap specific variables
	type_t* bayespot;
	int nparam, nvalues;

	// This module function reads in the potential form and its parameters (e.g. NFW)
	// Input: input file
	// Output: Potentials and its parameters

	module_readParameters_PotentialSOA_direct(inputFile, &lenses_SOA, runmode.nhalos, runmode.n_tot_halos, cosmology);
	module_readParameters_debug_potential_SOA(1, lenses_SOA, runmode.nhalos);
	module_readParameters_limit(inputFile, host_potentialoptimization, runmode.nhalos );
#if 0
	for(int ii = 0; ii <runmode.nhalos ; ii++){
		module_readParameters_debug_limit(1, host_potentialoptimization[ii]);
	}
#endif
	if (runmode.potfile == 1 )
	{
		module_readParameters_readpotfiles_param(inputFile, potfile, cosmology);
		module_readParameters_debug_potfileparam(1, &potfile[0]);
		module_readParameters_debug_potfileparam(1, &potfile[1]);
		module_readParameters_readpotfiles_SOA(&runmode, &cosmology,potfile,&lenses_SOA);
		module_readParameters_debug_potential_SOA(0, lenses_SOA, runmode.n_tot_halos);

	}

	module_readParameters_lens_dslds_calculation(&runmode,&cosmology,&lenses_SOA);

	// This module function reads in the grid form and its parameters
	// Input: input file
	// Output: grid and its parameters
	//
	module_readParameters_Grid(inputFile, &frame);
	//std::cerr <<frame.xmin <<std::endl;
	//
	std::cout << "--------------------------" << std::endl << std::endl; fflush(stdout);
	//
	double t_lt, t_lt_total;
	int turn = 0;
#ifdef __WITH_LENSTOOL
	double **array; // contains the bayes.dat data
	int nParam;
    long int iVal, nVal;  // size of array
	char fname[50]; // <map><ival>.fits
	char fname2[50]; // <map><ival>.fits
	FILE *pFile;
	int i;
	double *index;  // list of bayes.dat lines
	int    seed;   // random seed
	int tmp;

	printf("Setting up lenstool using %d lenses...", runmode.n_tot_halos); fflush(stdout);
	//convert_to_LT(&lenses_SOA, runmode.nhalos+runmode.npotfile);
	// Read the .par file

	init_grille(argv[1], 1);
	// remove the .fits extension tcpo filename
	if( M.imass ) M.massfile[strlen(M.massfile)-5]=0;
	if( M.ishear ) M.shearfile[strlen(M.shearfile)-5]=0;
	if( M.iampli ) M.amplifile[strlen(M.amplifile)-5]=0;
	if( M.idpl )
	{
	       M.dplxfile[strlen(M.dplxfile)-5]=0;
	       M.dplyfile[strlen(M.dplyfile)-5]=0;
	}
    if( M.pixel ) M.pixelfile[strlen(M.pixelfile)-5]=0;
    if( M.iclean ) ps.pixfile[strlen(ps.pixfile)-5]=0;

	// Read catalog of multiple images
	readConstraints();

	// Initialise the grid
	if( G.pol != 0 )
		gridp();
	else
		grid();

	// Switch to silent mode
	M.verbose = 0;
	printf("ok\n");
	std::cerr << " Read Bayes models"  << std::endl;

	// Read the bayes.dat file
	array = readBayesModels(&nParam, &nVal);
	if( array == NULL )
	{
		fprintf(stderr, "ERROR: bayes.dat file not found\n");
		return -1;
	}

#if 0
	for(int i = 0; i < G.nlens; i++){
		printf("Lenstool Potential[%d]: x = %f, y = %f, vdisp = %f, type = %d \n \t ellipticity = %f, ellipticity_pot = %f, ellipticity angle (radians) = %f, rcore = %f, rcut = %f,\n z = %f\n", i,lens[i].C.x, lens[i].C.y, lens[i].sigma, lens[i].type, lens[i].emass, lens[i].epot, lens[i].theta, lens[i].rc, lens[i].rcut, lens[i].z);
	}
#endif
	// Create the ./tmp directory
	i = system("mkdir -p tmp");

	// Prepare the index list
	index = (double *) malloc((unsigned) nVal*sizeof(double));
	for( i = 0 ; i < nVal ; i++ ) index[i]=i;
	seed = -2;

	std::cerr << " Finished setting up"  << std::endl;

	//Defining maps
	int ampli = 1;
	t_lt_total = -myseconds();
	// Loop over each line
	for( i = 0; i < nVal && i < 2000; i++ )
	{
		// Randomly draw a line from index array
		tmp = (int) floor(d_random(&seed) * (nVal - i));
		iVal = index[i+tmp];
		// and swap the indexes in the index list
		index[i+tmp] = index[i];

		// Set the lens parameters from <array>
		setBayesModel( iVal, nVal, array );
#if 0
		for(int i = 0; i < G.nlens; i++){
			printf("Lenstool Potential[%d]: x = %f, y = %f, vdisp = %f, type = %d \n \t ellipticity = %f, ellipticity_pot = %f, ellipticity angle (radians) = %f, rcore = %f, rcut = %f,\n z = %f\n", i,lens[i].C.x, lens[i].C.y, lens[i].sigma, lens[i].type, lens[i].emass, lens[i].epot, lens[i].theta, lens[i].rc, lens[i].rcut, lens[i].z);
		}
#endif
		if( M.imass != 0 )
		{
			sprintf( fname, "tmp/%s%ld.fits",M.massfile, iVal );
			printf("INFO: Compute file %d/%ld : %s [CTRL-C to interrupt]\n",i+1, nVal,fname);
			fflush(stdout);
			pFile = fopen( fname, "r" );
			if( pFile == NULL )
			{
				pFile = fopen( fname, "w");
				fprintf( pFile, "busy\n" );
				fclose(pFile);
				t_lt = -myseconds();
				g_mass( M.imass, M.nmass, M.zmass, S.zs, fname );
				t_lt += myseconds();
				std::cout << " Time  = " << std::setprecision(15) << t_lt << " " << turn <<std::endl;
				//std::cerr <<" para : " << M.zmass << " " << S.zs << " " <<  distcosmo2(M.zmass, S.zs) << distcosmo1(S.zs) << std::endl;
			}
			else
				fclose(pFile);
		}

		if( M.iampli != 0 )
		{
			sprintf( fname, "tmp/%s%ld.fits","Amplif_", iVal );
			printf("INFO: Compute file %d/%ld : %s [CTRL-C to interrupt]\n",i+1, nVal,fname);
			fflush(stdout);
			pFile = fopen( fname, "r" );
			if( pFile == NULL )
			{
				pFile = fopen( fname, "w");
				fprintf( pFile, "busy\n" );
				fclose(pFile);
				//std::cerr <<  runmode.amplif<< runmode.amplif_gridcells<< runmode.z_amplif << std::endl;
				t_lt = -myseconds();
				g_ampli( M.iampli, M.nampli, M.zampli, fname );
				//g_ampli( runmode.amplif, runmode.amplif_gridcells, runmode.z_amplif, fname );
				t_lt += myseconds();
				//
				turn += 1;
				std::cout << " Time  = " << std::setprecision(15) << t_lt << " " << turn <<std::endl;
			}
			else
				fclose(pFile);
		}

	}
	t_lt_total += myseconds();

#endif

#ifdef __WITH_GPU
	double t_1, t_2;
	//struct matrix *grid_gradient2_cpu;
	//grid_gradient2_cpu = (struct matrix *) malloc((int) (runmode.amplif_gridcells) * (runmode.amplif_gridcells) * sizeof(struct matrix));
	// Bayes Map specific functions
	////read bayes lines
	module_readParameters_preparebayes(nparam, nvalues);
	//std::cerr << nparam << "BLA" << nvalues << std::endl;
	bayespot = (type_t *) malloc((int) (nparam) * (nvalues) * sizeof(type_t));
	module_readParameters_bayesmodels(bayespot, nparam, nvalues);
	////read bayes lines
	//std::cerr <<  "BLA" << std::endl;
	t_1 = -myseconds();

		if (runmode.mass > 0){
			//Allocation
			type_t* mass_GPU = (type_t *) malloc((int) (runmode.mass_gridcells) * (runmode.mass_gridcells) * sizeof(type_t));

			for(int ii = 0; ii < nvalues; ii++){
				////calculate maps
				std::cout << " GPU launching for map mass " << ii << std::endl;
				t_2 = -myseconds();
				////set bayes potential
				module_readParameters_setbayesmapmodels(&runmode, &cosmology, host_potentialoptimization, potfile, &lenses_SOA,bayespot,nparam, ii);
				//std::cerr << runmode.n_tot_halos << std::endl;
				module_readParameters_debug_potential_SOA(0, lenses_SOA, runmode.n_tot_halos);
				//Init
				memset(mass_GPU, 0, (runmode.mass_gridcells) * (runmode.mass_gridcells) * sizeof(type_t));

				//Choosing Function definition
				map_gpu_function_t map_gpu_func;
				map_gpu_func = select_map_mass_function(&runmode);

				//calculating map using defined function
				map_grid_mass_GPU(map_gpu_func,mass_GPU,&cosmology, &frame, &lenses_SOA, runmode.n_tot_halos, runmode.mass_gridcells ,runmode.mass, runmode.z_mass, runmode.z_mass_s);

				std::cerr <<" para : " << runmode.z_mass << " " << runmode.z_mass_s << std::endl;
				//writing
				//std::cerr << runmode.amplif_name << std::endl;
				module_writeFits(path,runmode.mass_name,ii,mass_GPU,&runmode,runmode.mass_gridcells,&frame, runmode.ref_ra, runmode.ref_dec );
				t_2 += myseconds();
				std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;
				//std::cerr << "**" << mass_GPU[0] << std::endl;
			}
			//std::cerr << "**" << ampli_GPU[0] << std::endl;
			free(mass_GPU);
		}
		if (runmode.amplif > 0){
			//Allocation
			type_t* ampli_GPU = (type_t *) malloc((int) (runmode.amplif_gridcells) * (runmode.amplif_gridcells) * sizeof(type_t));
			for(int ii = 0; ii < nvalues; ii++){
				////calculate maps
				std::cout << " GPU launching for map amplif " << ii << std::endl;
				t_2 = -myseconds();
				////set bayes potential
				module_readParameters_setbayesmapmodels(&runmode, &cosmology, host_potentialoptimization, potfile, &lenses_SOA,bayespot,nparam, ii);
				module_readParameters_debug_potential_SOA(0, lenses_SOA, runmode.n_tot_halos);
				//Init
				memset(ampli_GPU, 0, (runmode.amplif_gridcells) * (runmode.amplif_gridcells) * sizeof(type_t));

				//Choosing Function definition
				map_gpu_function_t map_gpu_func;
				map_gpu_func = select_map_ampli_function(&runmode);

				//calculating map using defined function
				map_grid_ampli_GPU(map_gpu_func,ampli_GPU,&cosmology, &frame, &lenses_SOA, runmode.n_tot_halos, runmode.amplif_gridcells ,runmode.amplif, runmode.z_amplif);

				//writing
				//std::cerr << runmode.amplif_name << std::endl;
				module_writeFits(path,runmode.amplif_name,ii,ampli_GPU,&runmode,runmode.amplif_gridcells,&frame, runmode.ref_ra, runmode.ref_dec );
				t_2 += myseconds();
				std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;
				//std::cerr << "**" << ampli_GPU[0] << std::endl;
			}
			//std::cerr << "**" << ampli_GPU[0] << std::endl;
			free(ampli_GPU);
		}
		if (runmode.shear > 0){
			//Allocation
			type_t* shear_GPU = (type_t *) malloc((int) (runmode.shear_gridcells) * (runmode.shear_gridcells) * sizeof(type_t));
			for(int ii = 0; ii < nvalues; ii++){
				////calculate maps
				std::cout << " GPU launching for map shear " << ii << std::endl;
				t_2 = -myseconds();
				////set bayes potential
				module_readParameters_setbayesmapmodels(&runmode, &cosmology, host_potentialoptimization, potfile, &lenses_SOA,bayespot,nparam, ii);
				module_readParameters_debug_potential_SOA(0, lenses_SOA, runmode.n_tot_halos);
				//Init
				memset(shear_GPU, 0, (runmode.shear_gridcells) * (runmode.shear_gridcells) * sizeof(type_t));

				//Choosing Function definition
				map_gpu_function_t map_gpu_func;
				map_gpu_func = select_map_shear_function(&runmode);

				//calculating map using defined function
				map_grid_shear_GPU(map_gpu_func,shear_GPU,&cosmology, &frame, &lenses_SOA, runmode.n_tot_halos, runmode.shear_gridcells ,runmode.shear, runmode.z_shear);

				//writing
				module_writeFits(path,runmode.shear_name,ii,shear_GPU,&runmode,runmode.shear_gridcells,&frame, runmode.ref_ra, runmode.ref_dec );
				t_2 += myseconds();
				std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;
				//std::cerr << "**" << shear_GPU[0] << std::endl;
			}
			free(shear_GPU);
		}
		if (runmode.dpl > 0){
			//Allocation
			type_t* dpl_x = (type_t *) malloc((int) (runmode.dpl_gridcells) * (runmode.dpl_gridcells) * sizeof(type_t));
			type_t* dpl_y = (type_t *) malloc((int) (runmode.dpl_gridcells) * (runmode.dpl_gridcells) * sizeof(type_t));
			for(int ii = 0; ii < nvalues; ii++){

				////calculate maps
				std::cout << " GPU launching for map dpl " << ii << std::endl;
				t_2 = -myseconds();
				////set bayes potential
				module_readParameters_setbayesmapmodels(&runmode, &cosmology, host_potentialoptimization, potfile, &lenses_SOA,bayespot,nparam, ii);
				module_readParameters_debug_potential_SOA(0, lenses_SOA, runmode.n_tot_halos);
				//Init
				memset(dpl_x, 0, (runmode.dpl_gridcells) * (runmode.dpl_gridcells) * sizeof(type_t));
				memset(dpl_y, 0, (runmode.dpl_gridcells) * (runmode.dpl_gridcells) * sizeof(type_t));

				//Choosing Function definition
				map_gpu_function_t map_gpu_func;
				//map_gpu_func = select_map_dpl_function(&runmode);

				//calculating map using defined function
				map_grid_dpl_GPU(map_gpu_func, dpl_x, dpl_y, &cosmology, &frame, &lenses_SOA, runmode.n_tot_halos, runmode.dpl_gridcells ,runmode.dpl, runmode.z_dpl);

				std::string file_x, file_y;
				file_x = runmode.dpl_name1;
				file_x.append("_x");
				file_y = runmode.dpl_name2;
				file_y.append("_y");

				//writing
				module_writeFits(path,file_x,dpl_x,&runmode,runmode.dpl_gridcells,&frame, runmode.ref_ra, runmode.ref_dec );
				module_writeFits(path,file_y,dpl_y,&runmode,runmode.dpl_gridcells,&frame, runmode.ref_ra, runmode.ref_dec );
				t_2 += myseconds();
				std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;
				//std::cerr << "**" << dpl_x[0] << std::endl;
			}

			//std::cerr << "**" << ampli_GPU[0] << std::endl;
			free(dpl_x);
			free(dpl_y);
		}
		if (runmode.potential > 0){
			//Allocation
			type_t* pot_GPU = (type_t *) malloc((int) (runmode.pot_gridcells) * (runmode.pot_gridcells) * sizeof(type_t));
			for(int ii = 0; ii < nvalues; ii++){
				////calculate maps
				std::cout << " GPU launching for map potential " << ii << std::endl;
				t_2 = -myseconds();
				////set bayes potential
				module_readParameters_setbayesmapmodels(&runmode, &cosmology, host_potentialoptimization, potfile, &lenses_SOA,bayespot,nparam, ii);
				module_readParameters_debug_potential_SOA(0, lenses_SOA, runmode.n_tot_halos);
				//Init
				memset(pot_GPU, 0, (runmode.pot_gridcells) * (runmode.pot_gridcells) * sizeof(type_t));

				//Choosing Function definition
				map_pot_function_t map_pot_function;
				map_pot_function = select_map_potential_function(&runmode);

				//calculating map using defined function
				map_grid_potential_GPU(map_pot_function, pot_GPU, &cosmology, &frame, &lenses_SOA, runmode.n_tot_halos, runmode.pot_gridcells ,runmode.potential, runmode.z_pot);

				//writing
				module_writeFits(path,runmode.pot_name,ii,pot_GPU,&runmode,runmode.pot_gridcells,&frame, runmode.ref_ra, runmode.ref_dec );
				t_2 += myseconds();
				std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;
				//std::cerr << "**" << pot_GPU[0] << std::endl;
			}
			free(pot_GPU);
		}
		//

	t_1 += myseconds();
	std::cout << "Lenstool Total Time " << std::setprecision(15) << t_lt_total << std::endl;
	std::cout << "HPC Total Time  " << std::setprecision(15) << t_1 << std::endl;
#endif

}

