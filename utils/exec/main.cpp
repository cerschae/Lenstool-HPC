/**
Lenstool-HPC: HPC based massmodeling software and Lens-map generation
Copyright (C) 2017  Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

@brief: Main Lenstool-HPC executable

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
#include "grid_map_mass_GPU.cuh"
#include "grid_map_dpl_GPU.cuh"
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

	// Checks the terminal input parameters and if the mentioned folder exist.
	// Otherwise it exits LENSTOOL.
	char cwd[1024];
	if (getcwd(cwd, sizeof(cwd)) != NULL)
		fprintf(stdout, "Current working dir: %s\n", cwd);
	//
	std::string path;
	module_readCheckInput_readInput(argc, argv, &path);
	//

	// Reading cosmology parameters from the parameter file
	// 		Input: struct cosmologicalparameters cosmology, parameter file
	// 		Output: Initialized cosmology struct
	cosmo_param cosmology;  				// Cosmology struct to store the cosmology data from the file
	std::string inputFile = argv[1];   		// Input file
	module_readParameters_readCosmology(inputFile, cosmology);
	//
	// Reading the runmode paragraph and the number of sources, arclets, etc. in the parameter file.
	// The runmode_param stores the information of what exactly the user wants to do with lenstool.
	struct runmode_param runmode;
	module_readParameters_readRunmode(inputFile, &runmode);
	module_readParameters_debug_cosmology(runmode.debug, cosmology);
	module_readParameters_debug_runmode(1, runmode);
	//
	// Variable Declaration
	struct grid_param frame;
	struct galaxy images[runmode.nimagestot];
	struct galaxy sources[runmode.nsets];
	struct Potential_SOA lenses_SOA_table[NTYPES];
	struct Potential_SOA lenses_SOA;
	struct cline_param cline;
	struct potfile_param potfile[runmode.Nb_potfile];
	struct potentialoptimization host_potentialoptimization[runmode.nhalos];
	int nImagesSet[runmode.nsets]; 							// Contains the number of images in each set of images
	//Bayesmap specific variables
	type_t* bayespot;
	int nparam, nvalues;

	// Reading major potential Information for the .par file and minor potentials from the potfiles
	module_readParameters_PotentialSOA_direct(inputFile, &lenses_SOA, runmode.nhalos, runmode.n_tot_halos, cosmology);
	module_readParameters_debug_potential_SOA(0, lenses_SOA, runmode.nhalos);
	module_readParameters_limit(inputFile, host_potentialoptimization, runmode.nhalos );

	if (runmode.potfile == 1 )
	{
		module_readParameters_readpotfiles_param(inputFile, potfile, cosmology);
		module_readParameters_debug_potfileparam(0, &potfile[0]);
		module_readParameters_debug_potfileparam(0, &potfile[1]);
		module_readParameters_readpotfiles_SOA(&runmode, &cosmology,potfile,&lenses_SOA);
		module_readParameters_debug_potential_SOA(1, lenses_SOA, runmode.n_tot_halos);

	}
	//Updating cosmological calculation of the potential
	module_readParameters_lens_dslds_calculation(&runmode,&cosmology,&lenses_SOA);

	//Reading in image information and calculating cosmological distance information
	if (runmode.image == 1 or runmode.inverse == 1 or runmode.time > 0){

		// This module function reads in the strong lensing images
		module_readParameters_readImages(&runmode, images, nImagesSet);
		//runmode.nsets = runmode.nimagestot;
		for(int i = 0; i < runmode.nimagestot; ++i){

			images[i].dls = module_cosmodistances_objectObject(lenses_SOA.z[0], images[i].redshift, cosmology);
			images[i].dos = module_cosmodistances_observerObject(images[i].redshift, cosmology);
			images[i].dr = module_cosmodistances_lensSourceToObserverSource(lenses_SOA.z[0], images[i].redshift, cosmology);

		}
		module_readParameters_debug_image(0, images, nImagesSet,runmode.nsets);

	}

	// Reading in the grid type and its parameters
	module_readParameters_Grid(inputFile, &frame);
	//
	std::cout << "--------------------------" << std::endl << std::endl; fflush(stdout);

	//This is the lenstool computation area. Only compiled when linked to the lenstool library
	//It is essentialy used to compute correspond lenstool result for benchmarks and tests
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
#if 1
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

    /* antecedant d'une grille source ou d'une grille image */
    if ( M.grille != 0 ){
    	t_lt = -myseconds();
        g_grid(M.grille, M.ngrille, M.zgrille);
        t_lt += myseconds();}

    /* grille du amplification (keyword ampli)*/
    if ( M.iampli != 0 ){
    	t_lt = -myseconds();
        g_ampli(M.iampli, M.nampli, M.zampli, M.amplifile);
        t_lt += myseconds();}

    /* grille du potential */
    if ( M.ipoten != 0 ){
    	t_lt = -myseconds();
    	std::cerr << " PRINT " << M.ipoten << M.npoten << M.zpoten << M.potenfile << std::endl;
        g_poten(M.ipoten, M.npoten, M.zpoten, M.potenfile);
        t_lt += myseconds();}

    /* grille du mass */
    if ( M.imass != 0 ){
    	t_lt = -myseconds();
        g_mass(M.imass, M.nmass, M.zmass, S.zs, M.massfile);
        t_lt += myseconds();}

    /* grille du dpl */
    if ( M.idpl != 0 ){
    	t_lt = -myseconds();
    	//std::cerr << " PRINT " << M.idpl << M.ndpl << M.zdpl << M.dplxfile <<  M.dplyfile << std::endl;
        g_dpl(M.idpl, M.ndpl, M.zdpl, M.dplxfile, M.dplyfile);
        t_lt += myseconds();}

    /* grille du curv */
    if ( M.icurv != 0 )
        g_curv(M.icurv, M.ncurv, M.zcurv, M.cxxfile, M.cxyfile, M.cyyfile);

    /* grille du shear */
    if ( M.ishear != 0 )
        g_shear(M.ishear, M.nshear, M.zshear, M.shearfile);

    /* grille du time-delay */
    if ( M.itime != 0 )
        g_time(M.itime, M.ntime, M.ztime, M.timefile);

    /* grille du shear_field */
    if ( M.ishearf != 0 )
        g_shearf(M.ishearf, M.zshearf, M.shearffile, M.nshearf);

    /* amplification field grid */
    if ( M.iamplif != 0 )
        g_amplif(M.iamplif, M.zamplif, M.ampliffile);



#endif

	//This is the lenstool-HPC image computation area.
	if (runmode.image == 1){
		galaxy predicted_images[runmode.nimagestot*MAXIMPERSOURCE];
		//predicting_images_bruteforce_barycentersource(&predicted_images, &runmode, &lenses_SOA, &frame, nImagesSet, images);

	}

	//This is the lenstool-HPC map computation area. The following computes the maps for
    //every mode using GPUs.
#ifdef __WITH_GPU
	double t_1, t_2;

	t_1 = -myseconds();
	//std::cerr << runmode.mass << " " << runmode.amplif << " "  << std::endl;
		//Mode: Mass map computation
		if (runmode.mass > 0){
			//Allocation
			type_t* mass_GPU = (type_t *) malloc((int) (runmode.mass_gridcells) * (runmode.mass_gridcells) * sizeof(type_t));

			////calculate maps
			std::cout << " GPU launching for map mass " << std::endl;
			t_2 = -myseconds();
			module_readParameters_debug_potential_SOA(0, lenses_SOA, runmode.n_tot_halos);
			//Init
			memset(mass_GPU, 0, (runmode.mass_gridcells) * (runmode.mass_gridcells) * sizeof(type_t));

			//Choosing Function definition
			map_gpu_function_t map_gpu_func;
			map_gpu_func = select_map_mass_function(&runmode);

			//calculating map using defined function
			map_grid_mass_GPU(map_gpu_func,mass_GPU,&cosmology, &frame, &lenses_SOA, runmode.n_tot_halos, runmode.mass_gridcells ,runmode.mass, runmode.z_mass, runmode.z_mass_s);

			//writing
			module_writeFits(path,runmode.mass_name,mass_GPU,&runmode,runmode.mass_gridcells,&frame, runmode.ref_ra, runmode.ref_dec );
			t_2 += myseconds();
			std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;
			std::cerr << "**" << mass_GPU[0] << std::endl;

			free(mass_GPU);
		}
		//Mode: Amplification map computation
		if (runmode.amplif > 0){
			//Allocation
			type_t* ampli_GPU = (type_t *) malloc((int) (runmode.amplif_gridcells) * (runmode.amplif_gridcells) * sizeof(type_t));

			////calculate maps
			std::cout << " GPU launching for map amplif " << std::endl;
			t_2 = -myseconds();
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
			module_writeFits(path,runmode.amplif_name,ampli_GPU,&runmode,runmode.amplif_gridcells,&frame, runmode.ref_ra, runmode.ref_dec );
			t_2 += myseconds();
			std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;
			std::cerr << "**" << ampli_GPU[0] << std::endl;

			//std::cerr << "**" << ampli_GPU[0] << std::endl;
			free(ampli_GPU);
		}
		//Mode: Shear map computation
		if (runmode.shear > 0){
			//Allocation
			type_t* shear_GPU = (type_t *) malloc((int) (runmode.shear_gridcells) * (runmode.shear_gridcells) * sizeof(type_t));

			////calculate maps
			std::cout << " GPU launching for map shear " << std::endl;
			t_2 = -myseconds();
			module_readParameters_debug_potential_SOA(0, lenses_SOA, runmode.n_tot_halos);
			//Init
			memset(shear_GPU, 0, (runmode.shear_gridcells) * (runmode.shear_gridcells) * sizeof(type_t));

			//Choosing Function definition
			map_gpu_function_t map_gpu_func;
			map_gpu_func = select_map_shear_function(&runmode);

			//calculating map using defined function
			map_grid_shear_GPU(map_gpu_func,shear_GPU,&cosmology, &frame, &lenses_SOA, runmode.n_tot_halos, runmode.shear_gridcells ,runmode.shear, runmode.z_shear);

			//writing
			//std::cerr << runmode.amplif_name << std::endl;
			module_writeFits(path,runmode.shear_name,shear_GPU,&runmode,runmode.shear_gridcells,&frame, runmode.ref_ra, runmode.ref_dec );
			t_2 += myseconds();
			std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;
			std::cerr << "**" << shear_GPU[0] << std::endl;

			//std::cerr << "**" << ampli_GPU[0] << std::endl;
			free(shear_GPU);
		}
		//Mode: Displacement map computation
		if (runmode.dpl > 0){
			//Allocation
			type_t* dpl_x = (type_t *) malloc((int) (runmode.dpl_gridcells) * (runmode.dpl_gridcells) * sizeof(type_t));
			type_t* dpl_y = (type_t *) malloc((int) (runmode.dpl_gridcells) * (runmode.dpl_gridcells) * sizeof(type_t));

			////calculate maps
			std::cout << " GPU launching for map dpl " << runmode.dpl_gridcells << runmode.dpl << runmode.z_dpl <<std::endl;
			t_2 = -myseconds();
			module_readParameters_debug_potential_SOA(0, lenses_SOA, runmode.n_tot_halos);
			//Init
			memset(dpl_x, 0, (runmode.dpl_gridcells) * (runmode.dpl_gridcells) * sizeof(type_t));
			memset(dpl_y, 0, (runmode.dpl_gridcells) * (runmode.dpl_gridcells) * sizeof(type_t));

			//Choosing Function definition
			map_gpu_function_t map_gpu_func;
			//map_gpu_func = select_map_shear_function(&runmode);

			//calculating map using defined function
			map_grid_dpl_GPU(map_gpu_func, dpl_x, dpl_y, &cosmology, &frame, &lenses_SOA, runmode.n_tot_halos, runmode.dpl_gridcells ,runmode.dpl, runmode.z_dpl);

			std::string file_x, file_y;
			file_x = runmode.dpl_name1;
			//file_x.append("_x");
			file_y = runmode.dpl_name2;
			//file_y.append("_y");

			//writing
			module_writeFits(path,file_x,dpl_x,&runmode,runmode.dpl_gridcells,&frame, runmode.ref_ra, runmode.ref_dec );
			module_writeFits(path,file_y,dpl_y,&runmode,runmode.dpl_gridcells,&frame, runmode.ref_ra, runmode.ref_dec );
			t_2 += myseconds();
			std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;
			std::cerr << "**" << dpl_x[0] << std::endl;

			//std::cerr << "**" << ampli_GPU[0] << std::endl;
			free(dpl_x);
			free(dpl_y);
		}
		if (runmode.potential > 0){
			//Allocation
			type_t* pot_GPU = (type_t *) malloc((int) (runmode.pot_gridcells) * (runmode.pot_gridcells) * sizeof(type_t));

			////calculate maps
			std::cout << " GPU launching for map potential ATT: If pot ellip = 0, NAN are imminent " << std::endl;
			t_2 = -myseconds();
			module_readParameters_debug_potential_SOA(0, lenses_SOA, runmode.n_tot_halos);
			//Init
			memset(pot_GPU, 0, (runmode.pot_gridcells) * (runmode.pot_gridcells) * sizeof(type_t));

			//Choosing Function definition
			map_pot_function_t map_pot_function;
			map_pot_function = select_map_potential_function(&runmode);

			//calculating map using defined function
			map_grid_potential_GPU(map_pot_function, pot_GPU, &cosmology, &frame, &lenses_SOA, runmode.n_tot_halos, runmode.pot_gridcells ,runmode.potential, runmode.z_pot);

			std::string file;
			file = runmode.pot_name;
			//writing
			module_writeFits(path,file,pot_GPU,&runmode,runmode.pot_gridcells,&frame, runmode.ref_ra, runmode.ref_dec );
			t_2 += myseconds();
			std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;
			std::cerr << "**" << pot_GPU[0] << std::endl;

			//std::cerr << "**" << ampli_GPU[0] << std::endl;
			free(pot_GPU);
		}
		//

	t_1 += myseconds();
	std::cout << "Lenstool 1 Map Time " << std::setprecision(15) << t_lt << std::endl;
	//std::cout << "HPC Total Time  " << std::setprecision(15) << t_1 << std::endl;
#endif

}

