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
//
#include <mm_malloc.h>
//
#include "structure.h"
#include "timer.h"
#include "gradient.hpp"
#include "chi.hpp"
#include "module_cosmodistances.h"
#include "module_readParameters.hpp"


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

int main(int argc, char *argv[])
{


	// Setting Up the problem
	//===========================================================================================================

	// This module function reads the terminal input when calling LENSTOOL and checks that it is correct
	// Otherwise it exits LENSTOOL
	module_readCheckInput_readInput(argc, argv);

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

	module_readParameters_debug_cosmology(runmode.debug, cosmology);
	module_readParameters_debug_runmode(runmode.debug, runmode);


	//=== Declaring variables
	struct grid_param frame;
	struct galaxy images[runmode.nimagestot];
	struct galaxy sources[runmode.nsets];
	struct Potential lenses[runmode.nhalos+runmode.npotfile-1];
	struct Potential_SOA lenses_SOA[NTYPES];
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
	module_readParameters_PotentialSOA(inputFile, lenses, lenses_SOA, runmode.Nlens);
	module_readParameters_debug_potential(runmode.debug, lenses, runmode.nhalos);
	std::cerr << lenses_SOA[1].b0[0] << " " << lenses[0].b0  << std::endl;
	// This module function reads in the potfiles parameters
	// Input: input file
	// Output: Potentials from potfiles and its parameters

	if (runmode.potfile == 1 ){
		module_readParameters_readpotfiles_param(inputFile, &potfile);
		module_readParameters_debug_potfileparam(runmode.debug, &potfile);
		module_readParameters_readpotfiles(&runmode,&potfile,lenses);
		module_readParameters_debug_potential(runmode.debug, lenses, runmode.nhalos+runmode.npotfile);

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

			images[i].dls = module_cosmodistances_objectObject(lenses[0].z, images[i].redshift, cosmology);
			images[i].dos = module_cosmodistances_observerObject(images[i].redshift, cosmology);
			images[i].dr = module_cosmodistances_lensSourceToObserverSource(lenses[0].z, images[i].redshift, cosmology);

		}
		module_readParameters_debug_image(runmode.debug, images, nImagesSet,runmode.nsets);

	}

	if (runmode.inverse == 1){

		// This module function reads in the potential optimisation limits
		module_readParameters_limit(inputFile,host_potentialoptimization,runmode.nhalos);
		module_readParameters_debug_limit(runmode.debug, host_potentialoptimization[0]);
	}


	if (runmode.source == 1){
		//Initialisation to default values.(Setting sources to z = 1.5 default value)
		for(int i = 0; i < runmode.nsets; ++i){
			sources[i].redshift = 1.5;
		}
		// This module function reads in the strong lensing sources
		module_readParameters_readSources(&runmode, sources);
		//Calculating cosmoratios
		for(int i = 0; i < runmode.nsets; ++i){

			sources[i].dls = module_cosmodistances_objectObject(lenses[0].z, sources[i].redshift, cosmology);
			sources[i].dos = module_cosmodistances_observerObject(sources[i].redshift, cosmology);
			sources[i].dr = module_cosmodistances_lensSourceToObserverSource(lenses[0].z, sources[i].redshift, cosmology);
		}
		module_readParameters_debug_source(runmode.debug, sources, runmode.nsets);

	}

	std::cout << "--------------------------" << std::endl << std::endl;

	// Lenstool-GPU Bruteforce
	//===========================================================================================================
	double chi2(0);
	int error(0);

	double t_1(0),t_2(0),t_3(0);
#if 0
	t_1 = -myseconds();
	chi_bruteforce(&chi2,&error,&runmode,lenses,&frame,nImagesSet,images);
	t_1 += myseconds();

	std::cout << " Chi Brute Force Benchmark " << std::endl;
	std::cout << " Chi : " << std::setprecision(15) << chi2 <<  std::endl;
	std::cout << " Time  " << std::setprecision(15) << t_1 << std::endl;
#endif

#if 0
	t_2 = -myseconds();
	chi_bruteforce_SOA(&chi2,&error,&runmode,lenses_SOA,&frame,nImagesSet,images);
	t_2 += myseconds();

	std::cout << " Chi Brute Force SOA Benchmark " << std::endl;
	std::cout << " Chi : " << std::setprecision(15) << chi2 <<  std::endl;
	std::cout << " Time  " << std::setprecision(15) << t_2 << std::endl;
	std::cout << " Gain  " << std::setprecision(15) << t_1/t_2 << std::endl;
#endif

#if 1
	t_3 = -myseconds();
	chi_bruteforce_SOA_AVX(&chi2, &error, &runmode,lenses_SOA, &frame, nImagesSet, images);
	t_3 += myseconds();

	std::cout << " Chi Brute Force SOA AVX Benchmark " << std::endl;
	std::cout << " Chi : " << std::setprecision(15) << chi2 <<  std::endl;
	std::cout << " Time  " << std::setprecision(15) << t_3 << std::endl;
	std::cout << " Gain  " << std::setprecision(15) << t_1/t_3 << std::endl;
#endif


}
