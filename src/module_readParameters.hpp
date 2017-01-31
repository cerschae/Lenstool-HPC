/**
* @file   module_readParameters.h
* @Author Thomas Jalabert, EPFL (me@example.com), Christoph Schaaefer, EPFL (christophernstrerne.schaefer@epfl.ch)
* @date   July 2015
* @version 0,1
* @brief  Header file for readParameter module
*
* read parameter file
* 
*/

// header guard
#ifndef MODULE_READPARAMETERS_H
#define MODULE_READPARAMETERS_H




// Include
#include <iostream>
//#include "module_readParameters.hpp"
#include <cstring>
#include <sstream>
#include <stdlib.h>
#include <structure_hpc.h>




// Function declarations
void module_readParameters_readCosmology(std::string infile, cosmo_param &Cosmology);
void module_readParameters_readRunmode(std::string infile, struct runmode_param *Runmode_param);
void module_readParameters_readImages(const struct runmode_param *runmode, galaxy image[], int nImagesSet[]);
void module_readParameters_readSources(struct runmode_param *runmode, struct galaxy source[]);
void module_readParameters_Grid(std::string infile, grid_param *grid);
void module_readParameters_CriticCaustic(std::string infile, cline_param *cline);
void module_readParameters_arclets(std::string arclets_filename, point arclets_position[], ellipse arclets_shape[], double arclets_redshift[]);
void module_readParameters_limit(std::string infile, struct potentialoptimization host_potentialoptimization[], int nhalos );
void module_readParameters_Potential(std::string infile, Potential lens[], int nhalos);
void module_readParameters_PotentialSOA(std::string infile, Potential *lens, Potential_SOA *lens_SOA, int Nset[]);
void module_readParameters_PotentialSOA_nonsorted(std::string infile, Potential *lens, Potential_SOA *lens_SOA, int nhalos);
void module_readParameters_calculatePotentialparameter(Potential *ilens);
void module_readParameters_SingleLensingSourcesNumberSets(std::string infile, int &nsetofimages_cleanlens );
void module_readParameters_SingleLensingSources(std::string infile, point sources[], ellipse sources_shape[], double redshift[], int nimages_cleanlens[], int nsetofimages_cleanlens );
void module_readParameters_readpotfiles_param(std::string infile, potfile_param *potfile);
void module_readParameters_readpotfiles(const runmode_param *runmode, potfile_param *potfile, Potential *lens);

void module_readParameters_debug_cosmology(int DEBUG, cosmo_param cosmology);
void module_readParameters_debug_runmode(int DEBUG, runmode_param runmode);
void module_readParameters_debug_image(int DEBUG, galaxy image[],int nImagesSet[],int nsets);
void module_readParameters_debug_source(int DEBUG, galaxy source[], int nsets);
void module_readParameters_debug_potential(int DEBUG, Potential potential[], int nhalos);
void module_readParameters_debug_potfileparam(int DEBUG, potfile_param *potfile);
void module_readParameters_debug_criticcaustic(int DEBUG, cline_param cline);
void module_readParameters_debug_limit(int DEBUG, struct potentialoptimization host_potentialoptimization);




#endif // enf of header guard


