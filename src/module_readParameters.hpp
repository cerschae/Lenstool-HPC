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
#include "module_cosmodistances.hpp"
#include <cstring>
#include <sstream>
#include <stdlib.h>
#include <structure_hpc.hpp>




////Function (Master) declarations: Call Slave functions for reading
void module_readParameters_readCosmology(std::string infile, cosmo_param &Cosmology);
void module_readParameters_readRunmode(std::string infile, struct runmode_param *Runmode_param);
void module_readParameters_readImages(const struct runmode_param *runmode, galaxy image[], int nImagesSet[]);
void module_readParameters_readSources(struct runmode_param *runmode, struct galaxy source[]);
void module_readParameters_Grid(std::string infile, grid_param *grid);
void module_readParameters_CriticCaustic(std::string infile, cline_param *cline);
void module_readParameters_arclets(std::string arclets_filename, point arclets_position[], ellipse arclets_shape[], double arclets_redshift[]);
void module_readParameters_limit(std::string infile, struct potentialoptimization host_potentialoptimization[], int nhalos );
void module_readParameters_Potential(std::string infile, Potential lens[], int nhalos);
void module_readParameters_PotentialSOA_2(std::string infile, Potential *lens, Potential_SOA *lens_SOA, int Nset[]);
void module_readParameters_PotentialSOA(std::string infile, Potential *lens, Potential_SOA *lens_SOA, int nhalos);
void module_readParameters_PotentialSOA_direct(std::string infile, Potential_SOA *lens_SOA, int nhalos, int n_tot_halos, cosmo_param cosmology);
//void module_readParameters_PotentialSOA_nonsorted(std::string infile, Potential *lens, Potential_SOA *lens_SOA, int nhalos);
void module_readParameters_calculatePotentialparameter(Potential *ilens);
void module_readParameters_calculatePotentialparameter_SOA(Potential_SOA *lens, int ind);
void module_readParameters_SingleLensingSourcesNumberSets(std::string infile, int &nsetofimages_cleanlens );
void module_readParameters_SingleLensingSources(std::string infile, point sources[], ellipse sources_shape[], double redshift[], int nimages_cleanlens[], int nsetofimages_cleanlens );
void module_readParameters_readpotfiles_param(std::string infile, potfile_param potfile[], cosmo_param cosmology);
void module_readParameters_readpotfiles(const runmode_param *runmode, potfile_param potfile[], Potential *lens);
void module_readParameters_readpotfiles_SOA(const runmode_param *runmode, const cosmo_param *cosmology, potfile_param potfile[], Potential_SOA *lens);
//bayesmap specific functions
void module_readParameters_preparebayes(int &nparam, int &nvalues);
void module_readParameters_bayesmodels(double * bayespot, int nparam, int nvalues);
void setScalingRelations(const runmode_param *runmode, const cosmo_param *cosmology, potfile_param *pot, Potential_SOA* lenses, int index);
//void module_readParameters_setbayesmapmodels( Potential_SOA* lenses, const potentialoptimization* limit, double * bayespot, int nparam, int index, int nhalos);
//void module_readParameters_setbayesmapmodels( Potential_SOA* lenses, const runmode_param* runmode, const potentialoptimization* limit, const potfile_param* potfile, double * bayespot, int nparam, int index);
void module_readParameters_setbayesmapmodels(const runmode_param* runmode, const cosmo_param* cosmology, const potentialoptimization* limit, potfile_param* potfile, Potential_SOA* lenses, double * bayespot, int nparam, int index);
void module_readParameters_lens_dslds_calculation(const runmode_param* runmode, const cosmo_param* cosmo, Potential_SOA* lens);

////Function (Slave) declarations
//for read_Runmode
void read_runmode(std::istream &IN, struct runmode_param *runmode);
void read_runmode_potential(std::istream &IN, int &numberPotentials);
void read_runmode_image(std::istream &IN, struct runmode_param *runmode);
void read_runmode_potfile(std::istream &IN, struct runmode_param *runmode);
void read_runmode_countimages(struct runmode_param *runmode);
void read_runmode_countsources(struct runmode_param *runmode);
void read_runmode_countpotfile(struct runmode_param *runmode);


//Debug and printing functions
void module_readParameters_debug_cosmology(int DEBUG, cosmo_param cosmology);
void module_readParameters_debug_runmode(int DEBUG, runmode_param runmode);
void module_readParameters_debug_image(int DEBUG, galaxy image[],int nImagesSet[],int nsets);
void module_readParameters_debug_source(int DEBUG, galaxy source[], int nsets);
void module_readParameters_debug_potential(int DEBUG, Potential potential[], int nhalos);
void module_readParameters_debug_potential_SOA(int DEBUG, Potential_SOA lenses, int nhalos);
void module_readParameters_debug_potfileparam(int DEBUG, potfile_param *potfile);
void module_readParameters_debug_criticcaustic(int DEBUG, cline_param cline);
void module_readParameters_debug_limit(int DEBUG, struct potentialoptimization host_potentialoptimization);




#endif // enf of header guard


