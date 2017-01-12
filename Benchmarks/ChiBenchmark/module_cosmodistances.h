
/**
* @file   module_cosmodistances.h
* @date   July 2015
* @version 0,1
* @brief  Header file for cosmoditace calculation module
*
* header file
* 
*/




// Header guard
#ifndef MODULE_COSMODISTANCES_H
#define MODULE_COSMODISTANCES_H




// Include
#include <structure_hpc.h>
#include <string>





// Function definitions
void module_cosmodistances_lensSourceToSource( const int nsetofimages, int nImagesSet[], double z_lens, galaxy image[], double cosmoratio[], cosmo_param cosmopar);
double module_cosmodistances_observerObject(double z, cosmo_param cosmopar);
double module_cosmodistances_objectObject(double z1, double z2, cosmo_param cosmopar);
double module_cosmodistances_lensSourceToObserverSource(double zl, double zs, cosmo_param cosmopar);
int module_cosmodistances_debug(int runmode[], double strongLensingRatios_lensSourceToSource[], double weakLensingRatios_lensSourceToSource[], double weakLensing_observerSource[], int numberCleanLens, double cleanlensRatios_lensSourceToSource[], double cleanlens_observerSource[], std::string DEBUG );

 
#endif // end header guard
