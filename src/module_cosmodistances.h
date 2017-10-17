
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
void module_cosmodistances_lensSourceToSource( const int nsetofimages, int nImagesSet[], type_t z_lens, galaxy image[], type_t cosmoratio[], cosmo_param cosmopar);
type_t module_cosmodistances_observerObject(type_t z, cosmo_param cosmopar);
type_t module_cosmodistances_objectObject(type_t z1, type_t z2, cosmo_param cosmopar);
type_t module_cosmodistances_lensSourceToObserverSource(type_t zl, type_t zs, cosmo_param cosmopar);
int module_cosmodistances_debug(int runmode[], type_t strongLensingRatios_lensSourceToSource[], type_t weakLensingRatios_lensSourceToSource[], type_t weakLensing_observerSource[], int numberCleanLens, type_t cleanlensRatios_lensSourceToSource[], type_t cleanlens_observerSource[], std::string DEBUG );

 
#endif // end header guard
