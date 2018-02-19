/**
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)
* @date   July 2017
* @version 0,1
* 
*/



// Header guard
#ifndef MODULE_COSMODISTANCES_H
#define MODULE_COSMODISTANCES_H




// Include
#include <structure_hpc.hpp>
#include <string>





// Function definitions
void module_cosmodistances_lensSourceToSource( const int nsetofimages, int nImagesSet[], type_t z_lens, galaxy image[], type_t cosmoratio[], cosmo_param cosmopar);
type_t module_cosmodistances_observerObject(type_t z, cosmo_param cosmopar);
type_t module_cosmodistances_objectObject(type_t z1, type_t z2, cosmo_param cosmopar);
type_t module_cosmodistances_lensSourceToObserverSource(type_t zl, type_t zs, cosmo_param cosmopar);
int module_cosmodistances_debug(int runmode[], type_t strongLensingRatios_lensSourceToSource[], type_t weakLensingRatios_lensSourceToSource[], type_t weakLensing_observerSource[], int numberCleanLens, type_t cleanlensRatios_lensSourceToSource[], type_t cleanlens_observerSource[], std::string DEBUG );

void module_cosmodistances_relativecoordinates_XY( type_t *x, type_t *y, int iref, type_t ref_ra, type_t ref_dec );

 
#endif // end header guard
