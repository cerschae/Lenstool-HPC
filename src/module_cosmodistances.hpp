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

@brief: Function for cosmological distance computations

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
