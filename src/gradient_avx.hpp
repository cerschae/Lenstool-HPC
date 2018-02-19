/**
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)
* @date   July 2017
* @version 0,1
*
*/
#ifndef __GRADAVX_HPP__
#define __GRADAVX_HPP__
//
/** for both gradient and second derivatives **/
//static struct point rotateCoordinateSystem(struct point P, double theta);
#include "gradient.hpp"

/** gradient **/
struct point module_potentialDerivatives_totalGradient_8_SOA_AVX(const struct point *pImage, const struct Potential_SOA *lens, const int shalos, const int nhalos);
//
struct point module_potentialDerivatives_totalGradient_81_SOA_AVX(const struct point *pImage, const struct Potential_SOA *lens, const int shalos, const int nhalos);
//
struct point module_potentialDerivatives_totalGradient_SOA_AVX(const struct point *pImage, const struct Potential_SOA *lens, int nhalos);
//
struct point module_potentialDerivatives_totalGradient_SOA_novec(const struct point *pImage, const struct Potential_SOA *lens, int nhalos);
//
#endif
