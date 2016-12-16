#ifndef __GRADAVX_HPP__
#define __GRADAVX_HPP__
//
/** for both gradient and second derivatives **/
//static struct point rotateCoordinateSystem(struct point P, double theta);
#include "gradient.hpp"

/** gradient **/
struct point module_potentialDerivatives_totalGradient_SOA_AVX(const struct point *pImage, const struct Potential_SOA *lens,const int nhalos);
//
struct point module_potentialDerivatives_totalGradient_81_SOA_AVX(const struct point *pImage, const struct Potential_SOA *lens, const int nhalos);
//
#endif
