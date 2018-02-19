/**
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)
* @date   July 2017
* @version 0,1
*
*/
#ifndef __GRAD512F_HPP__
#define __GRAD512F_HPP__
#ifdef __AVX512F__
/** for both gradient and second derivatives **/
static struct point rotateCoordinateSystem(struct point P, double theta);

/** gradient **/
struct point module_potentialDerivatives_totalGradient_8_SOA_AVX512(const struct point *pImage, const struct Potential_SOA *lens,const int nhalos);
struct point module_potentialDerivatives_totalGradient_81_SOA_AVX512(const struct point *pImage, const struct Potential_SOA *lens,const int nhalos);
//
/** PIEMD **/
static complex piemd_1derivatives_ci05(double x, double y, double eps, double rc);
//
/** Potential **/
/*
void module_readParameters_calculatePotentialparameter(Potential *lens);
void module_readParameters_calculatePotentialparameter_SOA(Potential_SOA *lens, int i);
*/
#endif
#endif
