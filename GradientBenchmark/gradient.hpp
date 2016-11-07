#pragma once
#ifndef __GRAD_HPP__
#define __GRAD_HPP__

#include "structure.h"

/** for both gradient and second derivatives **/
static struct point rotateCoordinateSystem(struct point P, double theta);

/** gradient **/
struct point module_potentialDerivatives_totalGradient(const runmode_param *runmode, const struct point *pImage, const struct Potential *lens);
struct point module_potentialDerivatives_totalGradient_SOA(const int *Nlens, const struct point *pImage, const struct Potential_SOA *lens);
struct point module_potentialDerivatives_totalGradient_SOA_AVX(const int *Nlens, const struct point *pImage, const struct Potential_SOA *lens);

static struct point grad_halo(const struct point *pImage, const struct Potential *lens);

struct point grad_halo_sis_SOA(const struct point *pImage, int iterator,Potential_SOA *lens);
struct point grad_halo_piemd_SOA(const struct point *pImage, int iterator,Potential_SOA *lens);

struct point grad_halo_sis_SOA_AVX(const struct point *pImage, int iterator,Potential_SOA *lens);
struct point grad_halo_piemd_SOA_AVX(const struct point *pImage, int iterator,Potential_SOA *lens);

/** PIEMD **/
static complex piemd_1derivatives_ci05(double x, double y, double eps, double rc);

/** Potential **/
void module_readParameters_calculatePotentialparameter(Potential *lens);
void module_readParameters_calculatePotentialparameter_SOA(Potential_SOA *lens, int i);
#endif
