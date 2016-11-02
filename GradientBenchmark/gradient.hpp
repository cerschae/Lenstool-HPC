#pragma once
#ifndef __GRAD_HPP__
#define __GRAD_HPP__
/** for both gradient and second derivatives **/
static struct point rotateCoordinateSystem(struct point P, double theta);

/** gradient **/
struct point module_potentialDerivatives_totalGradient(const runmode_param *runmode, const struct point *pImage, const struct Potential *lens);
struct point module_potentialDerivatives_totalGradient_SOA(const runmode_param *runmode, const struct point *pImage, const struct Potential_SOA *lens, int nhalos);
struct point module_potentialDerivatives_totalGradient_SOA_AVX(const runmode_param *runmode, const struct point *pImage, const struct Potential_SOA *lens, int nhalos);
//struct point potentialDerivatives_totalGradient(const runmode_param *runmode, const struct point *pImage, const struct Potential *lens, int nhalos);
//struct point potentialDerivatives_totalGradient_SOA(const runmode_param *runmode, const struct point *pImage, const struct Potential *lens, int nhalos);
//
static struct point grad_halo(const struct point *pImage, const struct Potential *lens);

/** PIEMD **/
static complex piemd_1derivatives_ci05(double x, double y, double eps, double rc);

/** Potential **/
void module_readParameters_calculatePotentialparameter(Potential *lens);
void module_readParameters_calculatePotentialparameter_SOA(Potential_SOA *lens, int i);
#endif
