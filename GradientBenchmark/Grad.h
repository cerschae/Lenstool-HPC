// header guard
#ifndef GRAD_H
#define GRAD_H

//include
#include <iostream>
#include <string.h>
#include "structure.h"
#include <math.h>

/** for both gradient and second derivatives **/
struct point rotateCoordinateSystem(struct point P, double theta);

/** gradient **/
struct point module_potentialDerivatives_totalGradient(const int *Nlens, const struct point *pImage, PotentialSet *lens );
struct point grad_halo_sis(const struct point *pImage, int iterator,PotentialSet *lens);
struct point grad_halo_piemd(const struct point *pImage, int iterator,PotentialSet *lens);

/** PIEMD **/
complex piemd_1derivatives_ci05(double x, double y, double eps, double rc);

/** Potential **/
void module_readParameters_calculatePotentialparameter(Potential *lens);

// end header guard
#endif
