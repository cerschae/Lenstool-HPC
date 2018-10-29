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

@brief: Handvectorized gradient function for AVX512

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
