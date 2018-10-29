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

@brief: Function for single order computation on CPUs

*/


#ifndef __GRAD_HPP__
#define __GRAD_HPP__
/** for both gradient and second derivatives **/
//inline
//struct point rotateCoordinateSystem(struct point P, double theta);
inline
struct point rotateCoordinateSystem(struct point P, type_t theta)
{
        struct  point   Q;

        Q.x = P.x*cos(theta) + P.y*sin(theta);
        Q.y = P.y*cos(theta) - P.x*sin(theta);

        return(Q);
}


/** gradient **/
struct point module_potentialDerivatives_totalGradient(const int nhalos, const struct point *pImage, const struct Potential *lens);
//
struct point module_potentialDerivatives_totalGradient_SOA(const struct point *pImage, const struct Potential_SOA *lens, int nhalos);
struct point module_potentialDerivatives_totalGradient_5_SOA(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);
struct point module_potentialDerivatives_totalGradient_8_SOA(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);
struct point module_potentialDerivatives_totalGradient_81_SOA(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);

//struct point module_potentialDerivatives_totalGradient_5_SOA_print(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos, int index);
//
struct point grad_halo(const struct point *pImage, const struct Potential *lens);
//
/** PIEMD **/
complex piemd_1derivatives_ci05(type_t x, type_t y, type_t eps, type_t rc);
//
//struct point (*halo_func[100])(const struct point *pImage, const struct Potential_SOA *lens, int nhalos);

/** Potential **/
/*
void module_readParameters_calculatePotentialparameter(Potential *lens);
void module_readParameters_calculatePotentialparameter_SOA(Potential_SOA *lens, int i);
*/


#define KERNEL(RCUT, ZRES) \
znum.re = cx1*x; \
znum.im = 2.*sqe*sqrt(RCUT*RCUT + rem2) - y/cx1; \
zden.im = 2.*RCUT*sqe - y; \
norm    = (x*x + zden.im*zden.im); \
zis.re  = (znum.re*x + znum.im*zden.im)/norm; \
zis.im  = (znum.im*x - znum.re*zden.im)/norm; \
norm    = zis.re; \
zis.re  = log(sqrt(norm*norm + zis.im*zis.im)); \
zis.im  = atan2(zis.im, norm); \
ZRES.re = - zci.im*zis.im; \
ZRES.im =   zci.im*zis.re; \
\

#endif
