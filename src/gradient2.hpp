#ifndef __GRAD2_HPP__
#define __GRAD2_HPP__
//inline
//struct point rotateCoordinateSystem(struct point P, double theta);

/** gradient **/

struct matrix module_potentialDerivatives_totalGradient2_SOA(const struct point *pImage, const struct Potential_SOA *lens, int nhalos);

#endif
