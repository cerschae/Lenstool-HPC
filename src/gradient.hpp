#ifndef __GRAD_HPP__
#define __GRAD_HPP__
/** for both gradient and second derivatives **/
//inline
//struct point rotateCoordinateSystem(struct point P, double theta);
inline
struct point rotateCoordinateSystem(struct point P, double theta)
{
        struct  point   Q;

        Q.x = P.x*cos(theta) + P.y*sin(theta);
        Q.y = P.y*cos(theta) - P.x*sin(theta);

        return(Q);
}

/** gradient **/
struct point module_potentialDerivatives_totalGradient(const int nhalos, const struct point *pImage, const struct Potential *lens);
//
struct point module_potentialDerivatives_totalGradient_5_SOA(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);
struct point module_potentialDerivatives_totalGradient_8_SOA(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);
struct point module_potentialDerivatives_totalGradient_81_SOA(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);
//
struct point grad_halo(const struct point *pImage, const struct Potential *lens);
//
/** PIEMD **/
complex piemd_1derivatives_ci05(double x, double y, double eps, double rc);

/** Potential **/
/*
void module_readParameters_calculatePotentialparameter(Potential *lens);
void module_readParameters_calculatePotentialparameter_SOA(Potential_SOA *lens, int i);
*/
#endif
