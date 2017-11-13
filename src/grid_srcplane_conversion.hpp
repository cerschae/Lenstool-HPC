/*
 * grid_srcplane_conversion.hpp
 *
 *  Created on: Dec 21, 2016
 *      Author: cerschae
 */

#ifndef GRID_SRCPLANE_CONVERSION_HPP_
#define GRID_SRCPLANE_CONVERSION_HPP_

#include <math.h>
#include <structure_hpc.hpp>


void grid_srcplane_conversion_CPU(double *grid_srcplane_x, double *grid_srcplane_y, const struct grid_param *frame, const struct Potential_SOA *lens,const double dlsds, int *Nlens, int nbgridcells);

void grid_srcplane_conversion_sis_CPU(double *grid_srcplane_x, double *grid_srcplane_y, const struct grid_param *frame, int Nlens, int nbgridcells,const double dlsds, const struct Potential_SOA *lens);
void grid_srcplane_conversion_piemd_CPU(double *grid_srcplane_x, double *grid_srcplane_y, const struct grid_param *frame, int Nlens, int nbgridcells,const double dlsds, const struct Potential_SOA *lens);

static struct point rotateCoordinateSystem(struct point P, double theta);


#endif /* GRID_SRCPLANE_CONVERSION_HPP_ */
