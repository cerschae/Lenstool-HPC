

#ifndef SETUP_HPP_
#define SETUP_HPP_

#include <math.h>
#include <structure_hpc.h>
#include <string.h>
#include <omp.h>

#ifdef __WITH_LENSTOOL
#warning "linking with libtool..."
#include <fonction.h>
#include <constant.h>
#include <dimension.h>
#include <structure.h>
#include <setup.hpp>
#endif

void setup_lenstool();


#endif /* SETUP_HPP_ */
