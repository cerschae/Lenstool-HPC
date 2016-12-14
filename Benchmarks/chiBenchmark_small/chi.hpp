#pragma once
#ifndef __CHI_HPP__
#define __CHI_HPP__

#include "structure.h"
#include "gradient.hpp"
//#ifdef __AVX512F__
#include "gradient_avx512f.hpp"
//#endif

void chi_bruteforce(double *chi, int *error, runmode_param *runmode, const struct Potential *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images);
void chi_bruteforce_SOA(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images);
void chi_bruteforce_SOA_AVX(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images);


void chi_transformImageToSourcePlane(const runmode_param *runmode, const struct point *image_point, double dlsds, const struct Potential *lens, struct point *source_point);
void chi_transformImageToSourcePlane_SOA(const int *Nlens, const struct point *image_point, double dlsds, const struct Potential_SOA *lens, struct point *source_point);
void chi_transformImageToSourcePlane_SOA_AVX(const int *Nlens, const struct point *image_point, double dlsds, const struct Potential_SOA *lens, struct point *source_point);

void chi_transformtriangleImageToSourcePlane(const runmode_param *runmode, struct triplet *I, double dlsds, const struct Potential *lens, struct triplet *S);
void chi_transformtriangleImageToSourcePlane_SOA(const int *Nlens, struct triplet *I, double dlsds, const struct Potential_SOA *lens, struct triplet *S);
void chi_transformtriangleImageToSourcePlane_SOA_AVX(const int *Nlens, struct triplet *I, double dlsds, const struct Potential_SOA *lens, struct triplet *S);

double chi_determinant(const struct point *A,const struct point *B,const struct point *C);
int chi_inside(const struct point *P, struct triplet *T);
int chi_insideborder(const struct point *P, struct triplet *T);
struct  point   chi_barycenter(struct triplet *A);
double chi_dist(struct point A, struct point B);
#endif
