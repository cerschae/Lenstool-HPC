#pragma once
#ifndef __CHI_CPU_HPP__
#define __CHI_CPU_HPP__

#include "structure_hpc.hpp"
#include "gradient_avx.hpp"
//#include "grid_srcplane_conversion.hpp"
//#include "grid_gradient_CPU.hpp"
#ifdef __AVX512F__
#include "gradient_avx512f.hpp"
#endif

void mychi_bruteforce_SOA_CPU_grid_gradient(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images);
void mychi_bruteforce_SOA_CPU_grid_gradient_orig(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images);
void mychi_bruteforce_SOA_CPU_grid_gradient_barycentersource(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images);

void mychi_transformImageToSourcePlane(const runmode_param *runmode, const struct point *image_point, double dlsds, const struct Potential *lens, struct point *source_point);
void mychi_transformImageToSourcePlane_SOA(const int Nlens, const struct point *image_point, double dlsds, const struct Potential_SOA *lens, struct point *source_point);
void mychi_transformImageToSourcePlane_SOA_AVX(const int *Nlens, const struct point *image_point, double dlsds, const struct Potential_SOA *lens, struct point *source_point);
void mychi_transformImageToSourcePlane_SOA_Packed(const int *Nlens, const struct point *image_point, double dlsds, struct point *source_point, double *grad_x, double * grad_y, int grad_id);

void mychi_transformtriangleImageToSourcePlane(const runmode_param *runmode, struct triplet *I, double dlsds, const struct Potential *lens, struct triplet *S);
void mychi_transformtriangleImageToSourcePlane_SOA( struct triplet *I, double dlsds, const struct Potential_SOA *lens, struct triplet *S);
void mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_upper( struct triplet *I, double dlsds, struct triplet *S, double *grad_x, double * grad_y, int grad_id, int nbgridcell);
void mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_lower( struct triplet *I, double dlsds, struct triplet *S, double *grad_x, double * grad_y, int grad_id, int nbgridcell);

double mychi_determinant(const struct point *A,const struct point *B,const struct point *C);
int mychi_inside(const struct point *P, struct triplet *T);
int mychi_insideborder(const struct point *P, struct triplet *T);
struct  point   mychi_barycenter(struct triplet *A);
double mychi_dist(struct point A, struct point B);
#endif
