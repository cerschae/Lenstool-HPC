#pragma once
#ifndef __CHI_HPP__
#define __CHI_HPP__

#include "structure.h"
#include "gradient.hpp"

void chi_bruteforce(double *chi, int *error, runmode_param *runmode, const struct Potential *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images);
void module_chiClassic_transformImageToSourcePlane(const runmode_param *runmode, const struct point *image_point, double dlsds, const struct Potential *lens, struct point *source_point);
void module_chiClassic_transformtriangleImageToSourcePlane(const runmode_param *runmode, struct triplet *I, double dlsds, const struct Potential *lens, struct triplet *S);
double module_chiClassic_determinant(const struct point *A,const struct point *B,const struct point *C);
int module_chiClassic_inside(const struct point *P, struct triplet *T);
int module_chiClassic_insideborder(const struct point *P, struct triplet *T);
struct  point   module_chiClassic_barycenter(struct triplet *A);
double module_chiClassic_dist(struct point A, struct point B);
#endif
