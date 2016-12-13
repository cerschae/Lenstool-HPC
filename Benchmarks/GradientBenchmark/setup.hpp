#pragma once
#ifndef __SETUP__
//
void
setup_jauzac_SOA(Potential_SOA *lens_soa, int* nlenses, double* x, double* y, double* sol_grad_x, double* sol_grad_y);
//
void
setup_jauzac(Potential** lens, int* nlenses, double* x, double* y, double* sol_grad_x, double* sol_grad_y);
//
#endif
