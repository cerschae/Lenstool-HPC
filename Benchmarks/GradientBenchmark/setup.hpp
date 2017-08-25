#pragma once
#ifndef __SETUP__
//
void
setup_jauzac_SOA(Potential_SOA *lens_soa, int* nlenses, type_t* x, type_t* y, type_t* sol_grad_x, type_t* sol_grad_y);
//
void
setup_jauzac(Potential** lens, int* nlenses, type_t* x, type_t* y, type_t* sol_grad_x, type_t* sol_grad_y);
void
//setup_jauzac_LT(pot** lens, int* nlenses, type_t* x, type_t* y, type_t* sol_grad_x, type_t* sol_grad_y);
setup_jauzac_LT(int* nlenses, type_t* x, type_t* y, type_t* sol_grad_x, type_t* sol_grad_y);

//
#endif
