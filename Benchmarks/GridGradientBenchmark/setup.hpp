#include <iostream>
#include <string.h>
//#include <cuda_runtime.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <immintrin.h>

//#include "simd_math.h"

#ifdef __WITH_LENSTOOL
#include "structure.h"
#include "structure_hpc.h"

extern struct pot lens[NLMAX];
void
//setup_jauzac_LT(struct pot** lens, int* nlenses, double* x, double* y, double* sol_grad_x, double* sol_grad_y)
convert_to_LT(Potential_SOA* lenses, int nlenses)
{
        //
        //*lens = (struct pot*) malloc(sizeof(struct pot)*(*nlenses));
        //
        for (int i = 0; i < nlenses; ++i)
        {
                //
                lens[i].C.x            	      = lenses->position_x[i];
                lens[i].C.y                   = lenses->position_y[i];
                //
                lens[i].sigma                 = 1.; 			// disp
                lens[i].type                  = lenses->type[i];
                lens[i].emass                 = 0.11;
                lens[i].epot                  = lenses->ellipticity_potential[i];
                lens[i].theta                 = lenses->theta[i];
                lens[i].rcut                  = lenses->rcut[i];
                lens[i].rc                    = lenses->rcore[i];
                lens[i].b0                    = lenses->b0[i];
                lens[i].masse                 = 0;			// weight
                //lens[i].rc  	                 = 0;			// rscale
//               (&lens)[i].exponent              = 0;
                lens[i].alpha                 = 0.;
//                (&lens)[i].einasto_kappacritic   = 0;
                lens[i].z                     = 0.4;
        }
        std::cout << "Setup done..." << std::endl;
}
#endif








