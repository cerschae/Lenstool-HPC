#include <structure.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>

// Append a sky background to the pixels of image **z
 
void d_bruiter_omp(double **z, int nx, int ny)
{
   const extern struct g_observ O;
   
   //at this point we mustn't be in parallel
   check_not_in_parallel();
   
   int i;   
#pragma omp parallel for schedule(static)
   for (i = 0; i < ny; i++)
     {	
	int j;	
	for (j = 0; j < nx; j++)
         z[i][j] += O.SKY;
     }
}
