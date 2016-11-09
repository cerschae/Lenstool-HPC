//omp_funs.c
//functions usefull for OpenMP paralellisation

#ifdef _OPENMP
#include "omp.h"
#endif
#include <stdio.h>

void check_not_in_parallel(const char* s)
{
#ifdef _OPENMP
   if (omp_in_parallel())
     {
	fprintf(stderr, "Error | omp_funs.c/check_not_in_parallel | %s", s);
     }
#endif
}

