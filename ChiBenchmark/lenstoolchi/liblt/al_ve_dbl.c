/**********************************************************/
/*                                                        */
/* FUNCTION:  alloc_vector_double                          */
/*                                                        */
/* PURPOSE: Allocates a vector of double in memory         */
/*                                                        */
/* INPUT:  nbr_lin = number of lines                      */
/*                                                        */
/* RETURN: vector = pointer to matrix of double 	          */
/*                (NULL if memory allocation failure)     */
/*                                                        */
/* VERSION: 1.1  March  1992                              */
/*                                                        */
/* AUTHOR: Karim BOUYOUCEF                                */
/*                                                        */
/**********************************************************/
#include <stdlib.h>
#include "lt.h"

/**********************************************************/
double	*alloc_vector_double(int nbr_lin)
{
	auto double   *vector;
	register int      i;
	
	vector = (double *) malloc((unsigned) nbr_lin*sizeof(double));
	for (i=0; i<nbr_lin; i++)
		vector[i] = 0.0;
	return(vector);
}
