/**********************************************************/
/*                                                        */
/* FUNCTION: free_cubic_double                            */
/*                                                        */
/* PURPOSE: Frees a cube of double allocated dynamically  */
/* by alloc_cubic_double()                                */
/*                                                        */
/* INPUT: cube = pointer to matrix of pointers            */
/*        nbr_lin = number of lines                       */
/*        nbr_col = number of columns                     */
/*                                                        */
/* VERSION: 1.0    June  2012                             */
/*                                                        */
/* AUTHOR: Karim BOUYOUCEF                                */
/*                                                        */
/**********************************************************/
#include <stdlib.h>
#include "lt.h"

/**********************************************************/
void    free_cubic_double(double ***cube,int nbr_lin, int nbr_col)
{
	register  int    i, j;
	
	for (i=0; i<nbr_lin; i++)
	{
		for(j=0; j<nbr_col; j++)
			free(cube[i][j]);
		free((double *) cube[i]);
	}
	free((double **) cube);
}
