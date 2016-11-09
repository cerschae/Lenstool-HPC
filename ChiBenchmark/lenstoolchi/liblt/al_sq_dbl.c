/**********************************************************/
/*                                                        */
/* FUNCTION:  alloc_square_double                          */
/*                                                        */
/* PURPOSE: Allocates a square of double in memory         */
/*                                                        */
/* INPUT:  nbr_lin = number of lines                      */
/*         nbr_col = number of columns                    */
/*                                                        */
/* RETURN: square = pointer to matrix of pointers         */
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
/* Allocate an array of double values of size[nbr_lin][nbr_col]
 */
double	**alloc_square_double(int nbr_lin,int nbr_col)
{
	auto     double   **square;
	register int      i, j;
	
	square = (double **) malloc((unsigned) nbr_lin*sizeof(double *));
	if (square != 0)
	{
		for (i=0; i<nbr_lin; i++)
	    {
	    	square[i] = (double *)malloc((unsigned) nbr_col*sizeof(double));
	    	if (square[i] == 0) square = 0;
	    }
	}
	
	for (i=0; i<nbr_lin; i++)
		for (j=0; j<nbr_col; j++)
			square[i][j]=0.0;
	
	return(square);
}
