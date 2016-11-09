/**********************************************************/
/*                                                        */
/* FUNCTION: free_square_double                            */
/*                                                        */
/* PURPOSE: Frees a square of double allocated dynamically */
/* by alloc_square_double()                                */
/*                                                        */
/* INPUT: square = pointer to matrix of pointers          */
/*        nbr_lin = number of lines                       */
/*                                                        */
/* VERSION: 1.1  March  1992                              */
/*                                                        */
/* AUTHOR: Karim BOUYOUCEF                                */
/*                                                        */
/**********************************************************/
#include <stdlib.h>
#include "lt.h"

/**********************************************************/
void    free_square_double(double **square,int nbr_lin)
{
	register  int    i;
	if(square!=NULL)
	{
	    for (i=0; i<nbr_lin; i++)
		free(square[i]);
		
         	free((double *) square);
        }
}
