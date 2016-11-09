/**********************************************************/
/*                                                        */
/* FUNCTION: free_square_int                              */
/*                                                        */
/* PURPOSE: Frees a square of int allocated dynamically   */
/* by alloc_square_int()                                  */
/*                                                        */
/* INPUT: square = pointer to matrix of pointers          */
/*        nbr_lin = number of lines                       */
/*                                                        */
/* VERSION: 1.1  March  1992                              */
/*                                                        */
/* AUTHOR: Karim  BOUYOUCEF                               */
/*                                                        */
/**********************************************************/
#include <stdlib.h>
#include "lt.h"

/**********************************************************/
void    free_square_int(int **square,int nbr_lin)
{
	register  int    i;
	
	for (i=0; i<nbr_lin; i++)
		free((int *) square[i]);
		
	free((int *) square);
}
