/**********************************************************/
/*                                                        */
/* FUNCTION: free_cubic_int                               */
/*                                                        */
/* PURPOSE: Frees a cube of int allocated dynamically     */
/* by alloc_cubic_int()                                   */
/*                                                        */
/* INPUT: cube = pointer to matrix of pointers            */
/*        nbr_lin = number of lines                       */
/*        nbr_col = number of columns                     */
/*                                                        */
/* VERSION: 1.0    October  2012                          */
/*                                                        */
/* AUTHOR: Johan RICHARD                                  */
/*                                                        */
/**********************************************************/
#include <stdlib.h>
#include "lt.h"

/**********************************************************/
void    free_cubic_int(int ***cube,int nbr_lin, int nbr_col)
{
	register  int    i, j;
	
	for (i=0; i<nbr_lin; i++)
	{
		for(j=0; j<nbr_col; j++)
			free(cube[i][j]);
		free((int *) cube[i]);
	}
	free((int **) cube);
}
