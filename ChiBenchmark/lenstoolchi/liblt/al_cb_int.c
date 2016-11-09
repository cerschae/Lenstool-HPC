/**********************************************************/
/*                                                        */
/* FUNCTION:  alloc_cubic_int                          */
/*                                                        */
/* PURPOSE: Allocates a cube of integer in memory         */
/*                                                        */
/* INPUT:  nbr_lin = number of lines                      */
/*         nbr_col = number of columns                    */
/*         nbr_slice = number of slices                   */
/*                                                        */
/* RETURN: cube = pointer to matrix of pointers           */
/*                of matix of pointers                    */
/*                (NULL if memory allocation failure)     */
/*                                                        */
/* VERSION: 1.0    October 2012                           */
/*                                                        */
/* AUTHOR: Johan RICHARD                                  */
/*                                                        */
/**********************************************************/
#include <stdlib.h>
#include "lt.h"

/**********************************************************/
/* Allocate an array of int values of size[nbr_lin][nbr_col][nbr_slices]
 */
int	***alloc_cubic_int(int nbr_lin,int nbr_col,int nbr_slice)
{
	auto     int   ***cube;
	register int      i, j, k;
	
	cube = (int ***) malloc((unsigned) nbr_lin*sizeof(int **));
	if (cube != 0)
	{
		for (i=0; i<nbr_lin; i++)
	    {
		cube[i] = (int **) malloc((unsigned) nbr_col*sizeof(int *));
		if (cube !=0)
		{
			for (j=0;j<nbr_col;j++)
			{
				cube[i][j] = (int *) malloc((unsigned) nbr_slice*sizeof(int));
				if (cube[i][j] ==0) cube=0;
			}
		}
	    }
	}

	for (i=0; i<nbr_lin; i++)
		for (j=0; j<nbr_col; j++)
			for(k=0;k<nbr_slice;k++)
				cube[i][j][k]=0;
	
	return(cube);
}
