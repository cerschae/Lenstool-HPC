/**********************************************************/
/*                                                        */
/* FUNCTION:  alloc_cubic_double                          */
/*                                                        */
/* PURPOSE: Allocates a cube of double in memory          */
/*                                                        */
/* INPUT:  nbr_lin = number of lines                      */
/*         nbr_col = number of columns                    */
/*         nbr_slice = number of slices                   */
/*                                                        */
/* RETURN: cube = pointer to matrix of pointers           */
/*                of matix of pointers                    */
/*                (NULL if memory allocation failure)     */
/*                                                        */
/* VERSION: 1.0    June 2012                              */
/*                                                        */
/* AUTHOR: Vincent BINET                                  */
/*                                                        */
/**********************************************************/
#include <stdlib.h>
#include "lt.h"

/**********************************************************/
/* Allocate an array of double values of size[nbr_lin][nbr_col][nbr_slices]
 */
double	***alloc_cubic_double(int nbr_lin,int nbr_col,int nbr_slice)
{
	register int      i, j, k;
	
	double ***cube = (double ***) malloc((unsigned) nbr_lin*sizeof(double **));
	if (cube != 0)
	{
		for (i=0; i<nbr_lin; i++)
	    {
		cube[i] = (double **) malloc((unsigned) nbr_col*sizeof(double *));
		if (cube !=0)
		{
			for (j=0;j<nbr_col;j++)
			{
				cube[i][j] = (double *) malloc((unsigned) nbr_slice*sizeof(double));
				if (cube[i][j] ==0) cube=0;
			}
		}
	    }
	}

	for (i=0; i<nbr_lin; i++)
		for (j=0; j<nbr_col; j++)
			for(k=0;k<nbr_slice;k++)
				cube[i][j][k]=0.0;

	return((double ***)cube);
}
