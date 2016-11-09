/**********************************************************/
/*                                                        */
/* FUNCTION:  alloc_square_int                            */
/*                                                        */
/* PURPOSE: Allocates a square of int in memory           */
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
int	**alloc_square_int(int nbr_lin,int nbr_col)
{
auto     int   **square;
register int      i, j;

square = (int **) malloc((unsigned) nbr_lin*sizeof(int *));
if (square != 0)
  {
  for (i=0; i<nbr_lin; i++)
    {
    square[i] = (int *)malloc((unsigned) nbr_col*sizeof(int));
    if (square[i] == 0) square = 0;
    }
  }

for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
square[i][j]=0;

return(square);
}
