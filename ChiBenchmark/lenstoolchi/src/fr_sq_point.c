/**********************************************************/
/*                                                        */
/* FUNCTION: free_square_point                       */
/*                                                        */
/* PURPOSE: Frees a complex square of point allocated     */
/* dynamically by alloc_square_point()               */
/*                                                        */
/* INPUT: square = complex pointer to matrix of pointers  */
/*        nbr_lin = number of lines                       */
/*                                                        */
/* VERSION: 1.1  March  1992                              */
/*                                                        */
/* AUTHOR: Karim BOUYOUCEF                                */
/*                                                        */
/**********************************************************/
#include <stdlib.h>
#include "structure.h"

void    fr_sq_point(struct point **square, int nbr_lin)
{
    register  int    i;

    for (i = 0; i < nbr_lin; i++)
        free((struct point *) square[i]);

    free((struct point *) square);
}
