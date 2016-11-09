/**********************************************************/
/*                                                        */
/* FUNCTION:  alloc_square_point                          */
/*                                                        */
/* PURPOSE: Allocates a square of points in memory        */
/*                                                        */
/* INPUT:  nbr_lin = number of lines                      */
/*         nbr_col = number of columns                    */
/*                                                        */
/* RETURN: square = pointer to matrix of points           */
/*                (NULL if memory allocation failure)     */
/*                                                        */
/* VERSION: 1.1  March  1992                              */
/*                                                        */
/* AUTHOR: JP KNEIB                                   */
/*                                                        */
/**********************************************************/
#include <stdlib.h>
#include "structure.h"

struct   point  **al_sq_point(int nbr_lin, int nbr_col)
{
    struct   point   **square;
    register int      i, j;

    square = (struct point **) malloc((unsigned) nbr_lin
                                      *sizeof(struct point *));

    if (square != 0)
    {
        for (i = 0; i < nbr_lin; i++)
        {
            square[i] = (struct point *)malloc((unsigned) nbr_col
                                               * sizeof(struct point ));
            if (square[i] == 0) square = 0;
        }
    }

    for (i = 0; i < nbr_lin; i++)
        for (j = 0; j < nbr_col; j++)
            square[i][j].x = square[i][j].y = 0.0;

    return( square);
}
