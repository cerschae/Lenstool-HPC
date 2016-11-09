/***************************************************************/
/*                                                             */
/* FUNCTION: split_image                                       */
/*                                                             */
/* PURPOSE: splits an image                        */
/*                                                             */
/* INPUT:  image    image to split                         */
/*         nbr_lin  number of lines of image                   */
/*         nbr_col  number of colums of image                  */
/*                                                             */
/* RETURN: 0 if everything is OK                               */
/*         -1 if memory allocation failure                     */
/*                                                             */
/* VERSION: December 1992                                      */
/*                                                             */
/* AUTHOR: Karim BOUYOUCEF                                     */
/*                                                             */
/***************************************************************/
#include <stdio.h>
#include "fonction.h"
#include "lt.h"

int  split_image(double **image, int nbr_lin, int nbr_col)
{
    /****************  declarations  ************************************/
    auto int     i, j;
    auto double  tmp;

    /***********  split of image  ********************/
    for (i = 0; i < nbr_lin / 2; i++)
        for (j = 0; j < nbr_col / 2; j++)
        {
            tmp = image[i+nbr_lin/2][j+nbr_col/2];
            image[i+nbr_lin/2][j+nbr_col/2] = image[i][j];
            image[i][j] = tmp;

            tmp = image[i][j+nbr_col/2];
            image[i][j+nbr_col/2] = image[i+nbr_lin/2][j];
            image[i+nbr_lin/2][j] = tmp;
        }

    return(0);
}
