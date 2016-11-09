#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        study           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    cor_seeing(long int n, struct galaxie *image, double seeing)
{
    long int i;

    for (i = 0; i < n; i++)
    {
        image[i].E.a = sqrt( image[i].E.a * image[i].E.a - seeing * seeing / 4.);
        image[i].E.b = sqrt( image[i].E.b * image[i].E.b - seeing * seeing / 4.);
    }
}
