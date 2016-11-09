#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "fonction.h"
#include "constant.h"
#include"dimension.h"
#include "structure.h"

/****************************************************************/
/*      nom:        e_test_P            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Test if it is possible to find 2 small triangles that contains the
 * source in the source plane and one arclet found previously in the
 * image plane.
 *
 * Fill the image list with the position of the arclets corresponding
 * to the source.
 *
 * The arclets must have their source close enough from
 * the effective source position and far enough from the previously
 * found arclet of the same familly.
 *
 * Return the number of arclet for this source.
 *
 * Parameters :
 * - Tsol : a list of arclets of a familly (see inverse())
 * - ni : number of arclets in Tsol
 * - ps : source position
 * - image : a list of images for this source
 * - dlsds : DLS/DS ratio for the source
 ****************************************************************/


int e_test_P(struct bitriplet Tsol[NIMAX], int ni,
             struct point *ps,
             struct point image[NIMAX],
             double dlsds, double err)
{
    int it;
    //struct    bitriplet   btriangle;
    struct  bitriplet   btp; /*the smallest bitriangle around the source point*/
    //struct    chaine *parent;
    double  D;
    int i, j;

    j = 0;

    for ( i = 0 ; i < ni ; i++ )
    {
        //btriangle.i=(Tsol->I);
        //btriangle.s=(Tsol->S);
        it = 0;
        e_im_prec(&Tsol[i], ps, dlsds, &it, &btp);

        image[j] = barycentre(&btp.i);
        D = dist2(barycentre(&btp.s), (*ps));
        /*
            if ((D<err)&&((j==0)||(dd=dist2(image[j],image[j-1])>err)))
        98/09/16: D test removed, this was checking the resolution in the source
        plane which is not important for this routine.
        06/10/12: D test introduced again because of merging with the e_testclean() function
        */
        if ((D < err) && ((j == 0) || (dist2(image[j], image[j-1]) > err)))
            j++;

        /*parent = Tsol;
        Tsol=(Tsol->F);
        free(parent);*/
    };

    return(j);
}
