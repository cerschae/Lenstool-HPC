#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "fonction.h"
#include "constant.h"
#include"dimension.h"
#include "structure.h"


/****************************************************************/
/*      nom:        e_testg             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * For all the arclets of the Tsol list, test if it's possible to find a small
 * bitriangle close to P in the source plane or far enough from the previous
 * arclet of the same familly in the image plane.
 * If so, set the gianti[i][arclet].n identifier to "1".
 *
 * The identifier n set to "" mark the last arclet of the familly.
 *
 * The NULL pointer is returned for Tsol.
 *
 * Parameters :
 * - i : current total number of giant arclets
 * - j : current number of giant arclets for this familly
 * - Tsol : list of arclets (in the source and image planes)
 * - P : point in source plane
 * - dlsds : Lens efficiency for this source*/
void    e_testg(int i, int j,
                struct bitriplet *Tsol, int ni,
                struct point *P,
                double dlsds)
{
    int it;
    extern  struct  pointgal    gianti[][NIMAX];
    //struct    bitriplet btriangle;
    struct  bitriplet btp; /*btp the smallest 2 triangles around the source and 1 arclet*/
    //struct    chaine *parent;
    double  D;

    for ( i = 0; i < ni; i++ ) //while(Tsol!=NULL)
    {
        //btriangle.i=(Tsol->I);
        //btriangle.s=(Tsol->S);
        it = 0;
        e_im_prec(&Tsol[i], P, dlsds, &it, &btp);
        D = dist(barycentre(&btp.s), (*P));

        if (D < 0.002)
        {
            if (j == 0)
            {
                gianti[i][j].C = barycentre(&btp.i);
                gianti[i][j].n = 1;
                j++;
                gianti[i][j].n = 0;
            }
            else
            {
                gianti[i][j].C = barycentre(&btp.i);

                if (dist(gianti[i][j].C, gianti[i][j-1].C) > 0.01)
                {
                    gianti[i][j].n = 1;
                    j++;
                    gianti[i][j].n = 0;
                };
            };
        };
        /*
                parent = Tsol;
                Tsol=(Tsol->F);
                free(parent);*/
    };
}
