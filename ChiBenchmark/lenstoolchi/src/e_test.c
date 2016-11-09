#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "fonction.h"
#include "constant.h"
#include"dimension.h"
#include "structure.h"

/****************************************************************/
/*      nom:        e_test              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse
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
 * - Tsol : a pointer to the last arclet of a familly (see inverse())
 * - source : a source
 * - image : a list of images for this source
 ****************************************************************/
int e_test( struct bitriplet Tsol[NIMAX], int ni,
            struct galaxie source,
            struct galaxie image[NIMAX])
{
    int it;
    //struct    bitriplet   btriangle;
    struct  bitriplet   btp; /*the smallest bitriangle around the source point*/
    struct  ellipse ampli;
    //struct    chaine *parent;
    double  D;

    int i, j;

    j = 0;

    /*Iterate over every arclet of a familly beginning by the last one*/
    for ( i = 0; i < ni; i++ ) //Tsol != NULL)
    {
        //btriangle.i=Tsol->I;
        //btriangle.s=Tsol->S;
        it = 0;

        /*Find the smallest bitriangle around the source*/
        e_im_prec(&Tsol[i], &source.C, source.dr, &it, &btp);

        image[j].C = barycentre(&btp.i);
        D = dist(source.C, barycentre(&btp.s));

        if (D < 0.02)
            if ( j == 0 || dist(image[j].C, image[j-1].C) > 0.001 )
            {
                strcpy(image[j].n, source.n);
                strcpy(image[j+1].n, "");
                image[j].z = source.z;
                image[j].dos = source.dos;
                image[j].dl0s = source.dl0s;
                image[j].dr = source.dr;
                image[j].c = source.c;
                ampli = e_mag_gal(&image[j]);
                if (image[j].c == 's')
                    image[j].E.a = image[j].E.b = image[j].E.theta = 0.;
                else
                {
                    isoima(&source.E, &ampli, &image[j].E);
                    image[j].A = fabs(ampli.a * ampli.b);
                    image[j].mag = source.mag - 2.5 * log10(fabs(ampli.a * ampli.b));
                }

                j++;
            }

        /*      parent = Tsol;
                Tsol=(Tsol->F);
                free(parent);*/
    }

    return(j);
}

/**************************************************************/
/* Return a triplet of points in the source plane corresponding to the triplet
 * of images. dlsds is the lens efficiency at the source redshift.
 *
 * Global variables used :
 * - in e_dpl() : G, lens, lens_table
 */
void e_transform(struct triplet *I, double dlsds, struct triplet *S)
{
    e_dpl(&I->a, dlsds, &S->a);
    e_dpl(&I->b, dlsds, &S->b);
    e_dpl(&I->c, dlsds, &S->c);
}
/**************************************************************/
/* Return the triangle inside the triangle E. Each corner is located
 * at equal distance of the opposite corners.
 *
 * Ia is in front of Ea, Ib is in front of Eb and Ic is in front of Ec.
 *
 *          Ea
 *         /  \
 *       Ic -- Ib
 *       / \  / \
 *     Eb---Ia---Ec
 *
 * Global variables used :
 * - none
 */
void interieur(const struct triplet *E, struct triplet *I)
{
    I->a.x = (E->b.x + E->c.x) / 2.;
    I->a.y = (E->b.y + E->c.y) / 2.;
    I->b.x = (E->c.x + E->a.x) / 2.;
    I->b.y = (E->c.y + E->a.y) / 2.;
    I->c.x = (E->a.x + E->b.x) / 2.;
    I->c.y = (E->a.y + E->b.y) / 2.;
}


/*********************************************************************
 * Copy one triplet to another one.
 *
 * Global variables used :
 * - none
 */
void copyTriplet(struct triplet *in, struct triplet *out)
{
    out->c = in->a;
    out->a = in->b;
    out->b = in->c;
}
/*********************************************************************/

