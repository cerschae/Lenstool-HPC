#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        o_dpl               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/*  Return an array of source position corresponding to the given images
 * Parameters :
 * n_im is the number of arclets for the source
 * gali is an array of arclets
 * ps   is the returned array of sources position. !!!ps[n_im]=ps[0]
 ****************************************************************
 *
 * Global variables used :
 * - in e_grad_gal() : G, lens, lens_table
 */
void    o_dpl(int n_im, struct galaxie gali[], struct point ps[], double *np_b0)
{
    struct point    Grad;
    struct galaxie  *pgali;
    double          dlsds;
    register int    i;

    dlsds = gali[0].dr;

    for (i = 0; i < n_im; i++)
    {
        pgali = &gali[i];
        /* the multiple images z is defined only once by the first image*/
        /* but it has to be well defined for each image in the multiple images file */
        Grad = e_grad_gal(pgali, np_b0);
        pgali->Grad = Grad;

        ps[i].x = pgali->C.x - pgali->dr * Grad.x;
        ps[i].y = pgali->C.y - pgali->dr * Grad.y;
    };

    ps[n_im] = ps[0];
    /* printf("dr=%.3lf Ps=%.3lf\n",gali[0].dr,ps[0].x); */

}
