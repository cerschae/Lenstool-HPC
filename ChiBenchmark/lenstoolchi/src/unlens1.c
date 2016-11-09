#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        unlens1             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

struct  galaxie unlens1(struct galaxie arclet, double dlsds)
{
    struct  galaxie source;
    struct  ellipse ampli;

    ampli = e_unmag_gal(&arclet);
    isoima(&arclet.E, &ampli, &source.E);
    e_dpl(&arclet.C, dlsds, &source.C);

    return(source);
}
