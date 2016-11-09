#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        updatecut           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       12/94               */
/*      place:      IoA Cambridge           */
/****************************************************************
 * Update the parameters psicut, psimcut and psiccut of the
 * type 41 potentials from the optimised parameters b0, rcut and rc.
 * Parameters :
 * - i : potential index in the global lens array.
 */

void    updatecut(int i)
{
    extern  struct  pot     lens[];

    double  rcut2, den;

    switch (lens[i].type)
    {
        case(41):
            rcut2 = lens[i].rcut * lens[i].rcut;
            den = pow(1. + rcut2, 1.5);
            lens[i].psicut = lens[i].b0 / lens[i].rc * (1. + .5 * rcut2) / den / 2.;
            lens[i].psimcut = lens[i].b0 * lens[i].rc * 0.5 * rcut2 * rcut2 / den;
            lens[i].psiccut = lens[i].b0 * lens[i].rc * (1. + 1.5 * rcut2 * (1. + rcut2 / 2.)) / den
                              - lens[i].psimcut * log(lens[i].rcut * lens[i].rc);
            break;

        default:
            break;
    }

}


void updatecut_ptr(struct pot *ilens)
{
    double  rcut2, den;

    switch (ilens->type)
    {
        case(41):
            rcut2 = ilens->rcut * ilens->rcut;
            den = pow(1. + rcut2, 1.5);
            ilens->psicut = ilens->b0 / ilens->rc * (1. + .5 * rcut2) / den / 2.;
            ilens->psimcut = ilens->b0 * ilens->rc * 0.5 * rcut2 * rcut2 / den;
            ilens->psiccut = ilens->b0 * ilens->rc * (1. + 1.5 * rcut2 * (1. + rcut2 / 2.)) / den
                              - ilens->psimcut * log(ilens->rcut * ilens->rc);
            break;

        default:
            break;
    }

}
