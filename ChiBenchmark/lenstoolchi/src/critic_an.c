#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        critic              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Write the ellipse parameters of the critical and radial lines
 * around the lens[0].
 *
 * Output format :
 * cr_an.dat : 1 C.x C.y aa bb theta
 *
 */

void    critic_an()
{
    const extern struct pot lens[];
    FILE *OUT;
    double aa, bb;
    int i = 0;

    /* ATTENTION: l'ellipticite ici est celle de la distribution de masse   */
    /*  et non pas celle de la densite          */

    OUT = fopen("cr_an.dat", "w");

    bb = lens[i].ct / sqrt(1. - 2.*lens[i].emass);
    aa = lens[i].ct / sqrt(1. + 2.*lens[i].emass);

    fprintf(OUT, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t\n", 1, lens[i].C.x, lens[i].C.y,
            aa, bb, 180. / PI*lens[i].theta);

    aa = lens[i].cr / sqrt(1. - 2. / 9.*lens[i].emass);
    bb = lens[i].cr / sqrt(1. + 2. / 9.*lens[i].emass);

    fprintf(OUT, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t\n", 1, lens[i].C.x, lens[i].C.y,
            aa, bb, 180. / PI*lens[i].theta);

    fclose(OUT);
}

