#include<math.h>
#include<structure.h>

/****************************************************************/
/*      nom:        update_epot     */
/*      auteur:     Eric Jullo          */
/*      date:       4/2007              */
/*      place:  ESO, Chile              */
/****************************************************************/

void update_epot(int i, double *epot)
{
    const extern struct pot lens[];

    if ( lens[i].type ==  8  ||
         lens[i].type == -2  ||
         ( lens[i].type > 80 && lens[i].type < 90 ) )
        *epot = (1. - sqrt(1 - lens[i].emass * lens[i].emass)) / lens[i].emass;
    else
        *epot = lens[i].emass / 3.;

}

void update_epot_ptr(struct pot *ilens, double *epot)
{
    if ( ilens->type ==  8  ||
         ilens->type == -2  ||
         ( ilens->type > 80 && ilens->type < 90 ) )
        *epot = (1. - sqrt(1 - ilens->emass * ilens->emass)) / ilens->emass;
    else
        *epot = ilens->emass / 3.;
}
