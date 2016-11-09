#include<fonction.h>
#include<structure.h>

/****************************************************************/
/*      nom:        e_dpl               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

/* Return the position of the source from the position of an image
 *
 * Global variables used :
 * - in e_grad() : G, lens, lens_table
 */
void e_dpl(const struct point *gi, double dlsds, struct point *gs)
{
    struct point Grad;

    Grad = e_grad(gi);
    gs->x = gi->x - dlsds * Grad.x;
    gs->y = gi->y - dlsds * Grad.y;
}
