#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        e_time              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Return the time delay at point <pi>
 */
double  e_time(struct point pi, double dlsds)
{
    const extern  struct  pot lens[];
    const extern  struct  g_cosmo C;

    double  cst, dl, time;
    struct  point   gr;

    dl = distcosmo1(lens[0].z);
    cst = (1 + lens[0].z) * dl / dlsds * th_a2_day / C.h;   /* time delay in days */

    gr = e_grad(&pi);
    gr.x *= dlsds;
    gr.y *= dlsds;
    time = cst * ((gr.x * gr.x + gr.y * gr.y) / 2. - e_pot(pi, dlsds));

    return(time);
}
