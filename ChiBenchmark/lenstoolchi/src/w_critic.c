#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        w_critic            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Write the radial and tangeant critical lines in the ce.dat (radial)
 * and ci.dat (tangeantial) ASCII files.
 *
 * Output format :
 *  - ce.dat : Id Image.x Image.y Source.x Source.y
 *  - ci.dat : Id Image.x Image.y Source.x Source.y
 *
 */

void    w_critic()
{
    const extern  struct  g_mode  M;
    const extern  struct  biline  radial[], tangent[];
    const extern  int nrline, ntline;

    FILE    *OUT;
    int i;

    OUT = fopen("ce.dat", "w");
    // Append the reference point at the beginning of the files
    fprintf(OUT, "#REFERENCE 3 %.7lf %.7lf\n", M.ref_ra, M.ref_dec);

    for (i = 0; i < nrline; i++)
    {
        fprintf(OUT, "%d\t%lf\t%lf\t%lf\t%lf\t\n", radial[i].i,
                radial[i].I.x, radial[i].I.y, radial[i].S.x, radial[i].S.y);
    };
    fclose(OUT);

    OUT = fopen("ci.dat", "w");
    // Append the reference point at the beginning of the files
    fprintf(OUT, "#REFERENCE 3 %.7lf %.7lf\n", M.ref_ra, M.ref_dec);

    for (i = 0; i < ntline; i++)
    {
        fprintf(OUT, "%d\t%lf\t%lf\t%lf\t%lf\t\n", tangent[i].i,
                tangent[i].I.x, tangent[i].I.y, tangent[i].S.x, tangent[i].S.y);
    };
    fclose(OUT);
}
