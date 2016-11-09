#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        trace               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/*                              */
/*   Modified :                         */
/*      EJ (30/08/2005)
 *
 * Write the pot.dat file on disk.
 ****************************************************************
 * Fill the file pot.dat.
 *
 * For each clump, write :
 * - 1 ellipse with the ct
 */

void tracepot()
{
    const extern struct pot      lens[];
    const extern struct g_grille G;
    const extern struct g_mode   M;

    FILE *OUT;
    int i;
    double aa, bb, Cx, Cy;

    NPRINTF(stderr, "INFO: Write the pot.dat file\n");

    OUT = fopen("pot.dat", "w");
    fprintf( OUT, "#REFERENCE 0 %.7f %.7f\n", M.ref_ra, M.ref_dec );

    for (i = 0; i < G.nlens; i++)
    {
        Cx = lens[i].C.x;
        Cy = lens[i].C.y;
        // define the center of the clump in WCS if defined
        if ( M.iref != 0 )
        {
            Cy = Cy / 3600 + M.ref_dec;
            Cx = -Cx / 3600 / cos(M.ref_dec * DTR) + M.ref_ra;
        }

        // create an ellipse for rcut or Re if PIEMD
        // ct is defined in set_lens.c ( = b0*dlsds)
        if (lens[i].ct == 0 && lens[i].type<81 && lens[i].type>89)
        {
            aa = lens[i].rcut / sqrt(1. - lens[i].emass);
            bb = lens[i].rcut / sqrt(1. + lens[i].emass);
        }
        else if ( lens[i].ct != 0 )
        {
            aa = lens[i].ct / sqrt(1. - lens[i].emass);
            bb = lens[i].ct / sqrt(1. + lens[i].emass);
        }
        else if (lens[i].type == 811)
        {
            // plot Re (ref Hjorth & Kneib Eq 13)
            aa = (0.75 * lens[i].rcut + 0.125 * lens[i].rc ); //sqrt(1.-lens[i].emass);
            bb = (0.75 * lens[i].rcut + 0.125 * lens[i].rc ); //sqrt(1.+lens[i].emass);
        }
        else if  (i > G.nmsgrid && i < G.nlens)
        {
            aa = lens[i].rc / sqrt(1. - lens[i].emass);
            bb = lens[i].rc / sqrt(1. + lens[i].emass);
        }
        else
        {
            aa = lens[i].ct / sqrt(1. - lens[i].emass);
            bb = lens[i].ct / sqrt(1. + lens[i].emass);
        };

        /* aa,bb, ct and rcut are in arcsec */
        // create 1 ellipse and 2 lines for the semi-axis
        // the ellipse is for rcut or Re if PIEMD
        fprintf(OUT, "%s %.7lf %.7lf %.4lf %.4lf %.2lf %.2lf\n",
                lens[i].n, Cx, Cy, aa, bb, RTD*lens[i].theta, 0.);
        fprintf(OUT, "%s %.6lf %.6lf %.4lf %.4lf %.2lf %.2lf\n",
                lens[i].n, Cx, Cy, aa, 0., lens[i].theta*RTD, 0.);
        fprintf(OUT, "%s %.7lf %.7lf %.4lf %.4lf %.2lf %.2lf\n",
                lens[i].n, Cx, Cy, 0., bb, lens[i].theta*RTD, 0.);

        // create an ellipse for b0/sigma
        if (lens[i].type > 0)
        {
            if (lens[i].type != 8 && lens[i].type<81 && lens[i].type>89)
            {
                bb = lens[i].rc / sqrt(1. - lens[i].emass);
                aa = lens[i].rc / sqrt(1. + lens[i].emass);
            }
            else
            {
                aa = sqrt(lens[i].b0) / sqrt(1. - lens[i].emass);
                bb = sqrt(lens[i].b0) / sqrt(1. + lens[i].emass);
            };

            //lens[i].masse =
            // plot an ellipse for b0/sigma
            fprintf(OUT, "%s %.7lf %.7lf %.4lf %.4lf %.2lf %.2lf\n",
                    lens[i].n, Cx, Cy, aa, bb, RTD*lens[i].theta, lens[i].masse);
        };
    };

    fclose(OUT);
}
