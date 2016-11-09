#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        ecrire              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/*                              */
/*   Modified :                         */
/*      EJ (01/09/05)                   */
/****************************************************************
 * Write a list of galaxie structures in a file.
 *
 * flag definition :
 * 0 : Do not print the shear info. The coordinates are absolute.
 * 1 : Print the shear info. The coordinates are absolute.
 * 2 : Do not print the shear info. The coordinates are relative.
 * 3 : Print the shear info. The coordinates are relative.
 * 4
 */
void ecrire_r(long int nstart, long int nstop, struct galaxie *liste, char *name, int flag)
{
    const extern struct g_mode    M;
    FILE    *OUT;
    long int i;
    double   Cx, Cy;

    OUT = fopen(name, "w");
    if ( OUT == NULL )
    {
        fprintf(stderr, "ERROR: Unable to open %s in writting mode\n",name);
        exit(1);
    }

    if ( flag & 2 )
        fprintf( OUT, "#REFERENCE 3 %.7f %.7f\n", M.ref_ra, M.ref_dec );
    else
        fprintf( OUT, "#REFERENCE 0 %.7f %.7f\n", M.ref_ra, M.ref_dec );

    for (i = nstart; i < nstop; i++)
    {
        if (liste[i].E.theta > PI)
            liste[i].E.theta -= PI;

        if (liste[i].c == 's')
        {
            liste[i].E.a = liste[i].E.b = 0.;
        };

        Cx = liste[i].C.x;
        Cy = liste[i].C.y;

        // Convert the coordinates to absolute WCS coordinates if possible
        if ( M.iref != 0 && !(flag & 2) )
        {
            Cy = Cy / 3600 + M.ref_dec;
            Cx = -Cx / 3600 / cos(M.ref_dec * DTR) + M.ref_ra;
        }

        fprintf( OUT, "%s %.7lf %.7lf %.6lf %.6lf %.5lf %.4lf %.2lf ",
                 liste[i].n, Cx, Cy, liste[i].E.a, liste[i].E.b,
                 liste[i].E.theta*RTD, liste[i].z, liste[i].mag );

        if ( flag & 1 )
            fprintf(OUT, "%.3lf %.3lf %.3lf %.3lf ",
                    liste[i].kappa, liste[i].gamma1, liste[i].gamma2,
                    sqrt(liste[i].gamma1*liste[i].gamma1 + liste[i].gamma2*liste[i].gamma2) );

        if ( flag & 4 )
            fprintf(OUT, "%.2lf %.2lf ", liste[i].var1, liste[i].var2 );

        fprintf(OUT, "\n");

    };

    fclose(OUT);
}
