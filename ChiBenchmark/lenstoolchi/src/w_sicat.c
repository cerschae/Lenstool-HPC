#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        w_sicat             */
/*      auteur:     Jean-Paul Kneib         */
/****************************************************************/

void    w_sicat(struct galaxie ima[NFMAX][NIMAX], struct galaxie ss[NFMAX])
{
    const extern struct g_source  S;
    const extern struct g_mode    M;
    long  int i;
    int j;
    FILE    *OUT;

    OUT = fopen("sipos.dat", "w");

    fprintf(OUT, "#REFERENCE 3 %.7lf %.7lf\n", M.ref_ra, M.ref_dec);

    for ( i = 0 ; i < S.ns ; i++ )
    {
        if ( ss[i].E.theta > PI )
            ss[i].E.theta -= PI;

        for (j = 0; j < NIMAX && strcmp(ima[i][j].n, "") ; j++ )
        {

            fprintf(OUT, "%s %.3lf %.3lf %.3lf %.3lf %.2lf %.3lf %.2lf %.3lf %.3lf %.2lf %.2lf\n",
                    ss[i].n, ss[i].C.x, ss[i].C.y, ss[i].E.a, ss[i].E.b,
                    RTD*ss[i].E.theta, ss[i].z, ss[i].mag,
                    ima[i][j].C.x, ima[i][j].C.y, ss[i].mag - ima[i][j].mag,
                    ima[i][j].mag);
        };
    };

    fclose(OUT);
}
