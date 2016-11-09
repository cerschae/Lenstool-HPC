#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*                 Program  : grille                */
/*                 Version  : 1 mai 1992            */
/*                 Location : Obs. Toulouse         */
/*                 Auteur   : jean-paul             */
/****************************************************************/

void r_grille(FILE *IN, FILE *OUT)
{
    extern struct g_mode    M;
    extern struct g_grille  G;

    char    second[20], third[FILENAME_SIZE+10];

    G.dx = G.dy = 5.;
    G.nlens_crit = -1;

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(FILENAME_SIZE+10)

        if (!strcmp(second, "nombre") || !strcmp(second, "number"))
        {
            sscanf(third, "%d", &G.ngrid);
            if (G.ngrid < 10)
                G.ngrid = 10;
            else if (G.ngrid > NGGMAX)
            {
                NPRINTF( stderr, "WARN: number=%d too large. Maximum value=%d\n", G.ngrid, NGGMAX);
                G.ngrid = NGGMAX;
            }
            G.ngrid = 2 * ((int)(G.ngrid / 2));
            fprintf(OUT, "\t%s\t\t%d\n", second, G.ngrid);
        }
        else if (!strcmp(second, "polaire") || !strcmp(second, "polar"))
        {
            sscanf(third, "%d", &G.pol);
            fprintf(OUT, "\t%s\t\t%d\n", second, G.pol);
        }
        else if (!strcmp(second, "nlentille") || !strcmp(second, "nlens"))
        {
            sscanf(third, "%ld", &G.nlens);
            fprintf(OUT, "\t%s\t%ld\n", second, G.nlens);
        }
        else if (!strcmp(second, "nlens_opt"))
        {
            sscanf(third, "%ld", &G.no_lens);
            fprintf(OUT, "\t%s\t%ld\n", second, G.no_lens);
        }
        else if (!strcmp(second, "nlens_crit"))
        {
            sscanf(third, "%ld", &G.nlens_crit);
            fprintf(OUT, "\t%s\t%ld\n", second, G.nlens_crit);
        }
        else if (!strcmp(second, "sp_exc"))
        {
            sscanf(third, "%lf%lf", &G.exc, &G.excmin);
            fprintf(OUT, "\t%s\t%lf %lf\n", second, G.exc, G.excmin);
        }
        else if (!strcmp(second, "echant"))
        {
            sscanf(third, "%d", &G.echant);
            fprintf(OUT, "\t%s\t%d\n", second, G.echant);
        }
        else if (!strcmp(second, "spline"))
        {
            sscanf(third, "%s", G.splinefile);
            fprintf(OUT, "\t%s\t%s\n", second, G.splinefile);
        }
        else if ( !strcmp(second, "nmsgrid"))
        {
            sscanf(third, "%ld", &G.nmsgrid);
            fprintf(OUT, "\t%s\t%ld\n", second, G.nmsgrid);
        }

        // Read the next line
        fmot(IN, second);
    }


    fprintf(OUT, "\t%s\n", second);

    if ((G.nlens_crit > G.nlens) || (G.nlens_crit == -1))
        G.nlens_crit = G.nlens;

    if (G.no_lens > G.nlens)
        G.no_lens = G.nlens;

}
