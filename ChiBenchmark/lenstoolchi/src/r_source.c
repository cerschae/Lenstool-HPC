#include<stdio.h>
#include<string.h>
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

void r_source(FILE *IN, FILE *OUT)
{
    extern  struct  g_source    S;

    char    second[20], third[FILENAME_SIZE+10];

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(FILENAME_SIZE+10)

        if (!strcmp(second, "grid"))
        {
            sscanf(third, "%d", &S.grid);
            fprintf(OUT, "\t%s\t\t%d\n", second, S.grid);
        }
        else if (!strcmp(second, "random"))
        {
            sscanf(third, "%d", &S.rand);
            fprintf(OUT, "\t%s\t\t%d\n", second, S.rand);
        }
        else if (!strcmp(second, "n_source"))
        {
            sscanf(third, "%ld", &S.ns);
            fprintf(OUT, "\t%s\t%ld\n", second, S.ns);
        }
        else if (!strcmp(second, "elip_max"))
        {
            sscanf(third, "%lf", &S.emax);
            fprintf(OUT, "\t%s\t%lf\n", second, S.emax);
        }
        else if (!strcmp(second, "lf_alpha"))
        {
            sscanf(third, "%lf", &S.lfalpha);
            fprintf(OUT, "\t%s\t\t%.3lf\n", second, S.lfalpha);
        }
        else if (!strcmp(second, "lf_m_star"))
        {
            sscanf(third, "%lf", &S.lfm_star);
            fprintf(OUT, "\t%s\t\t%.3lf\n", second, S.lfm_star);
        }
        else if (!strcmp(second, "lf_m_min"))
        {
            sscanf(third, "%lf", &S.lfm_min);
            fprintf(OUT, "\t%s\t\t%.3lf\n", second, S.lfm_min);
        }
        else if (!strcmp(second, "lf_m_max"))
        {
            sscanf(third, "%lf", &S.lfm_max);
            fprintf(OUT, "\t%s\t\t%.3lf\n", second, S.lfm_max);
        }
        else if (!strcmp(second, "dist_z"))
        {
            sscanf(third, "%d%lf%lf", &S.distz, &S.zsmin, &S.zsmax);
            fprintf(OUT, "\t%s\t\t%d %lf %lf\n", second, S.distz, S.zsmin, S.zsmax);
        }
        else if (!strcmp(second, "smail"))
        {
            sscanf(third, "%lf%lf%lf", &S.par1, &S.par2, &S.par3);
            fprintf(OUT, "\t%s\t\t%lf %lf %lf\n", second, S.par1, S.par2, S.par3);
        }
        else if (!strcmp(second, "z_source"))
        {
            sscanf(third, "%lf", &S.zs);
            fprintf(OUT, "\t%s\t%lf\n", second, S.zs);
        }
        else if (!strcmp(second, "taille") || !strcmp(second, "size") )
        {
            sscanf(third, "%lf", &S.taille);
            fprintf(OUT, "\t%s\t\t%lf\n", second, S.taille);
        }

        // Read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);
}
