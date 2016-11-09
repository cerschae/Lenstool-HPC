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

void r_large(FILE *IN, FILE *OUT)
{
    char    second[20], third[FILENAME_SIZE+10];
    extern struct g_large L;


    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(FILENAME_SIZE+10)

        if (!strcmp(second, "iso"))
        {
            sscanf(third, "%d%d%lf%lf%lf",
                   &L.iso, &L.nmaxiso, &L.scale, &L.zonex, &L.zoney);

            fprintf(OUT, "\t%s\t%d %d %lf %lf %lf\n", second,
                    L.iso, L.nmaxiso, L.scale, L.zonex, L.zoney);
        }
        else if (!strcmp(second, "name"))
        {
            sscanf(third, "%s", L.iname);
            fprintf(OUT, "\t%s\t%s\n", second, L.iname);
        }
        else if (!strcmp(second, "profil"))
        {
            sscanf(third, "%d%d", &L.profil, &L.pt);
            fprintf(OUT, "\t%s\t%d %d\n", second, L.profil, L.pt);
        }
        else if (!strcmp(second, "contour"))
        {
            sscanf(third, "%d%d", &L.ncourbe, &L.npt);
            fprintf(OUT, "\t%s\t%d %d\n", second, L.ncourbe, L.npt);
        }
        else if (!strcmp(second, "vitesse"))
        {
            sscanf(third, "%d", &L.vitesse);
            fprintf(OUT, "\t%s\t%d\n", second, L.vitesse);
        }
        else if (!strcmp(second, "large_dist"))
        {
            sscanf(third, "%lf", &L.dlarge);
            fprintf(OUT, "\t%s\t%lf\n", second, L.dlarge);
        }

        // Read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);

}
