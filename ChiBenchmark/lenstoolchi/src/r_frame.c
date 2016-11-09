#include<stdio.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*                 Program  : pictp             */
/*                 Version  : 1 mai 1992            */
/*                 Location : Obs. Toulouse         */
/*                 Auteur   : jean-paul             */
/****************************************************************/

void r_frame(FILE *IN, FILE *OUT)
{
    char    second[20], third[FILENAME_SIZE+10];
    double  reel;
    extern struct   g_frame F;

    F.rmax = 0.;
    F.lmin = 0.;
    F.lmax = 0.;

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(FILENAME_SIZE+10)

        if (!strcmp(second, "dmax"))
        {
            sscanf(third, "%lf", &reel);
            fprintf(OUT, "\t%s\t%lf\n", second, reel);
            F.xmin = -reel;
            F.xmax = reel;
            F.ymin = -reel;
            F.ymax = reel;
            F.rmax = reel * sqrt(2.);
        }
        else if (!strcmp(second, "xmin"))
        {
            sscanf(third, "%lf", &F.xmin);
            fprintf(OUT, "\t%s\t%lf\n", second, F.xmin);
        }
        else if (!strcmp(second, "xmax"))
        {
            sscanf(third, "%lf", &F.xmax);
            fprintf(OUT, "\t%s\t%lf\n", second, F.xmax);
        }
        else if (!strcmp(second, "ymin"))
        {
            sscanf(third, "%lf", &F.ymin);
            fprintf(OUT, "\t%s\t%lf\n", second, F.ymin);
        }
        else if (!strcmp(second, "ymax"))
        {
            sscanf(third, "%lf", &F.ymax);
            fprintf(OUT, "\t%s\t%lf\n", second, F.ymax);
        }
        else if (!strcmp(second, "lmin"))
        {
            sscanf(third, "%lf", &F.lmin);
            fprintf(OUT, "\t%s\t%lf\n", second, F.lmin);
        }
        else if (!strcmp(second, "lmax"))
        {
            sscanf(third, "%lf", &F.lmax);
            fprintf(OUT, "\t%s\t%lf\n", second, F.lmax);
        }

        // Read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);

    if (F.rmax == 0.)
        F.rmax = Max(Max(fabs(F.xmax), fabs(F.xmin)), Max(fabs(F.ymax), fabs(F.ymin))) * sqrt(2.);

    if (F.rmax == 0.)
        F.rmax = 30.;

}
