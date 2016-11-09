#include<stdio.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*                 Program  : lenstool              */
/*                 Version  : Feb 13, 2006          */
/*                 Location : ESO, Chile            */
/*                 Auteur   : Eric Jullo                *
 ****************************************************************
 * Read and set the limits on the cosmological parameters.
 */

void r_cosmolimit(FILE *IN, FILE *OUT)
{
    extern  struct  g_cosmo     clmin, clmax;
    extern  int             cblock[NPAMAX];

    char    second[20], third[FILENAME_SIZE+10];

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(FILENAME_SIZE+10)

        if (!strcmp(second, "omegaM") || !strcmp(second, "omega"))
        {
            sscanf(third, "%d %lf %lf", &cblock[OMEGAM],
                   &clmin.omegaM, &clmax.omegaM);

            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, cblock[OMEGAM],
                    clmin.omegaM, clmax.omegaM);
        }
        else if (!strcmp(second, "omegaX") ||
                 !strcmp(second, "lambda") ||
                 !strcmp(second, "omegaL")    )
        {
            sscanf(third, "%d %lf %lf", &cblock[OMEGAX],
                   &clmin.omegaX, &clmax.omegaX);

            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, cblock[OMEGAX],
                    clmin.omegaX, clmax.omegaX);
        }
        else if (!strcmp(second, "wX"))
        {
            sscanf(third, "%d %lf %lf", &cblock[WX],
                   &clmin.wX, &clmax.wX);

            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, cblock[WX],
                    clmin.wX, clmax.wX);
        }
        else if (!strcmp(second, "wa"))
        {
            sscanf(third, "%d %lf %lf", &cblock[WA],
                   &clmin.wa, &clmax.wa);

            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, cblock[WA],
                    clmin.wa, clmax.wa);
        };

        // Read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);

}


