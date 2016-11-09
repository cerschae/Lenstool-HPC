#include<stdio.h>
#include<string.h>
#include "fonction.h"
#include "structure.h"
#include<lt.h>

/****************************************************************/
/*                 Program  : lenstool              */
/*                 Version  : 22 february 2008          */
/*                 Location : OAMP          */
/*                 Auteur   : eric              */
/****************************************************************/

void r_msgrid(FILE *IN, FILE *OUT)
{
    extern struct g_msgrid  H;

    char    second[20], third[FILENAME_SIZE+10];

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(FILENAME_SIZE+10)

        if ( !strcmp(second, "threshold"))
        {
            sscanf(third, "%lf", &H.threshold);
            fprintf(OUT, "\t%s\t%lf\n", second, H.threshold);
        }
        else if ( !strcmp(second, "levels"))
        {
            sscanf(third, "%d", &H.levels);
            fprintf(OUT, "\t%s\t%d\n", second, H.levels);
        }
        else if ( !strcmp(second, "param"))
        {
            sscanf(third, "%lf", &H.param);
            fprintf(OUT, "\t%s\t%lf\n", second, H.param);
        }
        else if ( !strcmp(second, "gridfile"))
        {
            sscanf(third, "%s", H.gridfile);
            fprintf(OUT, "\t%s\t%s\n", second, H.gridfile);
        }

        // Read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);
}
