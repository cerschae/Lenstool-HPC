#include<stdio.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*                 Program  : r_cline               */
/*                 Version  : 1 mai 1992            */
/*                 Location : Obs. Toulouse         */
/*                 Auteur   : jean-paul             */
/****************************************************************
 * Read the input parameters for the critical lines
 * See set_default.c for the default values
 */
void r_cline(FILE *IN, FILE *OUT)
{
    extern struct g_cline CL;
    char second[20], third[FILENAME_SIZE+10];
    char *pch;
    register int i;

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(FILENAME_SIZE+10)

        if (!strcmp(second, "nplan"))
        {
            pch = strtok(third, " ");
            sscanf(pch, "%d", &CL.nplan);
            if ( CL.nplan > NPZMAX )
                CL.nplan = NPZMAX;

            for ( i = 0 ; i < CL.nplan ; i++ )
            {
                pch = strtok(NULL, " ");
                sscanf(pch, "%lf", &CL.cz[i]);
            }

            fprintf(OUT, "\t%s\t\t%d", second, CL.nplan);

            for ( i = 0 ; i < CL.nplan ; i++)
                fprintf(OUT, " %lf ", CL.cz[i]);

            fprintf(OUT, "\n");
        }
        else if (!strcmp(second, "zonemult"))
        {
            sscanf(third, "%d%d%s", &CL.zone, &CL.npzone, CL.zonefile);
            fprintf(OUT, "\t%s\t\t%d %d %s\n", second, CL.zone,
                    CL.npzone, CL.zonefile);
        }
        else if (!strcmp(second, "dmax"))
        {
            sscanf(third, "%lf", &CL.dmax);
            fprintf(OUT, "\t%s\t\t%lf\n", second, CL.dmax);
        }
        else if (!strcmp(second, "pas") || !strcmp(second, "limitLow") || !strcmp(second, "limitlow") )
        {
            sscanf(third, "%lf", &CL.cpas);
            fprintf(OUT, "\t%s\t\t%lf\n", second, CL.cpas);
        }
        else if (!strcmp(second, "limitHigh") )
        {
            sscanf(third, "%lf", &CL.limitHigh);
            fprintf(OUT, "\t%s\t\t%lf\n", second, CL.limitHigh);
        }
        else if (!strncmp(second, "algo", 4))
        {
            sscanf(third, "%s", CL.algorithm);
            if ( strstr(upcase(CL.algorithm), "SNAKE") )
                strcpy(CL.algorithm, "SNAKE");
            fprintf(OUT, "\t%s\t\t%s\n", second, CL.algorithm);
        }

        // Read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);
}
