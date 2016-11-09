#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*                 Program  : grille                */
/*                 Version  : 1 mai 1992            */
/*                 Location : Obs. Toulouse         */
/*                 Auteur   : jean-paul             */
/****************************************************************/

void r_cosmologie(FILE *IN, FILE *OUT)
{
    extern struct g_cosmo C;

    char     second[20], third[FILENAME_SIZE+10];
    double ok;

    ok = -1;

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(FILENAME_SIZE+10)
        if ( !strcmp(second, "model") )                            //TV
        {                                                          //TV
            sscanf(third, "%d", &C.model);                        //TV
            fprintf(OUT, "\t%s\t%d\n", second, C.model);          //TV
        }
        else if ( !strcmp(second, "H0") )
        {
            sscanf(third, "%lf", &C.H0);
            fprintf(OUT, "\t%s\t%lf\n", second, C.H0);
            C.h = C.H0 / h0;
        }
        else if ( !strcmp(second, "omegaM") || !strcmp(second, "omega") )
        {
            sscanf(third, "%lf", &C.omegaM);
            fprintf(OUT, "\t%s\t%lf\n", second, C.omegaM);
        }
        else if ( !strcmp(second, "omegaX") || !strcmp(second, "lambda") )
        {
            sscanf(third, "%lf", &C.omegaX);
            fprintf(OUT, "\t%s\t%lf\n", second, C.omegaX);
        }
        else if ( !strcmp(second, "wX") || !strcmp(second, "q") || !strcmp(second, "w0") )     //TV  "q" for Model 2, "wX" for Model 3, "w0" for Model 4
        {
            sscanf(third, "%lf", &C.wX);
            fprintf(OUT, "\t%s\t%lf\n", second, C.wX);
        }
        else if ( !strcmp(second, "wa") || !strcmp(second, "n") || !strcmp(second, "delta") || !strcmp(second, "w1") )     //TV  "n" for Model 2, "delta" for model 3, "w1" for model 4  
        {
            sscanf(third, "%lf", &C.wa);
            fprintf(OUT, "\t%s\t%lf\n", second, C.wa);
        }
        else if ( !strcmp(second, "omegaK") )
        {
            sscanf(third, "%lf", &ok);
            fprintf(OUT, "\t%s\t%lf\n", second, ok);
        }

        // Read the next line
        fmot(IN, second);
    }

    // if a flat Universe
    if ( ok == 0. )
    {
        C.omegaX = 1 - C.omegaM;
        C.kcourb = 0.;
    }
    else
        C.kcourb = C.omegaM + C.omegaX - 1;

    fprintf(OUT, "\tkcourb\t%lf\n", C.kcourb);
    fprintf(OUT, "\t%s\n", second);

}


