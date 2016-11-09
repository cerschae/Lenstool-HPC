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
/*                 Program  : r_observ              */
/*                 Version  : 16 septembre 1992         */
/*                 Location : Obs. Toulouse         */
/*                 Auteur   : Henri             */
/****************************************************************/

void r_observ(FILE *IN, FILE *OUT)
{
    extern struct g_mode     M;
    extern struct g_observ   O;
    char   second[20], third[FILENAME_SIZE+10];
    double disp;

    disp = 0.;
    O.gain = 0.;

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(FILENAME_SIZE+10)

        if (!strcmp(second, "bruit"))
        {
            sscanf(third, "%d", &O.bruit);
            fprintf(OUT, "\t%s\t\t%d\n", second, O.bruit);
        }
        else if ((!strcmp(second, "sky")) || (!strcmp(second, "SKY")))
        {
            sscanf(third, "%lf", &O.SKY);
            fprintf(OUT, "\t%s\t\t%lf\n", second, O.SKY);
        }
        else if (!strcmp(second, "idum"))
        {
            sscanf(third, "%d", &O.idum);
            fprintf(OUT, "\t%s\t\t%d\n", second, O.idum);
        }
        else if (!strcmp(second, "dispersion"))
        {
            sscanf(third, "%lf", &disp);
            fprintf(OUT, "\t%s\t\t%lf\n", second, disp);
        }
        else if (!strcmp(second, "gain"))
        {
            sscanf(third, "%lf", &O.gain);
            fprintf(OUT, "\t%s\t\t%lf\n", second, O.gain);
        }
        else if (!strcmp(second, "prec"))
        {
            sscanf(third, "%lf", &O.prec);
            fprintf(OUT, "\t%s\t\t%lf\n", second, O.prec);
        }
        else if (!strcmp(second, "binning"))
        {
            sscanf(third, "%d%d", &O.setbin, &O.bin);
            fprintf(OUT, "\t%s\t\t%d %d\n", second, O.setbin, O.bin);
        }
        else if (!strcmp(second, "seeing"))
        {
            sscanf(third, "%d%lf", &O.setseeing, &O.seeing);
            if (O.setseeing) {O.setseeing=1;}
            O.r0st = O.seeing / sqrt(log(2.));
            O.r0st = O.r0st * O.r0st;

            fprintf(OUT, "\t%s\t\t%d %lf\n", second, O.setseeing, O.seeing);
        }
        else if (!strcmp(second, "seeing_e"))
        {
            sscanf(third, "%d%lf%lf%lf", &O.setseeing, &O.seeing_a, &O.seeing_b, &O.seeing_angle);
            if (O.setseeing) {O.setseeing=2;}
            O.r0st_a = O.seeing_a / sqrt(log(2.));
            O.r0st_b = O.seeing_b / sqrt(log(2.));
            O.r0st_a = O.r0st_a * O.r0st_a;
            O.r0st_b = O.r0st_b * O.r0st_b;

            fprintf(OUT, "\t%s\t\t%d %lf %lf %lf\n", second, O.setseeing, O.seeing_a, O.seeing_b, O.seeing_angle);
        }
        else if (!strcmp(second, "psf"))
        {
            sscanf(third, "%d%s", &O.setseeing, O.psffile);
            if (O.setseeing) {O.setseeing=3;}
        }

        // Read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);

    if ((O.bruit != 0) && (O.gain == 0.))
    {
        if (disp != 0)
            O.gain = O.SKY / disp / disp;
        else
        {
            NPRINTF(stderr, "WARNING: ain or dispersion not defined, no noise will be added\n");
            O.bruit = 0;
        }
    }
}
