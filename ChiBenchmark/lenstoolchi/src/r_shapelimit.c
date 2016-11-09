#include<stdio.h>
#include<string.h>
#include<fonction.h>
#include<dimension.h>
#include<constant.h>
#include<structure.h>

/****************************************************************/
/*                 Program  : lenstool              */
/*                 Date     : 03/10/2011         */
/*                 Location : Marseille            */
/*                 Auteur   : Eric Jullo                *
 ****************************************************************
 * Read and set the limits on the source parameters.
 */

void r_shapelimit(FILE *IN, FILE *OUT, long int i)
{
    extern  struct  galaxie     smin[NFMAX], smax[NFMAX];
    extern  int             sblock[NFMAX][NPAMAX];

    char    second[20], third[FILENAME_SIZE];

    fprintf(OUT, "%ld\n", i);

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(FILENAME_SIZE+10)

        if (!strcmp(second,"s_center_x") ||
            !strcmp(second, "x_center"))
        {
            sscanf(third, "%d %lf %lf", &sblock[i][SCX],
                   &smin[i].C.x, &smax[i].C.x);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, sblock[i][SCX],
                    smin[i].C.x, smax[i].C.x);
        }
        else if ( !strcmp(second,"s_center_y") ||
                  !strcmp(second, "y_center"))
        {
            sscanf(third, "%d %lf %lf", &sblock[i][SCY],
                   &smin[i].C.y, &smax[i].C.y);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, sblock[i][SCY],
                    smin[i].C.y, smax[i].C.y);
        }
		else if (!strcmp(second,"s_sigx") || 
                 !strcmp(second,"a_arcsec") )
        {
            sscanf(third, "%d %lf %lf", &sblock[i][SA],
                   &smin[i].E.a, &smax[i].E.a);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, sblock[i][SA],
                    smin[i].E.a, smax[i].E.a);
        }
		else if (!strcmp(second,"s_sigy") ||
                 !strcmp(second,"b_arcsec") )
        {
            sscanf(third, "%d %lf %lf", &sblock[i][SB],
                   &smin[i].E.b, &smax[i].E.b);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, sblock[i][SB],
                    smin[i].E.b, smax[i].E.b);
        }
		else if (!strcmp(second,"s_eps") )
        {
            sscanf(third, "%d %lf %lf", &sblock[i][SEPS],
                   &smin[i].eps, &smax[i].eps);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, sblock[i][SEPS],
                    smin[i].eps, smax[i].eps);
        }
        else if (!strcmp(second,"s_angle") ||
                 !strcmp(second, "angle_pos"))
        {
            sscanf(third, "%d %lf %lf", &sblock[i][STHETA],
                   &smin[i].E.theta, &smax[i].E.theta);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, sblock[i][STHETA],
                    smin[i].E.theta, smax[i].E.theta);
            smin[i].E.theta *= DTR;
            smax[i].E.theta *= DTR;
        }
        else if (!strcmp(second,"s_mag") ||
                 !strcmp(second, "mag"))
        {
            sscanf(third, "%d %lf %lf", &sblock[i][SFLUX],
                   &smin[i].mag, &smax[i].mag);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, sblock[i][SFLUX],
                    smin[i].mag, smax[i].mag);
        }
        else if (!strcmp(second,"index")) 
        {
            sscanf(third, "%d %lf %lf", &sblock[i][SINDEX],
                   &smin[i].var1, &smax[i].var1);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, sblock[i][SINDEX],
                    smin[i].var1, smax[i].var1);
        }

        // Read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);

}


