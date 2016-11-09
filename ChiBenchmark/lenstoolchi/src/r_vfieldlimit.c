#include<stdio.h>
#include<string.h>
#include<fonction.h>
#include<dimension.h>
#include<constant.h>
#include<structure.h>

/****************************************************************/
/*                 Program  : lenstool              */
/*                 Date     : 15/05/2013         */
/*                 Location : Obs.Lyon             */
/*                 Auteur   : Johan Richard             *
 ****************************************************************
 * Read and set the limits on the velocity field parameters.
 */

void r_vfieldlimit(FILE *IN, FILE *OUT)
{
    extern  struct  vfield     vfmin, vfmax;
    extern  int             vfblock[NPAMAX];

    char    second[20], third[FILENAME_SIZE];

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(FILENAME_SIZE+10)

        if (!strcmp(second,"x_centre"))
        {
            sscanf(third, "%d %lf %lf", &vfblock[VFCX],
                   &vfmin.C.x, &vfmax.C.x);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, vfblock[VFCX],
                    vfmin.C.x, vfmax.C.x);
        }
        else if ( !strcmp(second,"y_centre"))
        {
            sscanf(third, "%d %lf %lf", &vfblock[VFCY],
                   &vfmin.C.y, &vfmax.C.y);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, vfblock[VFCY],
                    vfmin.C.y, vfmax.C.y);
        }
	else if (!strcmp(second,"vt"))
        {
            sscanf(third, "%d %lf %lf", &vfblock[VFVT],
                   &vfmin.vt, &vfmax.vt);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, vfblock[VFVT],
                    vfmin.vt, vfmax.vt);
        }
	else if (!strcmp(second,"rt"))
        {
            sscanf(third, "%d %lf %lf", &vfblock[VFRT],
                   &vfmin.rt, &vfmax.rt);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, vfblock[VFRT],
                    vfmin.rt, vfmax.rt);
        }
	else if (!strcmp(second,"i"))
        {
            sscanf(third, "%d %lf %lf", &vfblock[VFI],
                   &vfmin.i, &vfmax.i);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, vfblock[VFI],
                    vfmin.i, vfmax.i);
            vfmin.i *= DTR;
            vfmax.i *= DTR;
        }
        else if (!strcmp(second,"theta"))
        {
            sscanf(third, "%d %lf %lf", &vfblock[VFTHETA],
                   &vfmin.theta, &vfmax.theta);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, vfblock[VFTHETA],
                    vfmin.theta, vfmax.theta);
            vfmin.theta *= DTR;
            vfmax.theta *= DTR;
        }
        else if (!strcmp(second,"lcent"))
        {
            sscanf(third, "%d %lf %lf", &vfblock[VFLCENT],
                   &vfmin.lcent, &vfmax.lcent);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, vfblock[VFLCENT],
                    vfmin.lcent, vfmax.lcent);
        }
        else if (!strcmp(second,"sigma"))
        {
            sscanf(third, "%d %lf %lf", &vfblock[VFSIGMA],
                   &vfmin.sigma, &vfmax.sigma);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, vfblock[VFSIGMA],
                    vfmin.sigma, vfmax.sigma);
        }

        // Read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);

}


