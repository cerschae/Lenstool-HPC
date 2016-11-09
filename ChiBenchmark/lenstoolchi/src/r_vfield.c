#include<stdio.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*                 Program  : r_vfield             */
/*                 Version  : 15 may 2013            */
/*                 Location : Obs. Lyon             */
/*                 Auteur   : Johan Richard         */
/****************************************************************/

void r_vfield(FILE *IN, FILE *OUT)
{
    char    second[20], third[FILENAME_SIZE+10];
    double  reel;
    extern struct  vfield vf;

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(FILENAME_SIZE+10)

	double vt;
	double rt;
	struct point C;
	double i;
	double theta;
	double lcent;
	double sigma;
	int profile;
	double dlamb;
	double nb_slice;
	double lmin;

        if (!strcmp(second, "x_centre"))
        {
            sscanf(third, "%lf",&vf.C.x);
            fprintf(OUT, "\t%s\t%lf\n",second,vf.C.x);
        }
        else if (!strcmp(second, "y_centre"))
        {
            sscanf(third, "%lf",&vf.C.y);
            fprintf(OUT, "\t%s\t%lf\n",second,vf.C.y);
        }
        else if (!strcmp(second, "profile"))
        {
            sscanf(third, "%d",&vf.profile);
            fprintf(OUT, "\t%s\t%d\n",second,vf.profile);
        }
        else if (!strcmp(second, "vt"))
        {
            sscanf(third, "%lf",&vf.vt);
            fprintf(OUT, "\t%s\t%lf\n",second,vf.vt);
        }
        else if (!strcmp(second, "rt"))
        {
            sscanf(third, "%lf",&vf.rt);
            fprintf(OUT, "\t%s\t%lf\n",second,vf.rt);
        }
        else if (!strcmp(second, "i"))
        {
            sscanf(third, "%lf",&vf.i);
            fprintf(OUT, "\t%s\t%lf\n",second,vf.i);
            vf.i *= DTR;
        }
        else if (!strcmp(second, "theta"))
        {
            sscanf(third, "%lf",&vf.theta);
            fprintf(OUT, "\t%s\t%lf\n",second,vf.theta);
            vf.theta *= DTR;
        }
        else if (!strcmp(second, "lcent"))
        {
            sscanf(third, "%lf",&vf.lcent);
            fprintf(OUT, "\t%s\t%lf\n",second,vf.lcent);
        }
        else if (!strcmp(second, "sigma"))
        {
            sscanf(third, "%lf",&vf.sigma);
            fprintf(OUT, "\t%s\t%lf\n",second,vf.sigma);
        }

        // Read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);
}
