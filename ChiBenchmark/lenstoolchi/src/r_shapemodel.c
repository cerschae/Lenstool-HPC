#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*                 Program  : lenstool				*/
/*                 Version  : May, 14 2007			*/
/*                 Location : LAM, Marseille			*/
/*                 Auteur   : Benjamin Clément				*
 ****************************************************************
 * Read and set the source shape parameters. */

void r_shapemodel(FILE *IN,FILE *OUT, long int i)
{
    extern struct galaxie source[NFMAX];

	char	second[20],third[50];
	
	struct galaxie *sshape = &source[i];

	sshape->C.x = sshape->C.y = 0.;
	sshape->E.a = sshape->E.b = sshape->E.theta = 0.;
	sshape->mag = 0;
    sshape->z = 0;
    sshape->dl0s = sshape->dos = sshape->dr = -1;
    sshape->I0 = 50;
    sshape->c = 'g';
    sshape->type = 3;  // gaussian profile (default)
    sshape->var1 = 4;  // Sersic index (default: De Vaucouleur)
	sprintf(sshape->n, "S%ld", i);  // source name (default: S%d)

	fmot(IN,second);
	while(strcmp(second,"end"))
	{
		flire(IN,third); 

        if (!strcmp(second,"s_center_x") ||
            !strcmp(second, "x_center"))
		{
			sscanf(third,"%lf",&sshape->C.x);
			fprintf(OUT,"\t%s\t%lf\n",second,sshape->C.x);
		}
        else if ( !strcmp(second,"s_center_y") ||
                  !strcmp(second, "y_center"))
		{
			sscanf(third,"%lf",&sshape->C.y);
			fprintf(OUT,"\t%s\t%lf\n",second,sshape->C.y);
		}
		else if (!strcmp(second,"s_sigx") || 
                 !strcmp(second,"a_arcsec") )
	    {
			sscanf(third,"%lf",&sshape->E.a);
			fprintf(OUT,"\t%s\t%lf\n",second,sshape->E.a);
		} 
		else if (!strcmp(second,"s_sigy") ||
                 !strcmp(second,"b_arcsec") )
	    {
			sscanf(third,"%lf",&sshape->E.b);
			fprintf(OUT,"\t%s\t%lf\n",second,sshape->E.b);
		}
		else if (!strcmp(second,"s_eps") )
	    {
			sscanf(third,"%lf",&sshape->eps);
			fprintf(OUT,"\t%s\t%lf\n",second,sshape->eps);
		}
        else if (!strcmp(second,"s_angle") ||
                 !strcmp(second, "angle_pos"))
		{
			sscanf(third,"%lf",&sshape->E.theta);
			sshape->E.theta*=DTR;
			fprintf(OUT,"\t%s\t%lf\n",second,sshape->E.theta);
		}
        else if (!strcmp(second,"s_mag") ||
                 !strcmp(second, "mag"))
	    {
			sscanf(third,"%lf",&sshape->mag);
			fprintf(OUT,"\t%s\t%lf\n",second,sshape->mag);
		}
        else if (!strcmp(second,"index"))
	    {
			sscanf(third,"%lf",&sshape->var1);
			fprintf(OUT,"\t%s\t%lf\n",second,sshape->var1);
		}
        else if (!strcmp(second,"type"))
	    {
			sscanf(third,"%d",&sshape->type);
			fprintf(OUT,"\t%s\t%d\n",second,sshape->type);
		}
        else if (!strcmp(second,"z"))
	    {
			sscanf(third,"%lf",&sshape->z);
			fprintf(OUT,"\t%s\t%lf\n",second,sshape->z);
		}
        else if (!strcmp(second,"id") )
	    {
			sscanf(third,"%s",sshape->n);
			fprintf(OUT,"\t%s\t%s\n",second,sshape->n);
		}

		// Read the next line
		fmot(IN,second);
	}
			 
	fprintf(OUT,"\t%s\n",second);

}
