#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*                 Program  : grille                */
/*                 Version  : 1 mai 1992            */
/*                 Location : Obs. Toulouse         */
/*                 Auteur   : jean-paul             */
/****************************************************************/

void r_potentiel(FILE *IN, FILE *OUT, int i)
{
    extern struct g_mode   M;
    extern struct pot      lens[];

    struct pot  *ilens;
    char    second[100], third[255];

    ilens = &lens[i];

    ilens->C.x = ilens->C.y = 0.;
    ilens->emass = ilens->epot = 0.;
    ilens->alpha = ilens->beta = 0;
    ilens->theta = ilens->phi = 0.;
    ilens->mag = 0;
    ilens->rcut = ilens->rcutkpc = DBL_MAX;
    ilens->rc = ilens->rckpc = 0;
    ilens->masse = ilens->pmass = 0;
    ilens->z = 0;

    fprintf(OUT, "%d\n", i);

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(255)

        if ( !strcmp(second, "profil") ||
             !strcmp(second, "profile") )
        {
            sscanf(third, "%d", &ilens->type);
            fprintf(OUT, "\t%s\t\t%d\n", second, ilens->type);
        }
        else if (!strcmp(second, "nid"))
        {
            sscanf(third, "%s", ilens->n);
            fprintf(OUT, "\t%s\t\t%s\n", second, ilens->n);
        }
        else if (!strcmp(second, "x_centre"))
        {
            sscanf(third, "%lf", &ilens->C.x);
            fprintf(OUT, "\t%s\t%lf\n", second, ilens->C.x);
        }
        else if (!strcmp(second, "y_centre"))
        {
            sscanf(third, "%lf", &ilens->C.y);
            fprintf(OUT, "\t%s\t%lf\n", second, ilens->C.y);
        }
        else  if (!strcmp(second, "x_centre_wcs"))
        {
            sscanf(third, "%lf", &ilens->C.x);
            fprintf(OUT, "\t%s\t%lf\n", second, ilens->C.x);
            ilens->C.x -= M.ref_ra;
            ilens->C.x *= -3600 * cos(M.ref_dec * DTR);
        }
        else if (!strcmp(second, "y_centre_wcs"))
        {
            sscanf(third, "%lf", &ilens->C.y);
            fprintf(OUT, "\t%s\t%lf\n", second, ilens->C.y);
            ilens->C.y -= M.ref_dec;
            ilens->C.y *= 3600;
        }
        else if (!strcmp(second, "pmass"))
        {
            sscanf(third, "%lf", &ilens->pmass);
            fprintf(OUT, "\t%s\t%lf\n", second, ilens->pmass);
        }
        else if (!strcmp(second, "ellip_pot"))
        {
            sscanf(third, "%lf", &ilens->epot);
            fprintf(OUT, "\t%s\t%lf\n", second, ilens->epot);
        }
        else if ( !strcmp(second, "ellipticite") ||
                  !strcmp(second, "ellipticity") ||
                  !strcmp(second, "gamma") )
        {
            sscanf(third, "%lf", &ilens->emass);
            fprintf(OUT, "\t%s\t%lf\n", second, ilens->emass);
        }
        else if (!strcmp(second, "angle_pos"))
        {
            sscanf(third, "%lf", &ilens->theta);
            fprintf(OUT, "\t%s\t%lf\n", second, ilens->theta);
            ilens->theta *= DTR;
        }
		else if (!strcmp(second,"phi") )
		{
			sscanf(third,"%lf",&ilens->phi); 
			fprintf(OUT,"\t%s\t%lf\n",second,ilens->phi);
			ilens->phi*=DTR;
		}
        else if ( !strcmp(second, "core_radius") ||
                  !strcmp(second, "scale_radius") ||
                  !strcmp(second, "re") )
        {
            sscanf(third, "%lf", &ilens->rc);
            fprintf(OUT, "\t%s\t%lf\n", second, ilens->rc);
        }
        else if ( !strcmp(second, "core_radius_kpc") ||
                  !strcmp(second, "scale_radius_kpc") ||
                  !strcmp(second, "re_kpc") )
        {
            sscanf(third, "%lf", &ilens->rckpc);
            fprintf(OUT, "\t%s\t%lf\n", second, ilens->rckpc);
        }
        else if (!strcmp(second, "cut_radius") ||
                 !strcmp(second, "virial_radius") ||
                 !strcmp(second, "r200") )
        {
            sscanf(third, "%lf", &ilens->rcut);
            fprintf(OUT, "\t%s\t%lf\n", second, ilens->rcut);
        }
        else if (!strcmp(second, "cut_radius_kpc") ||
                 !strcmp(second, "virial_radius_kpc") ||
                 !strcmp(second, "r200_kpc") )
        {
            sscanf(third, "%lf", &ilens->rcutkpc);
            fprintf(OUT, "\t%s\t%lf\n", second, ilens->rcutkpc);
        }
        else if ( !strcmp(second, "v_disp") ||
                  !strcmp(second, "sigma_e") )
        {
            sscanf(third, "%lf", &ilens->sigma);
            fprintf(OUT, "\t%s\t\t%lf\n", second, ilens->sigma);
        }
        else if ( !strcmp(second, "exponent") ||
                  !strcmp(second, "alpha") ||
                  !strcmp(second, "exposant") ||
                  !strcmp(second, "n") )
        {
            sscanf(third, "%lf", &ilens->alpha);
            fprintf(OUT, "\t%s\t%lf\n", second, ilens->alpha);
        }
        else if (!strcmp(second, "beta") ||
                 !strcmp(second, "concentration") ||
                 !strcmp(second, "c") )
        {
            sscanf(third, "%lf", &ilens->beta);
            fprintf(OUT, "\t%s\t\t%lf\n", second, ilens->beta);
        }
        else if (!strcmp(second, "rc_slope"))
        {
            sscanf(third, "%lf", &ilens->rcslope);
            fprintf(OUT, "\t%s\t%lf\n", second, ilens->rcslope);
        }
        else if (!strcmp(second, "z_lens"))
        {
            sscanf(third, "%lf", &ilens->z);
            fprintf(OUT, "\t%s\t\t%lf\n", second, ilens->z);
        }
        else if (!strncmp(second, "mag", 3))
        {
            sscanf(third, "%lf", &ilens->mag);
            fprintf(OUT, "\t%s\t\t%lf\n", second, ilens->mag);
        }
        else if ( !strncmp(second, "virial_mass", 6) ||
                  !strcmp(second, "masse") ||
                  !strcmp(second, "m200") ||
                  !strcmp(second, "mass") )
        {
            sscanf(third, "%lf", &ilens->masse);
            fprintf(OUT, "\t%s\t\t%le\n", second, ilens->masse);
        }
        else if ( !strcmp(second, "rhos") )
        {
            sscanf(third, "%lf", &ilens->pmass);
            fprintf(OUT, "\t%s\t\t%le\n", second, ilens->pmass);
        }

        // Read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);
    if ( ilens->z == 0. )
    {
        fprintf(stderr, "ERROR: No redshift defined for potential %d\n", i);
        exit(-1);
    }
}
