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
/*                 Program  : grille                */
/*                 Version  : 1 mai 1992            */
/*                 Location : Obs. Toulouse         */
/*                 Auteur   : jean-paul             */
/****************************************************************/

void r_limit(FILE *IN, FILE *OUT, int i)
{
    char                    second[100], third[255];
    extern int              block[NLMAX][NPAMAX];
    extern struct   pot     lmin[], lmax[], prec[];

    fprintf(OUT, "%d\n", i);

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(255)

        if (!strcmp(second, "x_centre"))
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][CX],
                   &lmin[i].C.x, &lmax[i].C.x, &prec[i].C.x);
            fprintf(OUT, "\t%s\t%d %lf %lf %lf\n", second, block[i][CX],
                    lmin[i].C.x, lmax[i].C.x, prec[i].C.x);
        }
        else if (!strcmp(second, "y_centre"))
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][CY],
                   &lmin[i].C.y, &lmax[i].C.y, &prec[i].C.y);
            fprintf(OUT, "\t%s\t%d %lf %lf %lf\n", second, block[i][CY],
                    lmin[i].C.y, lmax[i].C.y, prec[i].C.y);
        }
        else if (!strcmp(second, "ellip_pot"))
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][EPOT],
                   &lmin[i].epot, &lmax[i].epot, &prec[i].epot);
            fprintf(OUT, "\t%s\t%d %lf %lf %lf\n", second, block[i][EPOT],
                    lmin[i].epot, lmax[i].epot, prec[i].epot);
        }
        else if ( !strcmp(second, "ellipticite") ||
                  !strcmp(second, "ellipticity") ||
                  !strcmp(second, "gamma") )
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][EMASS],
                   &lmin[i].emass, &lmax[i].emass, &prec[i].emass);
            fprintf(OUT, "\t%s\t%d %lf %lf %lf\n", second, block[i][EMASS],
                    lmin[i].emass, lmax[i].emass, prec[i].emass);
        }
        else if (!strcmp(second, "angle_pos"))
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][THETA],
                   &lmin[i].theta, &lmax[i].theta, &prec[i].theta);
            lmin[i].theta *= DTR;
            lmax[i].theta *= DTR;
            prec[i].theta *= DTR;
            fprintf(OUT, "\t%s\t%d %lf %lf %lf\n", second, block[i][THETA],
                    lmin[i].theta, lmax[i].theta, prec[i].theta);
        }
		else if (!strcmp(second,"phi"))
		{
			sscanf(third,"%d%lf%lf%lf",&block[i][PHI],
				&lmin[i].phi,&lmax[i].phi,&prec[i].phi); 
			
			lmin[i].phi*=DTR;
			lmax[i].phi*=DTR;
			prec[i].phi*=DTR;
			fprintf(OUT,"\t%s\t%d %lf %lf %lf\n",second,block[i][PHI],
				lmin[i].phi,lmax[i].phi,prec[i].phi);
		}
        else if ( !strcmp(second, "cut_radius") ||
                  !strcmp(second, "virial_radius") ||
                  !strcmp(second, "r200") )
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][RCUT],
                   &lmin[i].rcut, &lmax[i].rcut, &prec[i].rcut);
            fprintf(OUT, "\t%s\t%d %lf %lf %lf\n", second, block[i][RCUT],
                    lmin[i].rcut, lmax[i].rcut, prec[i].rcut);
        }
        else if ( !strcmp(second, "cut_radius_kpc") ||
                  !strcmp(second, "virial_radius_kpc") ||
                  !strcmp(second, "r200_kpc") )
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][RCUT],
                   &lmin[i].rcutkpc, &lmax[i].rcutkpc, &prec[i].rcutkpc);
            fprintf(OUT, "\t%s\t%d %lf %lf %lf\n", second, block[i][RCUT],
                    lmin[i].rcutkpc, lmax[i].rcutkpc, prec[i].rcutkpc);
        }
        else if ( !strcmp(second, "core_radius") ||
                  !strcmp(second, "scale_radius") ||
                  !strcmp(second, "eff_radius") )
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][RC],
                   &lmin[i].rc, &lmax[i].rc, &prec[i].rc);
            fprintf(OUT, "\t%s\t%d %lf %lf %lf\n", second, block[i][RC],
                    lmin[i].rc, lmax[i].rc, prec[i].rc);
        }
        else if ( !strcmp(second, "core_radius_kpc") ||
                  !strcmp(second, "scale_radius_kpc") ||
                  !strcmp(second, "eff_radius_kpc" ) )
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][RC],
                   &lmin[i].rckpc, &lmax[i].rckpc, &prec[i].rckpc);
            fprintf(OUT, "\t%s\t%d %lf %lf %lf\n", second, block[i][RC],
                    lmin[i].rckpc, lmax[i].rckpc, prec[i].rckpc);
        }
        else if ( !strcmp(second, "v_disp") ||
                  !strcmp(second, "sigma_e") )
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][B0],
                   &lmin[i].sigma, &lmax[i].sigma, &prec[i].sigma);
            fprintf(OUT, "\t%s\t\t%d %lf %lf %lf\n", second, block[i][B0],
                    lmin[i].sigma, lmax[i].sigma, prec[i].sigma);
        }
        else if (!strcmp(second, "pmass"))
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][PMASS],
                   &lmin[i].pmass, &lmax[i].pmass, &prec[i].pmass);
            fprintf(OUT, "\t%s\t\t%d %lf %lf %lf\n", second, block[i][PMASS],
                    lmin[i].pmass, lmax[i].pmass, prec[i].pmass);
        }
        else if ( !strcmp(second, "exponent") ||
                  !strcmp(second, "alpha") ||
                  !strcmp(second, "n") )
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][ALPHA],
                   &lmin[i].alpha, &lmax[i].alpha, &prec[i].alpha);
            fprintf(OUT, "\t%s\t%d %lf %lf %lf\n", second, block[i][ALPHA],
                    lmin[i].alpha, lmax[i].alpha, prec[i].alpha);
        }
        else if ( !strcmp(second, "beta") ||
                  !strcmp(second, "concentration") ||
                  !strcmp(second, "c") )
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][BETA],
                   &lmin[i].beta, &lmax[i].beta, &prec[i].beta);
            fprintf(OUT, "\t%s\t\t%d %lf %lf %lf\n", second, block[i][BETA],
                    lmin[i].beta, lmax[i].beta, prec[i].beta);
        }
        else if (!strcmp(second, "virial_mass") ||
                 !strcmp(second, "masse") ||
                 !strcmp(second, "m200") ||
                 !strcmp(second, "mass"))
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][MASSE],
                   &lmin[i].masse, &lmax[i].masse, &prec[i].masse);
            fprintf(OUT, "\t%s\t\t%d %le %le %le\n", second, block[i][MASSE],
                    lmin[i].masse, lmax[i].masse, prec[i].masse);
        }
        else if ( !strcmp(second, "rhos") )
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][PMASS],
                   &lmin[i].pmass, &lmax[i].pmass, &prec[i].pmass);
            fprintf(OUT, "\t%s\t%d %le %le %le\n", second, block[i][PMASS],
                    lmin[i].pmass, lmax[i].pmass, prec[i].pmass);
        }
        else if ( !strcmp(second, "z_lens") )
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][ZLENS],
                   &lmin[i].z, &lmax[i].z, &prec[i].z);
            fprintf(OUT, "\t%s\t%d %le %le %le\n", second, block[i][ZLENS],
                    lmin[i].z, lmax[i].z, prec[i].z);
        }
        else if ( !strcmp(second, "rc_slope") )
        {
            sscanf(third, "%d%lf%lf%lf", &block[i][RCSLOPE],
                   &lmin[i].rcslope, &lmax[i].rcslope, &prec[i].rcslope);
            fprintf(OUT, "\t%s\t%d %le %le %le\n", second, block[i][RCSLOPE],
                    lmin[i].rcslope, lmax[i].rcslope, prec[i].rcslope);
        }

        // Read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);
}
