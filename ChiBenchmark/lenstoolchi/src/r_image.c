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
/*                 Program  : r_image               */
/*                 Version  : 1 mai 1992            */
/*                 Location : Obs. Toulouse         */
/*                 Auteur   : jean-paul             */
/****************************************************************/

static void scanzmlimit(char *third, int *opt, char *name, int *bk, double *min, double *max, double *prec);

void r_image(FILE *IN, FILE *OUT)
{
    extern struct g_mode   M;
    extern struct g_image  I;
    extern struct MCarlo   mc;
    extern struct z_lim    zlim[];
    extern struct z_lim    zalim;
    extern struct cline    cl[];
    extern struct sigposStr sigposAs;
    char   second[20], third[FILENAME_SIZE+10];
    int    ii, j;
    double x;

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(FILENAME_SIZE+10)

        if (!strcmp(second, "arcletstat"))
        {
            sscanf(third, "%d%d%s", &I.stat, &I.statmode, I.arclet);

            fprintf(OUT, "\t%s\t%d %d %s\n",
                    second, I.stat, I.statmode, I.arclet);
        }
        else if ( !strcmp(second, "sigell"))
        {
            sscanf(third, "%lf %lf", &I.sigell, &I.dsigell);
            fprintf(OUT, "\t%s\t%lf ", second, I.sigell);
            if ( I.dsigell != -1. )
                fprintf(OUT, "%lf", I.dsigell);
            fprintf(OUT, "\n");
        }
        else if (!strcmp(second, "shearmap"))
        {
            sscanf(third, "%d%lf%s", &I.shmap, &I.zsh, I.shfile);

            fprintf(OUT, "\t%s\t%d %lf %s\n", second, I.shmap, I.zsh, I.shfile);
        }
        else if (!strcmp(second, "z_arclet"))
        {
            sscanf(third, "%lf", &I.zarclet);
            fprintf(OUT, "\t%s\t%lf\n", second,I.zarclet);
        }
        else if (!strcmp(second, "multfile"))
        {
            sscanf(third, "%d%s", &I.n_mult, I.multfile);
            fprintf(OUT, "\t%s\t%d %s\n", second, I.n_mult, I.multfile);
        }
        else if (!strcmp(second, "sourcefit"))
        {
            sscanf(third, "%d%s%s", &I.srcfit, I.srcfitFile, I.srcfitMethod);
            fprintf(OUT, "\t%s\t%d %s %s\n", second, I.srcfit, I.srcfitFile, I.srcfitMethod);
        }
        else if (!strcmp(second, "mult_wcs"))
        {
            sscanf(third, "%d", &I.mult_abs);
            fprintf(OUT, "\t%s\t%d\n", second, I.mult_abs);
        }
        else if (!strcmp(second, "sigpos"))
        {
            sscanf(third, "%lf", &x);
            sigposAs.min = sqrt(x);
            fprintf(OUT, "\t%s\t\t%lf\n", second, sigposAs.min);
        }
        else if (!strcmp(second, "sigposArcsec"))
        {
            j = getWords(third);
            if ( j == 3 )
            {
                if ( sscanf(third, "%d %lf %lf", &sigposAs.bk, &sigposAs.min, &sigposAs.max) == 3 )
                    fprintf(OUT, "\t%s\t\t%d %lf %lf\n", second, sigposAs.bk, sigposAs.min, sigposAs.max);
                else
                    fprintf(OUT, "\t%s\tERROR reading  '%s'!!\n", second, third);
            }
            else if ( j == 1 )
            {
                if ( sscanf(third, "%lf", &sigposAs.min) == 1 )
                    fprintf(OUT, "\t%s\t\t%lf\n", second, sigposAs.min);
                else
                    fprintf(OUT, "\t%s\tERROR reading  '%s'!!\n", second, third);
            }
            else
            {
                fprintf( stderr, "ERROR: Syntax mismatched for sigposArcsec keyword\n");
                exit(1);
            }
        }
        else if (!strcmp(second, "sigamp"))
        {
            sscanf(third, "%lf", &x);
            fprintf(OUT, "\t%s\t\t%lf\n", second, x);
            I.sig2amp = x * x;
        }
        else if (!strcmp(second, "Dmag"))
        {
            sscanf(third, "%lf", &I.Dmag);

            fprintf(OUT, "\t%s\t\t%lf\n", second, I.Dmag);
        }
        else if (!strcmp(second, "forme"))
        {
            sscanf(third, "%d", &I.forme);

            fprintf(OUT, "\t%s\t\t%d\n", second, I.forme);
        }
        else if (!strcmp(second, "n_MonteCarlo"))
        {
            sscanf(third, "%d", &mc.n_MonteCarlo);

            fprintf(OUT, "\t%s\t\t%d\n", second, mc.n_MonteCarlo);
        }
        else if (!strcmp(second, "optMC"))
        {
            sscanf(third, "%d", &mc.optMC);

            fprintf(OUT, "\t%s\t\t%d\n", second, mc.optMC);
        }
        else if (!strcmp(second, "iterations"))
        {
            sscanf(third, "%d", &mc.iterations);

            fprintf(OUT, "\t%s\t\t%d\n", second, mc.iterations);
        }
        else if (!strcmp(second, "tosses_sq"))
        {
            sscanf(third, "%d", &mc.tosses_sq);

            fprintf(OUT, "\t%s\t\t%d\n", second, mc.tosses_sq);
        }
        else if (!strcmp(second, "squares_par"))
        {
            sscanf(third, "%d", &mc.squares_par);

            fprintf(OUT, "\t%s\t\t%d\n", second, mc.squares_par);
        }
        else if (!strcmp(second, "z_m_limit"))
        {
            j = I.nzlim;
            scanzmlimit(third, &zlim[j].opt, zlim[j].n, &zlim[j].bk,
                        &zlim[j].min, &zlim[j].max, &zlim[j].dderr);

            fprintf(OUT, "\t%s\t %d %s %d %.3lf  %.3lf  %.4lf \n", second, zlim[j].opt,
                    zlim[j].n, zlim[j].bk, zlim[j].min, zlim[j].max, zlim[j].dderr);
            if (zlim[j].opt > 0)
                I.nzlim++;
        }
        else if (!strcmp(second, "z_opt"))
        {
            j = I.nzlim;
            sscanf(third, "%lf%d", &zlim[j-1].percent, &zlim[j-1].bk0);

            fprintf(OUT, "\t%s\t\t%.3lf%d\n", second, zlim[j-1].percent, zlim[j-1].bk0);
        }
        else if (!strcmp(second, "z_a_limit"))
        {
            sscanf(third, "%d%lf%lf", &zalim.bk, &zalim.min, &zalim.max);
            fprintf(OUT, "\t%s\t %d %.4lf %.4lf\n", second, zalim.bk, zalim.min, zalim.max);
        }
        else if (!strcmp(second, "critic"))
        {
            j = I.npcl;
            if( sscanf(third, "%d%lf%lf%lf%lf%lf", &ii, &cl[j].C.x, &cl[j].C.y,
                   &cl[j].phi, &cl[j].dl, &cl[j].z) != 6 )
            {
                fprintf(stderr, "ERROR: reading .par file\n>\t%s %s\n", second, third);
                exit(1);
            }

            fprintf(OUT, "\t%s\t\t%d %lf %lf  %lf  %lf  %lf\n", second, ii,
                    cl[j].C.x, cl[j].C.y, cl[j].phi, cl[j].dl, cl[j].z);
            cl[j].phi *= DTR;
            if (ii > 0)
            {
                cl[j].n = ii;
                I.npcl++;
            }
        }
        else if (!strcmp(second, "adjust"))
        {
            sscanf(third, "%d%s", &I.adjust, I.Afile);
            fprintf(OUT, "\t%s\t\t%d %s\n", second, I.adjust, I.Afile);
        }
        else if (!strcmp(second, "nfilt"))
        {
            sscanf(third, "%d", &I.Anfilt);
            fprintf(OUT, "\t%s\t\t%d\n", second, I.Anfilt);
        }
        else if (!strcmp(second, "npixseeing"))
        {
            sscanf(third, "%d", &I.Anpixseeing);
            fprintf(OUT, "\t%s\t\t%d\n", second, I.Anpixseeing);
        }
        else if (!strcmp(second, "seeing"))
        {
            sscanf(third, "%lf", &I.Aseeing);
            fprintf(OUT, "\t%s\t\t%lf\n", second, I.Aseeing);
        }

        // Read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);
}

static void scanzmlimit(char *third, int *opt, char *name, int *bk, double *min, double *max, double *prec)
{
    int i, n; // number of names
    char str[50], tmp[50];
    char *pch;

    n = 1;  // First assume only individual system
    // if n>1, n systems have the same redshift

    // Initialize values
    *opt = *bk = -1;
    *min = *max = *prec = -1.;

    strcpy( tmp, third );
    pch = strtok( tmp, " " );
    if ( pch != NULL )
        sscanf( pch, "%d", opt);
    else
        goto ERROR;

    pch = strtok( NULL, " " );
    if ( pch != NULL )
        sscanf( pch, "%s", name );
    else
        goto ERROR;

    // Read additional system names
    pch = strtok( NULL, " " );
    while( sscanf(pch, "%d", bk) != 1 && n < ZMBOUND)
    {
        sscanf(pch, "%s", str);
        name[strlen(name)] = ' ';
        strcat(name, str);
        n++;
        pch = strtok( NULL, " " );
    }

    if( n == ZMBOUND )
        goto ERROR;

    pch = strtok( NULL, " " );
    if ( pch != NULL )
        sscanf( pch, "%lf", min);
    else
        goto ERROR;

    pch = strtok( NULL, " " );
    if ( pch != NULL )
        sscanf( pch, "%lf", max);
    else
        goto ERROR;

    pch = strtok( NULL, " " );
    if ( pch != NULL )
        sscanf( pch, "%lf", prec);
    else
        goto ERROR;

    pch = strtok( NULL, " " );
    if ( pch == NULL || pch[0] == '#' )
        return;

ERROR:
    if ( n >= ZMBOUND )
    {
        fprintf(stderr, "ERROR: Too many bounded images %s (max %d)\n",
                third, ZMBOUND);
        exit(-1);
    }
    else if ( *prec == -1 )
    {
        fprintf(stderr, "ERROR: Missing argument in z_m_limit %s\n", third);
        exit(-1);
    }
    else 
    {
        fprintf(stderr, "ERROR: Reading line %s. Check the number of arguments\n", third);
        exit(-1);
    }
}
