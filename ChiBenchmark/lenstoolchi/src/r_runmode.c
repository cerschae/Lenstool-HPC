#include<stdio.h>
#include<stdlib.h>
#include<string.h>
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

static void scanInverse(char * third);

void r_runmode(FILE *IN, FILE *OUT)
{
    extern struct   g_mode  M;

    char   second[20], third[FILENAME_SIZE+10];
    char   ref1[20], ref2[20];
    double ss0, tt0;
    int    hh0, mm0, dd0, nn0;

    M.iref = 0;     // default : no WCS
    M.nshearf = 25; // default for shearfield
    
    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(FILENAME_SIZE+10)

        if ((!strcmp(second, "arclet")) || (!strcmp(second, "image")))
        {
            sscanf(third, "%d%s", &M.image, M.imafile);
            fprintf(OUT, "\t%s\t%d %s\n", second, M.image, M.imafile);

            if (strlen(M.imafile) < 2)
                strcpy(M.imafile, "image.ext");
        }
        else if (!strcmp(second, "source"))
        {
            sscanf(third, "%d%s", &M.source, M.sourfile);
            fprintf(OUT, "\t%s\t%d %s\n", second, M.source, M.sourfile);


            if (strlen(M.sourfile) < 2)
                strcpy(M.sourfile, "source.ext");
        }
        else if (!strcmp(second, "sourceof"))
        {
            sscanf(third, "%d%s%s", &M.sof, M.imsfile, M.sfile);
            fprintf(OUT, "\t%s\t%d %s %s\n", second, M.sof, M.imsfile, M.sfile);
        }
        else if (!strcmp(second, "corshear"))
        {
            sscanf(third, "%d%s", &M.icorshear, M.corshfile);
            fprintf(OUT, "\t%s\t%d %s\n", second, M.icorshear, M.corshfile);
        }
        else if (!strcmp(second, "local"))
        {
            sscanf(third, "%d%s", &M.local, M.localfile);
            fprintf(OUT, "\t%s\t%d %s\n", second, M.local, M.localfile);

            if (strlen(M.localfile) < 2)
                strcpy(M.localfile, "local.ext");
        }
        else if (!strcmp(second, "study"))
        {
            sscanf(third, "%d%s", &M.study, M.studyfile);
            fprintf(OUT, "\t%s\t%d %s\n", second, M.study, M.studyfile);

            if (strlen(M.studyfile) < 2)
                strcpy(M.studyfile, "study.ext");
        }
        else if (!strcmp(second, "fake"))
        {
            sscanf(third, "%d", &M.fake);
            fprintf(OUT, "\t%s\t%d\n", second, M.fake);
        }
        else if (!strcmp(second, "meanz"))
        {
            sscanf(third, "%d", &M.mean);
            fprintf(OUT, "\t%s\t%d\n", second, M.mean);
        }
        else if (!strcmp(second, "verbose"))
        {
            sscanf(third, "%d", &M.verbose);
            fprintf(OUT, "\t%s\t%d\n", second, M.verbose);
        }
        else if (!strcmp(second, "imseeing"))
        {
            sscanf(third, "%lf", &M.seeing);
            fprintf(OUT, "\t%s\t%lf\n", second, M.seeing);
        }
        else if (!strcmp(second, "grille"))
        {

            sscanf(third, "%d%d%lf", &M.grille, &M.ngrille, &M.zgrille);
            fprintf(OUT, "\t%s\t%d %d %lf\n", second,
                    M.grille, M.ngrille, M.zgrille);
        }
        else if (!strcmp(second, "sort"))
        {
            sscanf(third, "%d", &M.sort);
            fprintf(OUT, "\t%s\t%d\n", second, M.sort);
        }
        else if (!strcmp(second, "inverse"))
        {
            scanInverse(third);
            if ( M.inverse <= 2 )
                fprintf(OUT, "\t%s\t%d %d\n", second, M.inverse, M.itmax);
            else
                fprintf(OUT, "\t%s\t%d %lf %d\n", second, M.inverse, M.rate, M.itmax);
        }
        else if (!strcmp(second, "minchi0"))
        {
            sscanf(third, "%lf", &M.minchi0);
            fprintf(OUT, "\t%s\t\t%lf\n", second, M.minchi0);
        }
        else if (!strcmp(second, "shearfield"))
        {
            sscanf(third, "%d%lf%s%d", &M.ishearf, &M.zshearf, M.shearffile, &M.nshearf);
            fprintf(OUT, "\t%s\t%d %lf %s %d\n", second, M.ishearf,
                    M.zshearf, M.shearffile, M.nshearf);
        }
        else if (!strcmp(second, "amplifield"))
        {
            sscanf(third, "%d%lf%s", &M.iamplif, &M.zamplif, M.ampliffile);
            fprintf(OUT, "\t%s\t%d %lf %s\n", second, M.iamplif,
                    M.zamplif, M.ampliffile);
        }
        else if (!strcmp(second, "shear"))
        {
            sscanf(third, "%d%d%lf%s", &M.ishear, &M.nshear, &M.zshear, M.shearfile);
            fprintf(OUT, "\t%s\t%d %d %lf %s\n", second, M.ishear, M.nshear,
                    M.zshear, M.shearfile);
        }
        else if (!strcmp(second, "poten"))
        {
            sscanf(third, "%d%d%lf%s", &M.ipoten, &M.npoten, &M.zpoten, M.potenfile);
            fprintf(OUT, "\t%s\t%d %d %lf %s\n", second, M.ipoten, M.npoten,
                    M.zpoten, M.potenfile);
        }
        else if (!strcmp(second, "mass"))
        {
            if ( sscanf(third, "%d%d%lf%s", &M.imass, &M.nmass, &M.zmass, M.massfile) != 4 )
            {
                printf("ERROR: Syntax error in line:\n\t%s\t%s [%d] [%d] [%lf] [%s]\n", second, third,M.imass,M.nmass,M.zmass,M.massfile);
                exit(-1);
            }
            fprintf(OUT, "\t%s\t%d %d %lf %s\n", second, M.imass, M.nmass,
                    M.zmass, M.massfile);
        }
        else if (!strcmp(second, "dpl"))
        {
            sscanf(third, "%d%d%lf%s%s",
                   &M.idpl, &M.ndpl, &M.zdpl, M.dplxfile, M.dplyfile);
            fprintf(OUT, "\t%s\t%d %d %lf %s %s\n", second, M.idpl, M.ndpl,
                    M.zdpl, M.dplxfile, M.dplyfile);
        }
        else if (!strcmp(second, "curv"))
        {
            sscanf(third, "%d%d%lf%s%s%s",
                   &M.icurv, &M.ncurv, &M.zcurv, M.cxxfile, M.cxyfile, M.cyyfile);
            fprintf(OUT, "\t%s\t%d %d %lf %s %s %s\n", second, M.icurv, M.ncurv,
                    M.zcurv, M.cxxfile, M.cxyfile, M.cyyfile);
        }
        else if (!strcmp(second, "ampli"))
        {
            sscanf(third, "%d%d%lf%s", &M.iampli, &M.nampli, &M.zampli, M.amplifile);
            fprintf(OUT, "\t%s\t%d %d %lf %s\n", second, M.iampli, M.nampli,
                    M.zampli, M.amplifile);
        }
        else if (!strcmp(second, "time"))
        {
            sscanf(third, "%d%d%lf%s", &M.itime, &M.ntime, &M.ztime, M.timefile);
            fprintf(OUT, "\t%s\t%d %d %lf %s\n", second, M.itime, M.ntime,
                    M.ztime, M.timefile);
        }
        else if (!strcmp(second, "prop"))
        {
            sscanf(third, "%d%d%lf%s", &M.prop, &M.nprop, &M.zprop, M.propfile);
            fprintf(OUT, "\t%s\t%d %d %.3lf %s\n", second,
                    M.prop, M.nprop, M.zprop, M.propfile);
        }
        else if (!strcmp(second, "propradius"))
        {
            sscanf(third, "%lf", &M.radius);
            fprintf(OUT, "\t%s\t%lf\n", second, M.radius);
        }
        else if (!strcmp(second, "pixel"))
        {
            sscanf(third, "%d%d%s", &M.pixel, &M.npixel, M.pixelfile);
            fprintf(OUT, "\t%s\t%d %d %s\n", second, M.pixel, M.npixel,
                    M.pixelfile);
        }
	else if (!strcmp(second, "cube"))
	{
		sscanf(third, "%d%d%d%s", &M.cube, &M.npixel, &M.nslices, M.cubefile);
		fprintf(OUT, "\t%s\t%d \t%d %d %s\n", second, M.cube, M.npixel,
				M.nslices, M.cubefile);
	}
        else if (!strcmp(second, "marker"))
        {
            sscanf(third, "%d%lf%s", &M.marker, &M.zmarker, M.markfile);
            fprintf(OUT, "\t%s\t%d %lf %s\n", second, M.marker, M.zmarker,
                    M.markfile);
        }
        else if (!strcmp(second, "reference"))
        {
            sscanf(third, "%d%s%s", &M.iref, ref1, ref2);
            if (M.iref == 1)
            {
                sscanf(ref1, "%d:%d:%lf", &hh0, &mm0, &ss0);
                sscanf(ref2, "%d:%d:%lf", &dd0, &nn0, &tt0);

                NPRINTF(stderr, "%d:%d:%lf %d:%d:%lf\n", hh0, mm0, ss0, dd0, nn0, tt0);

                M.ref_ra = ((double)(hh0) + ((double)(mm0)) / 60 + ss0 / 3600) * 15;

                if (dd0 < 0)
                    M.ref_dec = ((double)(dd0)) - ((double)(nn0)) / 60 - tt0 / 3600;
                else
                    M.ref_dec = ((double)(dd0)) + ((double)(nn0)) / 60 + tt0 / 3600;
                M.iref = 3; // for compatibility in the best.par and bestopt.par
            }
            else if ( M.iref == 3 )
            {
                sscanf(ref1, "%lf", &M.ref_ra);
                sscanf(ref2, "%lf", &M.ref_dec);
            }

            fprintf(OUT, "\t%s\t%d %s %s %lf %lf\n", second, M.iref, ref1, ref2,
                    M.ref_ra, M.ref_dec);
        }
        else if (!strcmp(second, "radialprop") || !strcmp(second, "profil"))
        {
            sscanf(third, "%d%lf%lf", &M.radial, &M.zradial, &M.theta);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, M.radial, M.zradial, M.theta);
            M.theta *= DTR;
        }
        else if (!strcmp(second, "chi2"))
        {
            sscanf(third, "%d", &M.ichi2);
            fprintf(OUT, "\t%s\t%d\n", second, M.ichi2);
            M.theta *= DTR;
        }

        // read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);

}

static void scanInverse(char * third)
{
    extern struct g_mode   M;
    double inverse1; // rate in bayesian mode (3) or number of iterations (mode 1 or 2)
    int    inverse2; // optional number of iterations in bayesian inverse mode (3)
    char   *pch;

    pch = strtok( third, " ");
    sscanf( pch, "%d", &M.inverse );

    if ( M.inverse <= 2 )
    {
        M.itmax = 100;  // Default value
        pch = strtok( NULL, " ");
        if ( pch != NULL && sscanf(pch, "%lf", &inverse1) == 1 )
            M.itmax = (int) inverse1;
    }
    else 
    {
        M.itmax = 1000;  // Default number of posterior samples
        M.rate = 0.5;    // Default convergence speed/evidence precision(5)
        pch = strtok( NULL, " ");
        if ( pch != NULL &&  sscanf(pch, "%lf", &inverse1) == 1 )
        {
            M.rate = inverse1;
            pch = strtok( NULL, " ");
            if ( pch != NULL &&  sscanf(pch, "%d", &inverse2) == 1 )
                M.itmax = inverse2;

        }
    }
}
