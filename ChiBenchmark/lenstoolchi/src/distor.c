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
/*      nom:        distor              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

/* Write the <dist.dat> file.
 */
static long int ampli_error(struct galaxie *image, long int ni);

void    distor(struct galaxie *image, long int ni)
{
    const extern struct g_mode    M;
    const extern struct pot       lens[];

    long int     i;
    double  t0, r, R, E;
    FILE    *OUT;
    struct ellipse amp;
    char   id0[IDSIZE];

    int     errorsOK;
    char   parity1[NAMAX], parity2[NAMAX];
    struct point pref;  // point of REFERENCE in runmode

    for (i = 0; i < ni; i++)
    {
        r = image[i].E.a / image[i].E.b;
        image[i].q = 1. / r;
        R = r * r;
        E = image[i].eps = (R - 1.) / (R + 1.);
        image[i].tau = (R - 1.) / 2. / r;
        image[i].time = e_time(image[i].C, image[i].dr);
        amp = e_unmag_gal(&image[i]);
        image[i].A = 1 / fabs(amp.a * amp.b);
        parity1[i] = amp.a > 0 ? '+' : '-';
        parity2[i] = amp.b > 0 ? '+' : '-';
    }

    NPRINTF(stderr, "COMP: dist.dat file ");
    errorsOK = ampli_error(image, ni);
    NPRINTF(stderr, "\n");

    if (M.sort == 1)
        sort(ni, image, comparer_tau);
    else if (M.sort == 2)
        sort(ni, image, comparer_pos);

    OUT = fopen("dist.dat", "w");

    fprintf(OUT, "#REFERENCE 3 %.7f %.7f\n", M.ref_ra, M.ref_dec);
    fprintf(OUT, "   #ID       X         Y         R       EPS       TAU       AMP       ");
    if ( errorsOK )  fprintf(OUT, "E_AMP      ");
    fprintf(OUT, "DMAG  TIME[days]  DTIME    PARITY\n");

    t0 = image[0].time; strcpy(id0, image[0].n);
    pref.x = pref.y = 0.;
    for (i = 0; i < ni; i++)
    {
        if ( indexCmp(image[i].n, id0) )
        {
            t0 = image[i].time; strcpy(id0, image[i].n);
        }

        fprintf(OUT, "%6s  %8.3lf  %8.3lf  %8.3lf  %8.3lf  %8.3lf  %8.3lf  ",
                image[i].n, image[i].C.x, image[i].C.y, dist(image[i].C, pref),
                image[i].eps, image[i].tau, image[i].A);
        if ( errorsOK ) fprintf(OUT, "%8.3lf  ", image[i].mu);
        fprintf(OUT, "%8.2lf  %8.1lf  %8.1lf    %1c%1c\n",
                -2.5*log10(1 / image[i].A),
                image[i].time, image[i].time - t0, parity1[i], parity2[i]);
    };

    fclose(OUT);
}

/* Compute the mean and stddev value of the amplification for
 * a list of images.
 * mean value in image[i].A and stddev value in image[i].mu
 *
 * Return >1 if the error have been computed, 0 otherwise.
 */
static long int ampli_error(struct galaxie *image, long int ni)
{
    const extern struct   g_mode M;

    struct  ellipse amp;
    double **array;  //contains the bayes.dat file
    double  tmp;
    int     nParam;
    long int     i;
    long int nVal, iVal;

    array = readBayesModels(&nParam, &nVal);

    if ( nVal == 0 )
        return 0;

    // Set all amplifications and errors to 0
    for ( i = 0; i < ni; i++)
	{
        image[i].A = image[i].mu = 0;
		image[i].grad2.a = image[i].grad2.c = 0.;  // reset not optimized potentials
	}

    for ( iVal = 0; iVal < nVal; iVal++ )
    {
        NPRINTF(stderr, "\rCOMP: dist.dat file (Amplification errors %ld/%ld)",
                iVal + 1, nVal);
        setBayesModel(iVal, nVal, array);
        for ( i = 0; i < ni; i++)
        {
            amp = e_unmag_gal(&image[i]);
            tmp = 1. / fabs(amp.a * amp.b);
            image[i].A += tmp;
            image[i].mu += tmp * tmp;
        }
    }

    // finalize the mean and stddev computation
    for ( i = 0; i < ni; i++)
    {
        tmp = image[i].A /= nVal;       // mean
        image[i].mu = (image[i].mu - nVal * tmp * tmp) / (nVal - 1);  // variance
        image[i].mu = sqrt(image[i].mu);    // stddev
    }

    NPRINTF(stderr, "\n");
    free(array);
    return nVal;
}
