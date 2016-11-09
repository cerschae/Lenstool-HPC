#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"

static void medstat(double *data, int ndata, double *med, double *q, double *tq);
static void st_ez(  int mode, struct galaxie imas[NASMAX], int nima,
                    double *z, double *drz, double ts[NASMAX][200], double theta[NASMAX][200],
                    double amp[NASMAX][200],
                    double pz[NASMAX][200], double ipz[NASMAX][200],
                    double kappa[NASMAX], double gam[NASMAX], double thetap[NASMAX],
                    double sumpz[NASMAX], double tauix[NASMAX],
                    double tauiy[NASMAX], int jmax );
static void sp_dx(double *xx, double *yy, double *yy2, int nn, double *dydx);


/****************************************************************/
/*      nom:        study           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    study_pg(int type, double seeing, char studyfile[], int fake)
{
    extern  struct  g_mode  M;
    extern  struct  pot lens[];
    //extern    double  z_dlsds;

    struct  galaxie study[NASMAX]; //source,
    //struct    ellipse ampli;

    register int    i, j, k;
    int l, nst = 0, jmax;
    double  z[200], drz[200];
    double  ts[NASMAX][200], theta[NASMAX][200];
    double  amp[NASMAX][200], pz[NASMAX][200], ipz[NASMAX][200];
    double  zpm[NASMAX], nmag[NASMAX];
    double  zpmf[NASMAX], nmagf[NASMAX];
    double  gam[NASMAX], kappa[NASMAX], thetap[NASMAX], tauix[NASMAX], tauiy[NASMAX];
    double  sumpz[NASMAX];

    double   mo, **zmag;

    double  magbin[15];
    int nz[15], nzc[15];
    double  meanz[15], dispz[15];
    int nzf[15];
    double  meanzf[15], dispzf[15];
    double   datmed[15][100];
    double   median[15], quart[15], tquart[15];

    int nbinm = 6;
    double  sizebinm = 1;
    double  zccor = 0.3;
    double  dzz, zzmin, zzmax;
    double  mzmin = .25, mzmax = 4.5;

    FILE    *OUT;

    NPRINTF(stderr, "START\n");

    jmax = 200;

    /*
    * ALLOC: areas
    */

    NPRINTF(stderr, "ALLOC: areas \n");

    zmag = (double **) alloc_square_double(15, jmax);

    /*
    * read arclet catalogue
    */

    if (type == 1)
        f_shape((long int*)&nst, study, studyfile,0);
    else if (type == 2)
        f_shape2((long int*)&nst, study, studyfile);
    else if (type == 3)
        f_shape3((long int*)&nst, study, studyfile);
    else if (type == 4)
        f_shape4((long int*)&nst, study, studyfile);
    else
    {
        NPRINTF(stderr, "ERROR: Unrecognize study type %d\n", type);
        return;
    }

    /*
    * do a seeing correction if necessary
    */

    if (seeing > 0.)
    {
        NPRINTF(stderr, "COMP: seeing correction %lf\n", seeing);
        cor_seeing(nst, study, seeing);
    };


    /* definition of redshift intervals for inversion */
    /* do a linear scale in redshift */

    dzz = 0.05;
    zzmin = lens[0].z + dzz;
    zzmax = 5.;

    NPRINTF(stderr, "COMP: D ratio \n");
    z[0] = zzmin;
    for (j = 0; (j < jmax) && (z[j] <= zzmax); j++)
    {
        z[j] = zzmin + dzz * j;
        drz[j] = dratio(lens[0].z, z[j]);
    }

    /*
    * Compute the ellipticity variation with redshift
    */

    NPRINTF(stderr, "COMP: compute variation of ellipticities with redshift\n");
    st_ez(0, study, nst, z, drz, ts, theta, amp, pz, ipz, kappa, gam, thetap, sumpz, tauix, tauiy, jmax);
    NPRINTF(stderr, "COMP: done ...\n");

    /*
    * Automatic search of the most probable redshift using the z-probability
    */

    NPRINTF(stderr, "COMP: determination of the most probable redshift\n");
    st_opt(0, study, nst, z, ts, pz, ipz, amp, kappa, gam, thetap, sumpz, tauix, tauiy, jmax, zpm, nmag);


    /* compute mean/median and dispersion of the redshift distribution
    versus magnitude */


    /* initialize variables*/

    for (i = 0; i < nbinm; i++)
    {
        magbin[i] = 21. + i * sizebinm;
        meanz[i] = dispz[i] = 0.;
        nz[i] = nzc[i] = 0;
    }

    /* get z and mag and put them into zmag array */

    for (k = 0; k < nst; k++)
    {
        if ((sumpz[k] > 0.) && (zpm[k] > mzmin) && (zpm[k] < mzmax))
        {
            for (j = 0; (j < jmax) && (z[j] <= mzmax); j++)
            {
                for (i = 0; (i < nbinm); i++)
                {
                    if (amp[k][j] > 0.)
                        mo = study[k].mag - 2.5 * log10(amp[k][j]);
                    else
                    {
                        NPRINTF(stderr, "WARNING: log10: Domain error in study.c");
                        mo = study[k].mag;
                    }
                    if ((mo >= magbin[i]) && (mo < magbin[i] + 1.))
                    {
                        zmag[i][j] += pz[k][j];
                    }
                }
            }
        }
    }


    OUT = fopen("zmag.dat", "w");
    fprintf(OUT, "#i n magbin meanz mean+1s mean-1s median quart tquart\n");

    /* get z and mag and compute the mean-median and the dispersion */

    for (i = 0; i < nbinm; i++)
    {
        l = 0;
        for (j = 0; j < 100; j++)
            datmed[i][j] = 0.;

        for (j = 0; j < nst; j++)
        {
            if ((nmag[j] > magbin[i] - .5) && (nmag[j] <= magbin[i+1] + .5) && (zpm[j] > mzmin) && (zpm[j] < mzmax))
            {
                printf("%d %lf %lf\n", i, nmag[j], zpm[j]);
                datmed[i][l++] = zpm[j];
                meanz[i] += zpm[j];
                dispz[i] += zpm[j] * zpm[j];
                nz[i]++;
                if (zpm[j] < zccor)
                    nzc[i]++;
            }
        }

        medstat(datmed[i], l, &median[i], &quart[i], &tquart[i]);

        if (nz[i] != 0)
        {
            meanz[i] /= nz[i];
            dispz[i] /= nz[i];
            dispz[i] -= meanz[i] * meanz[i];
            dispz[i] = sqrt(dispz[i]);
            printf("%d %d %lf %lf %lf\n", i, nz[i], meanz[i], dispz[i], median[i]);
            fprintf(OUT, "%d %d %.1lf %.4lf %.4lf %.4lf %.3lf %.3lf %.3lf\n",
                    i, nz[i], magbin[i] + 0.5*sizebinm, meanz[i], meanz[i] + dispz[i], meanz[i] - dispz[i],
                    median[i], quart[i], tquart[i]);
        }
    }
    fclose(OUT);


    /*
    * do the same for the fake one if necessary
    */

    if (fake == 1)
    {
        NPRINTF(stderr, "COMP: compute variation of ellipticities with redshift\n");
        st_ez(1, study, nst, z, drz, ts, theta, amp, pz, ipz, kappa, gam, thetap, sumpz, tauix, tauiy, jmax);
        NPRINTF(stderr, "COMP: determination of the most probable redshift\n");
        st_opt(1, study, nst, z, ts, pz, ipz, amp, kappa, gam, thetap, sumpz, tauix, tauiy, jmax, zpmf, nmagf);


        for (i = 0; i < nbinm; i++)
        {
            magbin[i] = 21. + i * sizebinm;
            meanzf[i] = dispzf[i] = 0.;
            nzf[i] = 0;
        }


        for (k = 0; k < nst; k++)
        {
            if ((sumpz[k] > 0.) && (zpm[k] > mzmin) && (zpm[k] < mzmax))
            {
                for (j = 0; (j < jmax) && (z[j] < mzmax); j++)
                {
                    for (i = 0; (i < nbinm); i++)
                    {
                        if (amp[k][j] > 0.)
                            mo = study[k].mag - 2.5 * log10(amp[k][j]);
                        else
                        {
                            NPRINTF(stderr, "WARNING: log10: Domain error in study.c");
                            mo = study[k].mag;
                        }
                        if ((mo >= magbin[i]) && (mo < magbin[i] + 1.))
                        {
                            zmag[i][j] -= pz[k][j];
                        }
                    }
                }
            }
        }
        wrf_fits("zmag.fits", zmag, 15, jmax, z[0], z[jmax-1], 22., 30.);


        OUT = fopen("zmagf.dat", "w");
        fprintf(OUT, "#i n magbin meanz mean+1s mean-1s\n");
        for (i = 0; i < nbinm; i++)
        {
            for (j = 0; j < nst; j++)
            {
                if ((nmagf[j] > magbin[i] - .5) && (nmagf[j] <= magbin[i+1] + .5) && (zpmf[j] > mzmin) && (zpmf[j] < zccor))
                {
                    meanzf[i] += zpmf[j];
                    dispzf[i] += zpmf[j] * zpmf[j];
                    nzf[i] += 1;
                }
            }
            if (nzf[i] != 0)
            {
                meanzf[i] /= nzf[i];
                dispzf[i] /= nzf[i];
                dispzf[i] -= meanzf[i] * meanzf[i];
                dispzf[i] = sqrt(dispzf[i]);

                fprintf(OUT, "%d %d %.1lf %.4lf %.4lf %.4lf\n", i, nzf[i], magbin[i] + 0.5*sizebinm, meanzf[i],
                        meanzf[i] + dispzf[i], meanzf[i] - dispzf[i]);
            }
        }
        fclose(OUT);


        OUT = fopen("zcor.dat", "w");
        fprintf(OUT, "#bin ncor ntrue-nfake mag zmeancor zmedf(-0+) zmedc(-0+)\n");
        for (i = 0; i < nbinm; i++)
        {
            if (nz[i] - nzc[i] > 0)
                fprintf(OUT, "%d %d %d %.1lf %.4lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n",
                        i, nz[i] - nzc[i], nz[i] - nzf[i], magbin[i] + 0.5*sizebinm,
                        (nz[i]*meanz[i] - nzc[i]*meanzf[i]) / (nz[i] - nzc[i]),
                        datmed[i][((int) ((nz[i] - nzf[i])*.5) + nzf[i])],
                        datmed[i][((int) ((nz[i] - nzf[i])*.25) + nzf[i])],
                        datmed[i][((int) ((nz[i] - nzf[i])*.75) + nzf[i])],
                        datmed[i][((int) ((nz[i] - nzc[i])*.5) + nzc[i])],
                        datmed[i][((int) ((nz[i] - nzc[i])*.25) + nzc[i])],
                        datmed[i][((int) ((nz[i] - nzc[i])*.75) + nzc[i])]
                       );

        }
        fclose(OUT);

    }



}

static void medstat(double *data, int ndata, double *med, double *q, double *tq)
{
    sortf(ndata, data, comp_asc);

    *med = data[((int) (ndata*.5))];
    *q = data[((int) (ndata*.4))];
    *tq = data[((int) (ndata*.6))];

}

static void st_ez(  int mode, struct galaxie imas[NASMAX], int nima,
                    double *z, double *drz, double ts[NASMAX][200], double theta[NASMAX][200],
                    double amp[NASMAX][200],
                    double pz[NASMAX][200], double ipz[NASMAX][200],
                    double kappa[NASMAX], double gam[NASMAX], double thetap[NASMAX],
                    double sumpz[NASMAX], double tauix[NASMAX],
                    double tauiy[NASMAX], int jmax )
{
    extern struct   g_mode M;
    extern  double  z_dlsds;

    register int    i, j, k;

    double   dzz = 0.05, zzmax = 5.;
    double   qi, tp, dp, ds;
    double  di, taui2, taui, tausx[200], tausy;
    double  ga1, ga2, gg, g2;
    double  signa;
    double  sum_pz, sumpznj[NASMAX];
    double  tausx2[200]; //,jacsp[200];
    double  **pznj;
    double  **jacn;
    struct  matrix  MA;

    char    name[20];
    FILE    *OUT;

    pznj = (double **) alloc_square_double(NASMAX, jmax);
    jacn = (double **) alloc_square_double(NASMAX, jmax);

    /* if fake sample rotate it by Pi/2 */

    if (mode == 1)
    {
        for (i = 0; i < nima; i++)
            imas[i].E.theta += M_PI / 2.;
    }

    /* computing kappa gamma thetap tauix tauiy for each arclet */

    for (i = 0; i < nima; i++)
    {
        z_dlsds = 0.3;

        MA = e_grad2_gal(&imas[i], NULL);
        MA.a *= z_dlsds;
        MA.b *= z_dlsds;
        MA.c *= z_dlsds;
        MA.d *= z_dlsds;
        kappa[i] = (MA.a + MA.c) / 2. / z_dlsds;
        ga1 = (MA.a - MA.c) / 2. / z_dlsds;
        ga2 = MA.b / z_dlsds;
        gam[i] = sqrt(ga1 * ga1 + ga2 * ga2);
        thetap[i] = 0.5 * atan2(ga2, ga1);
        qi = imas[i].E.b / imas[i].E.a;
        taui = (1. - qi * qi) / 2. / qi;
        tauix[i] = taui * cos(2.*(imas[i].E.theta - thetap[i]));
        tauiy[i] = taui * sin(2.*(imas[i].E.theta - thetap[i]));
        taui2 = tauix[i] * tauix[i] + tauiy[i] * tauiy[i];
        di = sqrt(1. + taui2);


        /* computing p(z) */

        sumpz[i] = 0.;
        sumpznj[i] = 0.;
        for (j = 0; (j < jmax) && (z[j] <= zzmax); j++)
        {
            z_dlsds = drz[j];

            gg = z_dlsds * gam[i] / (1. - z_dlsds * kappa[i]);
            g2 = gg * gg;
            signa = sgn(1 - g2);
            tp = 2 * gg / (1 - g2);
            dp = sqrt(1 + tp * tp);
            amp[i][j] = fabs( (1. - z_dlsds * kappa[i]) * (1. - z_dlsds * kappa[i]) -
                              z_dlsds * gam[i] * z_dlsds * gam[i] );
            tausy = signa * tauiy[i];
            tausx[j] = signa * tauix[i] * dp - tp * di;
            ts[i][j] = sqrt(tausx[j] * tausx[j] + tausy * tausy);
            ds = sqrt(1. + tausx[j] * tausx[j] + tausy * tausy);
            theta[i][j] = RTD * (thetap[i] + .5 * atan2(tausy, tausx[j]));

            pznj[i][j] = ptau(tausx[j], tausy) / ptau(0, tausy);
            if (j != 0)
                sumpznj[i] += (pznj[i][j] + pznj[i][j-1]) / 2.*dzz;
        };

        /* compute the jacobian (with a spline fit) and the integral of the
           probability, in order to normalize it
        */

        spline(z, tausx, j, 1e30, 1e30, tausx2);
        sp_dx(z, tausx, tausx2, j, jacn[i]);

        ipz[i][0] = 0.;
        for (k = 0; k < j; k++)
        {
            pz[i][k] = pznj[i][k] * fabs(jacn[i][k]);
            if (k != 0)
                ipz[i][k] = ipz[i][k-1] +
                            (pz[i][k] + pz[i][k-1]) / 2.*dzz;   /* trapez integration */
        }
        sumpz[i] = ipz[i][j-1];


    };    /* end of loop over the images */


    /* Normalizing probability and writing files if nima<50 */

    for (i = 0; i < nima; i++)
    {
        if (nima < 50)
        {
            NPRINTF(stderr, "WRITE: sz%s.dat\n", imas[i].n);
            sprintf(name, "sz%s.dat", imas[i].n);
            OUT = fopen(name, "w");
            fprintf(OUT, "#id  z  tau_s(z) Dmag(z) p(z) oldp(z) dtauIx/dz Ip(z)\n");

            sum_pz = 0.;
            for (j = 0; (j < jmax) && (z[j] <= zzmax); j++)
            {
                pz[i][j] /= sumpz[i];
                pznj[i][j] /= sumpznj[i];
                ipz[i][j] /= sumpz[i];

                fprintf(OUT, "%s %.3lf %.3lf %.2lf %.3lf %.3lf %.3lf %.3lf\n",
                        imas[i].n, z[j], ts[i][j], -2.5*log10(amp[i][j]),
                        pz[i][j], pznj[i][j], jacn[i][j], ipz[i][j]);
            };

            fclose(OUT);
        }
        else
        {
            sum_pz = 0.;
            for (j = 0; (j < jmax) && (z[j] <= zzmax); j++)
            {
                pz[i][j] /= sumpz[i];
                pznj[i][j] /= sumpznj[i];
                ipz[i][j] /= sumpz[i];
            };
        }
    }


    free_square_double(pznj, jmax);
    free_square_double(jacn, jmax);

} /* end of st_ez */


static void sp_dx(double *xx, double *yy, double *yy2, int nn, double *dydx)
{
    int k;
    double  deltax, deltay;

    deltax = deltay = 0;
    for (k = 0; k < nn - 1; k++)
    {
        deltax = xx[k+1] - xx[k];
        deltay = yy[k+1] - yy[k];

        dydx[k] = deltay / deltax - deltax / 3.*(yy2[k] + yy2[k+1] / 2.);
    }

    dydx[nn-1] = deltay / deltax + deltax / 3.*(yy2[nn-2] / 2. + yy2[nn-1]);

}
