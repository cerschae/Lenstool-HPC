#include<stdio.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        st_opt          */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    st_opt( int mode,
                struct galaxie *imas,
                int nima,
                double *z,
                double ts[NASMAX][200],
                double pz[NASMAX][200],
                double ipz[NASMAX][200],
                double amp[NASMAX][200],
                double *kappa,
                double *gam,
                double *thetap,
                double *sumpz,
                double *tauix,
                double *tauiy,
                int jmax,
                double *zpm,
                double *nmag )
{
    const extern  struct  g_mode  M;
    const extern  struct  pot lens[];
    extern  double  z_dlsds;

    struct  galaxie     source;
    //struct    ellipse ampli;

    char    nsource[32], nimage[32], namez[32], nselect[32],
    nmult[32], ncomp[32], nttt[32], ncip[32];
    char    ext[4];

    register int    i, j;
    double  dmdz[NASMAX];

    double  peakmax, peak[5];
    double  taus, tausx, tausy;
    double  ts_mx = 0., ts_my = 0., ts_sx = 0., ts_sy = 0.;

    double  pzmax, zpzmax, zpzmax1, zpzmax2;
    int npz, jpzmax, jpzmax1, jpzmax2;
    double  t0, t1;
    double  zTmin, zTmin1, zTmin2;
    double  zTmean1, zTmean2;
    double  Tmin;
    int jTmin, jTmin1, jTmin2;
    int jTmean1, jTmean2;

    double  depsi, dthetai, epsi, gpot, dpot;
    double  pkratio = 0.5, zzmin = 0.3, zzmax = 5.0;
    int jzmax;
    double  z2[200];
    double  ci25, ci50, ci75;

    FILE    *OUTT, *COMP, *COMPT, *OUT_S, *OUT_I, *OUT_Z, *OUTC; //*OUT,
    FILE    *SELECT, *MULT;

    /*
    * Automatic search of the most probable redshift
    * using the z-probability
    */

    if (mode == 0)
        strcpy(ext, "zopt");
    else
        strcpy(ext, "zofk");

    strcpy(ncip, "cip.");
    strcat(ncip, ext);
    strcpy(nsource, "source.");
    strcat(nsource, ext);
    strcpy(nimage, "image.");
    strcat(nimage, ext);
    strcpy(namez, "zz.");
    strcat(namez, ext);
    strcpy(nselect, "select.");
    strcat(nselect, ext);
    strcpy(nmult, "mult.");
    strcat(nmult, ext);
    strcpy(ncomp, "comp.");
    strcat(ncomp, ext);
    strcpy(nttt, "ttt.");
    strcat(nttt, ext);

    OUTC = fopen(ncip, "w");
    OUTT = fopen(nttt, "w");
    OUT_S = fopen(nsource, "w");
    OUT_I = fopen(nimage, "w");
    OUT_Z = fopen(namez, "w");
    COMP = fopen(ncomp, "w");
    COMPT = fopen("compt.dat", "w");
    SELECT = fopen(nselect, "w");

    if (mode == 0)
        MULT = fopen(nmult, "w");

    fprintf(OUT_Z, "#na\t zm- zmin zm+ \tpz_max npz mag_true depsi dthetai depsi+dthetai \n");
    fprintf(COMP, "#na  kappa gamma tauix tauiy zmin zopt zmax err_z pb_max n_peak true_er_z true_z\n");
    fprintf(COMPT, "#na  kappa gamma zmin zopt zmax err_z eps_s n_peak true_er_z true_z\n");

    for (i = 0; i < nima; i++)
    {

        /* recherche de l'ellipticite min (Tmin,zTmin) */

        pzmax = pz[i][1];
        zpzmax = z[1];
        jpzmax = 1;
        npz = 0;
        peakmax = 0.;
        Tmin = fabs(ts[i][1]);
        zTmin = z[1];
        jTmin = 1;

        for (j = 1; (j < jmax - 1) && (z[j] <= zzmax); j++)
        {
            if (fabs(ts[i][j]) < Tmin)
            {
                Tmin = fabs(ts[i][j]);
                zTmin = z[j];
                jTmin = j;
            }


            /* peak finder, find a peak in the probability and keep the narrowest one */

            if ( ((pz[i][j] > pz[i][j-1]) && (pz[i][j] > pz[i][j+1]))
                 || ((j == 1) && (pz[i][1] < pz[i][0]))
                 || ((j == jmax - 2) && (pz[i][jmax-1] > pz[i][jmax-2])) )
            {
                if (j > 1)
                {
                    peak[npz] = 0.5 * (2.*pz[i][j] - pz[i][j-1] - pz[i][j+1]);
                }
                else if (j == 1)
                {
                    peak[npz] = pz[i][0] - pz[i][1];
                }
                else
                {
                    peak[npz] = pz[i][jmax-1] - pz[i][jmax-2];
                }

                npz++;

                if (npz == 1)
                {
                    peakmax = peak[npz-1];
                }
                else if (peak[npz-1] > peakmax)
                {
                    peakmax = peak[npz-1];
                    pzmax = pz[i][j];
                    zpzmax = z[j];
                    jpzmax = j;
                }
            }

            if (pz[i][j] > pzmax)
            {
                pzmax = pz[i][j];
                zpzmax = z[j];
                jpzmax = j;
            }
        }

        /* interpolate the confidence interval from the ipz array */

        for (j = 1; (j < jmax - 1) && (z[j] <= zzmax); j++);

        jzmax = j;

        spline(ipz[i], z, jzmax, 1e30, 1e30, z2);

        ci25 = interpol(0.25, ipz[i], z, z2, jzmax);
        ci50 = interpol(0.50, ipz[i], z, z2, jzmax);
        ci75 = interpol(0.75, ipz[i], z, z2, jzmax);

        fprintf(OUTC, "%d (%s) %.3lf %.3lf %.3lf  mintau:%.3lf maxpz%.3lf \n",
                jzmax, imas[i].n, ci25, ci50, ci75, zTmin, zpzmax);


        /* determination de l'ecart */
        for (j = jpzmax + 1; ((j < jmax) && (z[j] <= zzmax) && (pz[i][j] > pzmax*pkratio)); j++)
        {
            NPRINTF(stderr, "%s %d %.3lf %.3lf %.3lf\n", imas[i].n, j, z[j], pz[i][j], pzmax);
        }

        if (jpzmax + 1 == jmax || z[jpzmax+1] >= zzmax )
        {
            zpzmax2 = zpzmax;
            jpzmax2 = jpzmax;
        }
        else
        {
            zpzmax2 = z[j];
            jpzmax2 = j;
        }

        for (j = jpzmax - 1; ((j > 0) && (pz[i][j] > pzmax*pkratio)); j--);

        zpzmax1 = z[j];
        jpzmax1 = j;

        /*
        * most probable redshift
        * new magnitude
        * magnitude gradent error
        */
        zpm[i] = zpzmax;
        nmag[i] = imas[i].mag - 2.5 * log10(amp[i][jpzmax]);
        dmdz[i] = 2.5 * (log10(amp[i][jpzmax1]) - log10(amp[i][jpzmax2])) / (zpzmax2 - zpzmax1);

        t0 = sqrt(.2 * .2 + Tmin * Tmin);
        t1 = t0 + .17;

        for (j = jTmin; ((j < jmax) && (z[j] <= zzmax) && (z[j+1] != 0.) && (fabs(ts[i][j]) < t0)); j++);

        zTmean2 = z[j];
        jTmean2 = j;

        for (j = jTmin; ((j > 0) && (fabs(ts[i][j]) < t0)); j--);

        zTmean1 = z[j];
        jTmean1 = j;

        for (j = jTmin; ((j < jmax) && (z[j] <= zzmax) && (z[j+1] != 0.) && (fabs(ts[i][j]) < t1)); j++);

        zTmin2 = z[j];
        jTmin2 = j;

        for (j = jTmin; ((j > 0) && (fabs(ts[i][j]) < t1)); j--);

        zTmin1 = z[j];
        jTmin1 = j;


        if ((mode == 0) && (npz > 2))
            fprintf(MULT, "%s %.4lf %.4lf\t%.4lf %.4lf %.4lf\t%.3lf %.3lf %.3lf %.3lf\n", imas[i].n,
                    imas[i].C.x, imas[i].C.y,
                    imas[i].E.a, imas[i].E.b, imas[i].E.theta*RTD,
                    zpzmax, imas[i].mag - 2.5*log10(amp[i][jpzmax]), -2.5*log10(amp[i][jpzmax]), pzmax);

        z_dlsds = dratio(lens[0].z, zpzmax);
        source = unlens1(imas[i], z_dlsds);
        taus = 0.5 * (source.E.a / source.E.b - source.E.b / source.E.a);
        tausx = taus * cos(2.*source.E.theta);
        tausy = taus * sin(2.*source.E.theta);

        if ((source.E.a != 0.) && (source.E.b != 0.))
        {
            ts_mx += tausx;
            ts_my += tausy;
            ts_sx += tausx * tausx;
            ts_sy += tausy * tausy;
        }

        epsi = (imas[i].E.a * imas[i].E.a - imas[i].E.b * imas[i].E.b) /
               (imas[i].E.a * imas[i].E.a + imas[i].E.b * imas[i].E.b);
        gpot = z_dlsds * gam[i] / (1 - z_dlsds * kappa[i]);
        dpot = (1 + gpot * gpot) / (1 - gpot * gpot);

        fprintf(OUTT, "%s %.3lf %.3lf %.3lf %.3lf\n", imas[i].n, zpzmax,
                epsi, (1 - z_dlsds*kappa[i])*dpot, (1 - z_dlsds*kappa[i])*imas[i].var1 / epsi);

        if (epsi != 0.0)
            depsi = (1 - z_dlsds * kappa[i]) * dpot * imas[i].var1 / epsi;
        else
            depsi = 9.;

        dthetai = fabs(tan(2 * (imas[i].E.theta - thetap[i])) * imas[i].var2);

        fprintf(OUT_I, "%s %.4lf %.4lf\t%.4lf %.4lf %.4lf\t%.3lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n",
                imas[i].n, imas[i].C.x, imas[i].C.y,
                imas[i].E.a, imas[i].E.b, imas[i].E.theta*RTD,
                zpzmax, imas[i].mag - 2.5*log10(amp[i][jpzmax]),
                -2.5*log10(amp[i][jpzmax]), imas[i].mu, pzmax, tauix[i], tauiy[i],
                imas[i].E.a / imas[i].E.b,
                depsi, dthetai, depsi + dthetai);

        fprintf(OUT_S, "%s %.4lf %.4lf\t %.4lf %.4lf %.4lf\t%.4lf %.4lf %.4lf %.4lf\n",
                imas[i].n, source.C.x, source.C.y,
                source.E.a, source.E.b, source.E.theta*RTD,
                zpzmax, imas[i].mag - 2.5*log10(amp[i][jpzmax]), tausx, tausy);

        if ((zpzmax > zzmin) && (zpzmax < zzmax) && (tauix[i] > 0.))
            fprintf(SELECT, "%s %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf %.2lf %.2lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n",
                    imas[i].n, imas[i].C.x, imas[i].C.y,
                    imas[i].E.a, imas[i].E.b, imas[i].E.theta*RTD,
                    zpzmax, imas[i].mag - 2.5*log10(amp[i][jpzmax]),
                    -2.5*log10(amp[i][jpzmax]), imas[i].mu, pzmax, tauix[i], tauiy[i],
                    imas[i].E.a / imas[i].E.b,
                    depsi, dthetai, depsi + dthetai);

        fprintf(OUT_Z, "%s\t%.4lf %.4lf %.4lf\t%.4lf %d  %.3lf %.3lf %.3lf %.3lf\n",
                imas[i].n, zpzmax1, zpzmax, zpzmax2, pzmax, npz,
                imas[i].mag - 2.5*log10(amp[i][jpzmax]), depsi, dthetai,
                depsi + dthetai);

        fprintf(COMP, "%s %.3lf %.3lf %.2lf %.2lf %.4lf %.4lf %.4lf %.4lf %.4lf %d %.4lf %.4lf %.2lf %.2lf\n",
                imas[i].n, kappa[i], gam[i], tauix[i], tauiy[i],
                zpzmax1, zpzmax, zpzmax2, (zpzmax2 - zpzmax1) / 2.,
                pzmax, npz, fabs(zpzmax - imas[i].z), imas[i].z, dmdz[i],
                imas[i].mag - 2.5*log10(amp[i][jpzmax]));

        fprintf(COMPT, "%s %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %d %.4lf %.4lf\n",
                imas[i].n, kappa[i], gam[i],
                zTmean1, zTmin, zTmean2, (zTmean2 - zTmean1) / 2.,
                source.E.a / source.E.b, npz, fabs(zTmin - imas[i].z), imas[i].z);

    }

    ts_mx /= i;
    ts_my /= i;
    ts_sx /= i;
    ts_sy /= i;
    ts_sx -= ts_mx * ts_mx;
    ts_sy -= ts_my * ts_my;

    NPRINTF(stderr, "RES: mean tsx:%.4lf tsy:%.4lf disp tsx:%.4lf tsy:%.4lf\n",
            ts_mx, ts_my, sqrt(ts_sx), sqrt(ts_sy) );


    fclose(OUTT);
    fclose(OUT_S);
    fclose(OUT_I);
    fclose(OUT_Z);
    fclose(COMP);
    fclose(COMPT);
    fclose(SELECT);

    if (mode == 0)
        fclose(MULT);
}

