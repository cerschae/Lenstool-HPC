#include<stdio.h>
#include<math.h>
#include<float.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"


static void w_clump(int nlens);


/****************************************************************/
/*      nom:        set_lens_par            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    set_lens_par(FILE *OUT)
{
    extern struct g_mode    M;
    extern struct g_source  S;
    extern struct g_grille  G;
    extern struct pot       lens[];
    const extern struct g_cosmo   C;

    double q; // elliptical parameter b/a_m
    double GG = 10.867;
    double invGG = 2.325968e-7; // 1/G in 10^12 Msol/Kpc/(km/s)^2
    int i;
    double rhos, c, m200;

    // Set the clumps parameters for use in Lenstool
    set_lens();

    //********************************************************************
    // Display the parameter values.
    //
    fprintf(OUT, "For z_s = %.4lf  DLS/DS:%.4lf\n", S.zs, lens[0].dlsds);
    fprintf(OUT, "DLS:%.4lf(lt) %.2lf(Mpc) ", distcosmo2(lens[0].z, S.zs), D0Mpc / C.h*distcosmo2(lens[0].z, S.zs));
    fprintf(OUT, "DOS:%.4lf(lt) %.2lf(Mpc)\n", distcosmo1(S.zs), D0Mpc / C.h*distcosmo1(S.zs));
    fprintf(OUT, "DOL:%.4lf(lt) %.2lf(Mpc) ", distcosmo1(lens[0].z), D0Mpc / C.h*distcosmo1(lens[0].z));
    fprintf(OUT, "DOL_lum:%.4lf(lt) %.2lf(Mpc)\n", dlumcosmo1(lens[0].z), D0Mpc / C.h*dlumcosmo1(lens[0].z));
    fprintf(OUT, "Mcrit:%e (10^12 Msol/kpc^2)\n", cH0_4piG * C.h / distcosmo1(lens[0].z) / lens[0].dlsds);

    fprintf(OUT, "Conversion Factor @ z = %lf, 1 arcsec == %.3lf kpc\n", lens[0].z, d0 / C.h*distcosmo1(lens[0].z));
    fprintf(OUT, "Number of Clumps: %ld\n", G.nlens);


    for (i = 0; i < G.nlens; i++)
    {
        switch (lens[i].type)
        {
            case(0):
                NPRINTF(stderr, "Clump %s: SIS:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: Singular Isothermal Sphere \n", lens[i].n);
                break;
            case(1):
                NPRINTF(stderr, "Clump %s: SIE:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: Elliptical Singular Isothermal Sphere \n", lens[i].n);
                break;
            case(2):
                NPRINTF(stderr, "Clump %s: ISC:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: Isothermal Sphere with core radius\n", lens[i].n);
                break;
            case(3):
                NPRINTF(stderr, "Clump %s: EISC:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: Elliptical Isothermal Sphere with core radius\n", lens[i].n);
                break;
            case(12):
                NPRINTF(stderr, "Clump %s: NFW:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: NFW \n", lens[i].n);
                break;
            case(13):
                NPRINTF(stderr, "Clump %s: Sersic:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: Sersic \n", lens[i].n);
                break;
            case(14):
                NPRINTF(stderr, "Clump %s: External shear:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: External shear \n", lens[i].n);
                break;
  	    case(15):
  	    	NPRINTF(stderr, "Clump %s: Einasto \n", lens[i].n);
  	    	fprintf(OUT, "-------- Clump %s: Einasto \n", lens[i].n);
  	    	break;
            case(16):
                NPRINTF(stderr, "Clump %s: Hernquist:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: Hernquist \n", lens[i].n);
                break;
            case(-1):
                NPRINTF(stderr, "Clump %s: True Elliptical SIS:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: True Elliptical SIS \n", lens[i].n);
            case(-2):
                NPRINTF(stderr, "Clump %s: True Elliptical BBS model:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: True Elliptical BBS model \n", lens[i].n);
                break;
            case(7):
                NPRINTF(stderr, "Clump %s: Point Masse:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: Point Masse \n", lens[i].n);
                break;
            case(9):
                NPRINTF(stderr, "Clump %s: Plan Masse:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: Plan Masse \n", lens[i].n);
                break;
            case(5):
                NPRINTF(stderr, "Clump %s: Hubble Modified Law:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: Hubble Modified Law \n", lens[i].n);
                break;
            case(8):
                NPRINTF(stderr, "Clump %s: PIEMD Kovner:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: PIEMD Kovner\n", lens[i].n);
                break;
            case(81):
                lens[i].masse = 1.5 * M_PI * invGG * lens[i].sigma * lens[i].sigma *
                                lens[i].rcutkpc * lens[i].rcutkpc / (lens[i].rckpc + lens[i].rcutkpc);

                NPRINTF(stderr, "Clump %s: trunc. PIEMD Kovner:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: PIEMD Kovner, truncated\n", lens[i].n);
                fprintf(OUT, " Total mass: %lf(cor) %lf (10^12 M_sol)\n",
                        4*M_PI / 3*M_PI / GG*(lens[i].sigma / 1000)*
                        (lens[i].sigma / 1000)*lens[i].rcut*(d0 / C.h*distcosmo1(lens[i].z)),
                        lens[i].masse );
                fprintf(OUT, " rcut:%.2lf(\") %.2lf(kpc)\n",
                        lens[i].rcut, lens[i].rcut*(d0 / C.h*distcosmo1(lens[i].z)) );
                break;

            case(82):
                NPRINTF(stderr, "Clump %s: PIEMD Kovner, shallow center:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: PIEMD Kovner, shallow center\n", lens[i].n);
                fprintf(OUT, " Steep radius:%.2lf(\") %.2lf(kpc)\n",
                        lens[i].rcut, lens[i].rcut*(d0 / C.h*distcosmo1(lens[i].z)) );
                break;
            case(83):
                NPRINTF(stderr, "Clump %s: EMD Kovner, 3/2:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: PIEMD Kovner, shallow center\n", lens[i].n);
                fprintf(OUT, " Steep radius:%.2lf(\") %.2lf(kpc)\n",
                        lens[i].rcut, lens[i].rcut*(d0 / C.h*distcosmo1(lens[i].z)) );
                break;
            case(84):
                NPRINTF(stderr, "Clump %s: EMD Kovner, 0.5a-0.5s+1.5s:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: PIEMD Kovner, 0.5a-0.5s+1.5s\n", lens[i].n);
                fprintf(OUT, " Steep radius:%.2lf(\") %.2lf(kpc)\n",
                        lens[i].rcut, lens[i].rcut*(d0 / C.h*distcosmo1(lens[i].z)) );
                break;
            case(85):
                NPRINTF(stderr, "Clump %s: EMD Kovner, 1:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: PIEMD Kovner, 1a\n", lens[i].n);
                fprintf(OUT, " Steep radius:%.2lf(\") %.2lf(kpc)\n",
                        lens[i].rcut, lens[i].rcut*(d0 / C.h*distcosmo1(lens[i].z)) );
                break;
            case(86):
                NPRINTF(stderr, "Clump %s: EMD Kovner, 1a-1s:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: PIEMD Kovner, 1a-1s\n", lens[i].n);
                fprintf(OUT, " Steep radius:%.2lf(\") %.2lf(kpc)\n",
                        lens[i].rcut, lens[i].rcut*(d0 / C.h*distcosmo1(lens[i].z)) );
                break;
            case(87):
                NPRINTF(stderr, "Clump %s: EMD Kovner, 1a-1s+0.5a-0.5s:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: PIEMD Kovner, 1a-1s+0.5a-0.5s\n", lens[i].n);
                fprintf(OUT, " Steep radius:%.2lf(\") %.2lf(kpc)\n",
                        lens[i].rcut, lens[i].rcut*(d0 / C.h*distcosmo1(lens[i].z)) );
                break;
            case(88):
                NPRINTF(stderr, "Clump %s: EMD Kovner, 1a-1s+1.5s:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: PIEMD Kovner, 1a-1s+1.5s\n", lens[i].n);
                fprintf(OUT, " Steep radius:%.2lf(\") %.2lf(kpc)\n",
                        lens[i].rcut, lens[i].rcut*(d0 / C.h*distcosmo1(lens[i].z)) );
                break;
            case(89):
                NPRINTF(stderr, "Clump %s: EMD Kovner, 1a-1s+0.5a-0.5s:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: PIEMD Kovner, 1a-1s+0.5a-0.5s\n", lens[i].n);
                fprintf(OUT, " Steep radius:%.2lf(\") %.2lf(kpc)\n",
                        lens[i].rc*lens[i].beta, lens[i].rc*lens[i].beta*(d0 / C.h*distcosmo1(lens[i].z)) );
                break;
            case(10):
                NPRINTF(stderr, "Clump %s: Spline Potential:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: Spline Potential\n", lens[i].n);
            default:
                NPRINTF(stderr, "Clump %s: Pseudo-Elliptical Potential with Core Radius:", lens[i].n);
                fprintf(OUT, "-------- Clump %s: Pseudo-Elliptical Potential with Core Radius\n", lens[i].n);
                break;
        }

        /*
        *  elliptical parameters q, just to be printed
        */
        if ( lens[i].type == 0 || lens[i].type == 2 )
            q = 1.;
        else if ( lens[i].type ==  8  ||
                  lens[i].type == -2  ||
                  ( lens[i].type > 80 && lens[i].type < 90 ) )
            q = (1. - lens[i].epot) / (1. + lens[i].epot);
        else
            q = sqrt((1. - lens[i].epot) / (1. + lens[i].epot));

        switch ( lens[i].type )
        {
            case(14):
                fprintf(OUT, " gamma:%.4lf\n", lens[i].emass);
                break;
            case(12):
                fprintf(OUT, " e_m:%.4lf b/a_m:%.4lf e_p:%.4lf\n",
                        lens[i].emass, q, lens[i].epot);
                fprintf(OUT, " sigma_s:%.2lf(km/s) b0:%.4lf\n", lens[i].sigma, lens[i].b0);
                fprintf(OUT, " rs:%.2lf(\") %.2lf(kpc)\n", lens[i].rc, lens[i].rckpc);

                if ( lens[i].rcut != DBL_MAX )
                    fprintf(OUT, " r200:%.2lf(\") %.2lf(kpc)\n", lens[i].rcut, lens[i].rcutkpc);

                e_nfw_rs2c(lens[i].sigma, lens[i].rckpc, &rhos, &c, &m200, lens[i].z);
                fprintf(OUT, " rhos=%.1le c=%.1lf M200=%.1le\n", rhos, c, m200);
                NPRINTF(stderr, " rhos=%.1le c=%.1lf M200=%.1le e_m=%.3lf r_ct=%.1lf r_cr=%.1lf\n",
                        rhos, c, m200, lens[i].emass, lens[i].ct, lens[i].cr);
                break;
            case(13):
                fprintf(OUT, " sigmae:%.1le (Msol) b0:%.4lf\n", lens[i].sigma, lens[i].b0);
                fprintf(OUT, " Re:%.2lf(\") %.2lf(kpc)\n", lens[i].rc, lens[i].rckpc);
                fprintf(OUT, " n:%.2lf\n", lens[i].alpha);
                NPRINTF(stderr, " sigmae=%.1le Re=%.1lf n=%.1lf e_m=%.3lf\n",
                        lens[i].sigma, lens[i].rc, lens[i].alpha, lens[i].emass);
                break;
  	    case(15): //einasto
  	    	fprintf(OUT,"rhos/%.1le (Msol) b0:%.4lf\n",lens[i].pmass,lens[i].b0);
  	    	fprintf(OUT,"Rs:%.2lf(\")%.2lf(kpc)\n",lens[i].rc,lens[i].rckpc);
  	    	fprintf(OUT,"n/%.2lf\n",lens[i].alpha);
  	    	NPRINTF(stderr,"rhos=%.1le Rs=%.1lf n=%.1lf e_m=%.3lf\n",lens[i].pmass,lens[i].rc,lens[i].alpha,lens[i].emass);
  	    	break;
            case(16):
                fprintf(OUT, " e_m:%.4lf b/a_m:%.4lf e_p:%.4lf\n",
                        lens[i].emass, q, lens[i].epot);
                fprintf(OUT, " sigma_s:%.2lf(km/s) b0:%.4lf\n", lens[i].sigma, lens[i].b0);
                fprintf(OUT, " rs:%.2lf(\") %.2lf(kpc)\n", lens[i].rc, lens[i].rckpc);
                NPRINTF(stderr, " sigma_s=%.1le rs=%.2lf(\") e_m=%.3lf\n",
                        lens[i].sigma, lens[i].rc, lens[i].emass);
                break;
           case(9):
                fprintf(OUT, " pmass:%.2lf (g/cm2) b0:%.2lf \n", lens[i].pmass, lens[i].b0);
                NPRINTF(stderr, " pmass:%.2lf (g/cm2) b0:%.2lf \n", lens[i].pmass, lens[i].b0);
                break;
            default:
                fprintf(OUT, " vdisp:%.2lf(km/s) b0:%.4lf\n", lens[i].sigma, lens[i].b0);
                fprintf(OUT, " rc:%.2lf(\") %.2lf(kpc)\n", lens[i].rc, lens[i].rckpc);

                if ( lens[i].rcut != DBL_MAX )
                {
                    fprintf(OUT, " rt:%.2lf(\") %.2lf(kpc)\n", lens[i].rcut, lens[i].rcutkpc);
                    NPRINTF(stderr, " vdisp=%.0lf rc=%.1lf rt=%.1lf e_m=%.3lf r_ct=%.1lf r_cr=%.1lf\n",
                            lens[i].sigma, lens[i].rc, lens[i].rcut, lens[i].emass, lens[i].ct, lens[i].cr);
                }
                else
                {
                    NPRINTF(stderr, " vdisp=%.0lf rc=%.1lf rt=%.1lf e_m=%.3lf r_ct=%.1lf r_cr=%.1lf\n",
                            lens[i].sigma, lens[i].rc, 0., lens[i].emass, lens[i].ct, lens[i].cr);
                }
                fprintf(OUT, " r_ct:%.4lf r_cr:%.4lf\n", lens[i].ct, lens[i].cr);
                break;
        }
    }

    if (M.verbose > 0)
        w_clump(G.nlens);
}

/* Write the clump parameters values in clump.dat and optionally on stderr.
 * The a and b elliptical semi axis are in arcsec such that :
 * - b is rcore if rcore > 5, otherwise b is rcut
 * - a is b / q
 */
static void w_clump(int nlens)
{
    extern struct g_mode  M;
    extern struct pot       lens[];
    register int  i;
    double q, a, b;
    FILE    *CLUMP;


    CLUMP = fopen("clump.dat", "w");

    if ( M.iref != 2 )
        fprintf(CLUMP, "#REFERENCE 3 %.7lf %.7lf\n", M.ref_ra, M.ref_dec);
    else
        fprintf(CLUMP, "#REFERENCE 2 %.7lf %.7lf\n", M.ref_ra, M.ref_dec);

    for (i = 0; i < nlens; i++)
    {
        q = sqrt((1 - lens[i].emass) / (1 + lens[i].emass));
        if ( lens[i].rc > 5 && lens[i].rcut != DBL_MAX )
            b = lens[i].rc;
        else
            b = lens[i].rcut;

        a = b / q;
        if ( lens[i].rcut != DBL_MAX )
            q = lens[i].rcutkpc;
        else
            q = 0;

        fprintf(CLUMP,
                "%2d %7.2lf %7.2lf %6.2lf %6.2lf %7.3lf %7.2lf %7.3lf %7.3lf %6.1lf %5.3lf %6.2lf %s %d\n",
                lens[i].type, lens[i].C.x, lens[i].C.y,
                a, b, lens[i].theta*RTD, lens[i].emass, lens[i].rckpc,
                q, lens[i].sigma, lens[i].z, lens[i].mag, lens[i].n, i);

        if ( lens[i].type == 12 )
        {
            NPRINTF(stderr, "%s %d mag:%.2lf sig:%.2le L:%.3lf M/L:%.2lf c:%.1lf m200:%.2le\n",
                    lens[i].n, 12, lens[i].mag, lens[i].sigma,
                    lens[i].lum, lens[i].mtol, lens[i].beta, lens[i].masse);
        }
        else
        {
            NPRINTF(stderr, "%s %d mag:%.2lf sig:%.2lf L:%.3lf M/L:%.2lf rc:%.3lf rt:%.2lf\n",
                    lens[i].n, lens[i].type, lens[i].mag, lens[i].sigma,
                    lens[i].lum, lens[i].mtol, lens[i].rckpc, q);
        }
    };
    fclose(CLUMP);
}

