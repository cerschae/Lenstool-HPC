#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<float.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#undef CHIRES

static double chi2SglImage( struct galaxie *pima, struct point *ps );
static int chi2_img( double *chi2, double *lh0 );
static int chi2_src( double *chi2, double *lh0, double *np_b0 );
static void chi2_wl( double *chi2, double *lh0 );
static void srcfitLinFit( double *chi2, double *lh0 ) __attribute__((noinline));

/********************************************************/
/*      fonction: o_chi             */
/*      auteur: jpk             */
/********************************************************
 * Return the total chi2 value
 *
 * Global variables used :
 * - chip, chia, chix, chiy, chil, chis, chi_im, M, G, I, lens, shm, nshmap, arclet
 *   multi, cl, narclet, nwarn, optim_z, elix, SC, amplifi_mat, amplifi_matinv
 *   map_p, map_axx, map_ayy
 * - in sig2posS() : amplifi
 * - in sig2posSj() : amplifi
 * - in sig2posS4j() : G, lens, lens_table
 * - in distcosmo1() : C
 * - in e_zeroamp() : G, lens, lens_table
 * - in weight_baryc() : amplifi, G, lens, lens_table, multi, C
 * - in fmin_ell() : elix, SC, G, lens, lens_table
 * - in o_flux() :   multi, G, lens, lens_table
 * - in o_chi_flux() : I
 * - in amplif() : I, amplifi, G, lens, C, lens_table
 * - in amplif_mat() : I, multi, amplifi_mat, G, lens, lens_table, C
 * - in amplif_matinv() : I, multi, amplifi_matinv, G, lens, lens_table, C
 * - in dratio() : C
 * - in e_unmag() : G, lens, lens_table
 * - in chi_invim() : pi, ps, M, iline, radial, tangent, nrline, ntline,
 *                    CL, lens, flagr, flagt, G, lens_table
 * - in sp_set() : v_xx, v_yy
 * - in o_dpl() : G, lens, lens_table
 * - in o_mag_m() : G, lens, lens_table
 * - in unlens_bc() : nwarn, distmin, G, lens, lens_table
 */
#ifdef CHIRES
FILE *OUT;

void o_chires(char *filename)
#else
// Return 1 if error, 0 otherwise
int o_chi_lhood0(double *chi_ext, double *lhood0_ext, double *np_b0)
#endif
{
    /* variables externes */
    extern struct g_mode    M;
    const extern struct g_grille  G;
    const extern struct g_image   I;
	extern struct g_dyn     Dy;       //TV
    const extern struct pot       lens[NLMAX];
    extern struct shear     shm[];
    extern struct cline     cl[];
    extern double chil, chis, chi_im, chi_vel, chi_mass;   //TV
    extern int nshmap;

//  extern double amplifi[NFMAX][NIMAX];

    extern double **map_p, **map_axx, **map_ayy;
    extern int nplo, **imuo;
    extern double drima, **imo, **wo, **ero, **soo;
    extern struct pixlist plo[];

    /* local variables */
    struct point A, B, D;
    struct ellipse ampli;
    double chi_mult, chisrcfit, chish;
    double chi, lhood0;
    double x, qp, tp;
    int i, j;

#ifdef CHIRES
    OUT = fopen(filename, "w");
    double *np_b0 = NULL;
#endif

    /*variables initialisation*/
    chi = 0.;
    chish = 0.;
    chi_im = 0.;
    chis = 0.;
    chi_mult = 0.;
    chisrcfit = 0.;
    chil = 0.;
    lhood0 = 0.;
	chi_vel = 0.;
	chi_mass = 0.;

    /* spline mapping */
    if (lens[0].type == 10)
    {
        sp_set(map_p, G.nx, G.ny, map_axx, map_ayy);
        /*
            wr_pot("p.ipx",map_p);
            wr_mass("m.ipx",map_axx,map_ayy);
        */
    }

    /* if cleanset is true in cleanlens*/
    if (M.iclean != 0)
    {
        if (M.iclean == 2)
        {
            // Add up barycenter position to source positions
            extern struct g_source S;
            extern struct galaxie source[NFMAX];
            extern struct galaxie multi[NFMAX][NIMAX];
            const extern struct g_image  I;
            struct point Bs[NFMAX];  // list of source barycenters  
            struct point Ps[NIMAX];
            char str[IDSIZE], *pch;  // get image identifier

            if( I.n_mult > 0 )
            {
                for( i = 0; i < S.ns; i++ )
                {
                    // check if source is attached to an image
                    strcpy(str, source[i].n);
                    pch = strtok(str, "_");
                    pch = strtok(NULL, "_");
                    if( pch == NULL )    break;  // no attachment to an image

                    // look for the corresponding image family
                    j = 0; 
                    while ( indexCmp( pch, multi[j][0].n ) && j < I.n_mult ) j++;

                    // add barycenter position
                    if( j < I.n_mult )
                    {
                        o_dpl(I.mult[j], multi[j], Ps, np_b0);
                        Ps[I.mult[j]] = Ps[0];
                        Bs[i] = bcentlist(Ps, I.mult[j]);
                        source[i].C.x += Bs[i].x;
                        source[i].C.y += Bs[i].y;
                    }
                }
            }

            // Compute image
            extern struct g_pixel imFrame;
            extern struct galaxie source[NFMAX];
            const extern struct g_observ  O;
        
            double dx = imFrame.xmax - imFrame.xmin;
            double dy = imFrame.ymax - imFrame.ymin;

            int verbose_save = M.verbose;
            M.verbose = 0;
            
	    //2Eric: You assume that imFrame.pixelx == imFrame.pixely ?	   
	    if (fabs(( imFrame.pixelx - imFrame.pixely)/imFrame.pixelx) > 1e-4)
	     {
		fprintf(stderr, "Error | imFrame.pixelx(%f) != imFrame.pixely(%f) but we assume that they are equal\n", imFrame.pixelx, imFrame.pixely);
		exit(EXIT_FAILURE);
	     }
	    
	    o_pixel(ero, imFrame.nx, imFrame.ny, imFrame.pixelx, imFrame.xmin, imFrame.xmax, imFrame.ymin, imFrame.ymax, source, map_axx, map_ayy);
	   	   
        //wrf_fits_abs("simu.fits", ero, imFrame.nx, imFrame.ny, imFrame.xmin, imFrame.xmax, imFrame.ymin, imFrame.ymax, M.ref_ra, M.ref_dec);

            if (O.setseeing)
                d_seeing_omp(ero, imFrame.nx, imFrame.ny, imFrame.pixelx);

            if (O.bruit)
                d_bruiter_omp(ero, imFrame.nx, imFrame.ny);

            M.verbose = verbose_save;
	
            if (wo != NULL)
            {
                for (i = 0; i < imFrame.ny; i++ )
                    for (j = 0; j < imFrame.nx; j++ )
                    {
                        ero[i][j] -= imo[i][j];
                        chi_im += ero[i][j] * ero[i][j] * wo[i][j];
                        lhood0 -= log( 2.*M_PI*wo[i][j] );
                    }
            }
            else
            {
                for (i = 0; i < imFrame.ny; i++ )
                    for (j = 0; j < imFrame.nx; j++ )
                    {
                        ero[i][j] -= imo[i][j];
                        chi_im += ero[i][j] * ero[i][j] / fabs(imo[i][j] + 1.);
                        lhood0 += log( 2.*M_PI*fabs(imo[i][j] + 1.) );
                    }
            }

            // restore relative source positions
            if( I.n_mult > 0 )
                for( i = 0; i < S.ns; i++ )
                {
                    // check if source is attached to an image
                    strcpy(str, source[i].n);
                    pch = strtok(str, "_");
                    pch = strtok(NULL, "_");
                    if( pch == NULL )    break;  // no attachment to an image

                    // look for the corresponding image family
                    j = 0; 
                    while ( indexCmp( pch, multi[j][0].n ) && j < I.n_mult ) j++;

                    // add barycenter position
                    if( j < I.n_mult )
                    {
                        source[i].C.x -= Bs[i].x;
                        source[i].C.y -= Bs[i].y;
                    }
                }

        }
        else 

        /* chi_im is the error associated to the transformation image -> source*/
            chi_im = chi_invim(imo, plo, nplo, drima, soo, ero, imuo);
#ifdef CHIRES
        fprintf(OUT, "chi_im %.2lf\n", chi_im);
#endif
    }

    /* ----------------------------------------------------
     * Weak-Lensing chi2
     * if arcletstat is true in image section*/
    if (I.stat != 0)
    {
        chi2_wl(&chis, &lhood0);

#ifndef CHIRES
        // check for NaN. It can happen in some very particular conditions
        // - when for 1 arclet, grad2.a = -5.2e-08 (EJ 16/03/2011) 
        if( chis != chis ) return(1);
#endif
    }


    /* if shearmap is true in image section*/
    if (I.shmap != 0)
    {
        chish = 0;
        for (i = 0; i < nshmap; i++)
        {
            ampli = e_unmag(&shm[i].C, I.dl0ssh, I.dossh, I.zsh);
            qp = fabs(ampli.b / ampli.a);
            tp = fabs(1. - qp * qp) / 2. / qp;
            ampli.theta += M_PI / 2.;
            x = (shm[i].mx - tp * cos(2.*ampli.theta)) / shm[i].err;
            chish += x * x;
            x = (shm[i].my - tp * sin(2.*ampli.theta)) / shm[i].err;
            chish += x * x;
        }

        chish /= nshmap;
        chis += chish;
    } /*end if (I.shmap!=0)*/

    /* if I.n_mult is true in image keyword*/
    if (I.n_mult != 0)
    {
        if (I.forme >= 0)
            // Source plan chi2
            x = chi2_src( &chi_mult, &lhood0, np_b0 );
        else
        {
            // Image plane chi2
            tp = lhood0;
            x = chi2_img( &chi_mult, &lhood0 );
            /*
            //Uncomment if you want to chain image plane and source plane
                        if( x != 0.  )
                        {
                            chi_mult = 0.; lhood0 = tp;
                            chi2_src(&chi_mult, &lhood0 );
                        }
            */
        }
#ifndef CHIRES
        // This case should never happens in o_chires().
        // this #ifndef condition aims only at avoiding the compilation warning
        if ( x != 0 ) return(1);
#endif

    } /*end of optimization with the arclets*/


    // if there are points for source plane fitting
    if (I.srcfit != 0 )
    {
        srcfitLinFit(&chisrcfit, &lhood0);
#ifndef CHIRES
        if ( chisrcfit == -1. ) return(1);
#endif
    }

    /* if there are constraints on the critical lines positions in image keyword*/
    if (I.npcl != 0)
    {
#ifdef CHIRES
        fprintf(OUT, "chi critical ligne\n");
#endif
        double ampa, ampb, d;
        chil = 0.;
        for (i = 0; i < I.npcl; i++)
        {
            A.x = cl[i].C.x - cl[i].dl * cos(cl[i].phi);
            A.y = cl[i].C.y - cl[i].dl * sin(cl[i].phi);
            B.x = cl[i].C.x + cl[i].dl * cos(cl[i].phi);
            B.y = cl[i].C.y + cl[i].dl * sin(cl[i].phi);
            ampa = e_amp(&A, cl[i].dl0s, cl[i].dos, cl[i].z);
            ampb = e_amp(&B, cl[i].dl0s, cl[i].dos, cl[i].z);

            if (ampa*ampb < 0.)
                // look for the critical line position between A and B
                // the precision if fixed by PREC_ZERO
                D = e_zeroamp(A, B, cl[i].dl0s, cl[i].dos, cl[i].z);
            else
            {
                if ((B.x != A.x) && (ampb != ampa))
                    D.x = (ampb * A.x - ampa * B.x) / (ampb - ampa);
                else
                    D.x = A.x;
                if ((B.y != A.y) && (ampb != ampa))
                    D.y = (ampb * A.y - ampa * B.y) / (ampb - ampa);
                else
                    D.y = A.y;
            }

            d = dist(D, cl[i].C) / cl[i].dl;
            chil += d * d;

            // add parity chi2 (side A must be +-, side B must be ++)
            ampli = e_unmag(&A, cl[i].dl0s, cl[i].dos, cl[i].z);
            if( ampli.a < 0. ) chil += ampli.a * ampli.a;
            if( ampli.b > 0. ) chil += ampli.b * ampli.b;
            ampli = e_unmag(&B, cl[i].dl0s, cl[i].dos, cl[i].z);
            if( ampli.a < 0. ) chil += ampli.a * ampli.a;
            if( ampli.b < 0. ) chil += ampli.b * ampli.b;

#ifdef CHIRES
            fprintf(OUT, "%d  %.3lf %.3lf %.2lf\n", i, cl[i].C.x, cl[i].C.y, d*d);
        }

        fprintf(OUT, "chil %.2lf\n", chil);
#else
        }
#endif
    }
///////////////////////////////////////	
if (Dy.dynnumber == 1 || Dy.dynnumber == 2 || Dy.dynnumber == 3 || Dy.dynnumber == 4)  
{
	//extern struct pot  lens[NLMAX];
	double  vel_model;     //TV Oct2011
	double  mass_model;
	
	
	if (Dy.dynnumber == 1 || Dy.dynnumber == 2) 
		
	{
		
	    vel_model = lens[0].sigma/1.473972264;      //Transformation of the velocity in "our" velocity
		
	    chi_vel = (Dy.dynvel - vel_model) * (Dy.dynvel - vel_model) / Dy.dynevel / Dy.dynevel;
		
	    lhood0 +=  log( 2.*M_PI*(Dy.dynevel*Dy.dynevel) );    // NORMALIZATION
		
	}	
	
	
	if (Dy.dynnumber == 2) 
		
	{
		
		mass_model = mass2d_NFW(lens[0].sigma,Dy.refradius,lens[0].rckpc);  //Analitical mass
		
		chi_mass = (Dy.indmass - mass_model)*(Dy.indmass - mass_model) / Dy.indemass / Dy.indemass;
		
		lhood0 +=  log( 2.*M_PI*(Dy.indemass*Dy.indemass) );    // NORMALIZATION
		
	}
	
	if (Dy.dynnumber == 3) 
		
	{
		
		mass_model = mass3d_NFW(lens[0].sigma,Dy.refradius,lens[0].rckpc);  //Analitical mass
		
		chi_mass = (Dy.indmass - mass_model)*(Dy.indmass - mass_model) / Dy.indemass / Dy.indemass;
		
		lhood0 +=  log( 2.*M_PI*(Dy.indemass*Dy.indemass) );    // NORMALIZATION
		
	}
	
/*	if (Dy.dynnumber == 4)   //MAMPOST		
	{
      Full MAMPOSSt will be implemented soon
		
	}*/
#ifdef CHIRES
	fprintf(OUT,"chi_vel %7.4lf\n",chi_vel);
	fprintf(OUT,"chi_mass %7.4lf\n",chi_mass);
#endif
}
////////////////////////
    chi = chi_im + chis + chisrcfit + chil + chi_mult + chi_vel + chi_mass;  //TV //Min(chi_mult, 100000000.);
//    printf("chi_mult %lf chisrcfit %lf chil %lf\n",chi_mult,chisrcfit,chil);
#ifdef CHIRES
    fprintf(OUT, "\nchitot           %.2lf\n", chi);
    fprintf(OUT, "log(Likelihood)    %.2lf\n", -0.5*(chi + lhood0));
    fclose(OUT);
#endif

    /* NPRINTF(stderr,"chip=%.3lf",chip); */

    /* NPRINTF(stderr,f"chia=%.3lf  chip=%.3lf  chix=%.3lf  chiy=%.3lf\n",chia,chip,chix,chiy); */
#ifndef CHIRES
    *chi_ext = chi;
    *lhood0_ext = lhood0;
    return 0;    // OK
#endif
}

#ifndef CHIRES
// Return the chi2. Return -1 if error
double o_chi()
{
    int error;
    double chi2, lhood0;
    error = o_chi_lhood0(&chi2, &lhood0, NULL);
    if ( !error )
        return(chi2);
    else
        return(-1);
}

#ifdef __SEGER__O_CHI_TEST
static int o_lhood_count = 0;
struct timeval o_lhood_tv1;
#endif

double o_lhood(int *error)
{
    double chi2, lhood0;


    *error = o_chi_lhood0(&chi2, &lhood0, NULL);

#ifdef __SEGER__O_CHI_TEST
    o_lhood_count++;
    fprintf(stdout, "%d %g %g\n", o_lhood_count, chi2, lhood0);

    if (o_lhood_count == 1)
        gettimeofday(&o_lhood_tv1, NULL);

    if (o_lhood_count == 200)
    {
        struct timeval tv2;
        gettimeofday(&tv2, NULL);

        fprintf(stdout, "time %g\n", (double)(tv2.tv_sec - o_lhood_tv1.tv_sec) +
                1e-6 * (tv2.tv_usec - o_lhood_tv1.tv_usec));
        //just in case o_chires test
        o_chires("o_chires.txt");
        exit(0);
    }
#endif

    return(-0.5*(chi2 + lhood0));
}
#endif
/* Return the chi2 to a set of arclet according to the
 * I.stat chosen method.
 */
static void chi2_wl(double *chi2, double *lh0)
{
    extern struct g_mode    M;
    const extern struct g_image   I;
    extern struct g_source  S;
    extern long int narclet;
    extern struct galaxie   arclet[NAMAX];

    long int i;
    double ells, wht_ell, es1, es2; /*error and weight in ellipticity in the image plane*/
    double a, b, e, g, g1, g2, e1, e2;
    double chi2i;
    double x, qp, dp, tp, k;
    struct ellipse ell_sou, source, ampli; /* variables locales images */
    struct matrix grad2;

#ifdef CHIRES
    double ga; // shear,kappa,reduced shear estimator
    fprintf(OUT, "chis arclets\n");
    fprintf(OUT, "     N       ID      z        chi2      gamma         1-k          g            es1          es2\n");
#endif

    if ( I.stat == 1 || I.stat == 2 )
    {
        for (i = 0; i < narclet; i++)
        {
            /* optimisation du z */
            if (I.stat == 2)
            {
                printf("WARN: redshift optimisation with I.stat=2 not yet implemented.\n");
//                if (indexCmp(I.nza, arclet[i].n)) //TO CHECK (previous version : I.nza!=i+1)
//                {
//                    SC = arclet[i].C;
//                    ampli = e_unmag_gal(&arclet[i]);
//                    elix = fabs(arclet[i].eps *
//                                cos(2 * (arclet[i].E.theta - ampli.theta)));
//                    //arclet[i].dr =
//                    //    zero_t(arclet[i].dr, arclet[i].dr + .05, fmin_ell);
//                }
//                else
//                    // all arclets are optimised but arclet[i]
//                    arclet[i].dr = I.drarclet;
            }

            /* calcul de ts optimise */

            ampli = e_unmag_gal(&arclet[i]);
            isoima(&arclet[i].E, &ampli, &ell_sou);
            x = ell_sou.b / ell_sou.a;
            x = .5 * (1. / x - x);

            if (I.stat == 2)
            {
//                if (indexCmp(I.nza, arclet[i+1].n)) //TO CHECK (previous version : I.nza!=i+1)
//                    *chi2 += x * x;
//                else
//                {
//                    // less weight to strong lensing arcs
//                    *chi2 += 50.*x * x;
//                }
            }
            else
                *chi2 += x * x;

        };  /*end for (i=0;i<narclet;i++)*/
    } /* if ((I.stat==1)||(I.stat==2))*/

    else if (I.stat == 3) /* optimisation de l'orientation */
    {
        for (i = 0; i < narclet; i++)
        {
            ampli = e_unmag_gal(&arclet[i]);
            x = arclet[i].tau * sin( 2.*(arclet[i].E.theta - ampli.theta));
            *chi2 += x * x;
#ifdef CHIRES
            fprintf(OUT, "%ld %s %.2lf\n", i, arclet[i].n, x*x);
#endif

        }
    }
    else if (I.stat == 4) /* optimisation orientation + ellipticite*/
    {
        for (i = 0; i < narclet; i++)
        {
            ampli = e_unmag_gal(&arclet[i]);
            qp = fabs(ampli.b / ampli.a);
            dp = (1. + qp * qp) / 2. / qp;
            tp = fabs(1. - qp * qp) / 2. / qp;
            x = dp * arclet[i].tau * cos( 2.*(arclet[i].E.theta - ampli.theta))
                - tp * arclet[i].dis;
            *chi2 += x * x;
        }
        //chis /=narclet;
    }
    else if (I.stat == 5) /* optimisation orientation + ellipticite*/
    {
        for (i = 0; i < narclet; i++)
        {
            ampli = e_unmag_gal(&arclet[i]);
            isoima(&arclet[i].E, &ampli, &source);
            x = fabs( (source.a * source.a - source.b * source.b)
                      / 2. / source.a * source.b);
            *chi2 += x * x / 0.001;
        }
        //chis /=narclet;
    }
    else if (I.stat == 6) /* optimisation ellipticite a la Marshall*/
    {
        double lh0_add  = 0;
        double chi2_add = 0;
        int num_threads = 1;
#ifdef _OPENMP
        num_threads = omp_get_max_threads();
#endif
        if (num_threads > narclet / 100)
            num_threads = narclet / 100;

#ifdef CHIRES
        num_threads = 1;
#endif

#pragma omp parallel for if (num_threads > 1) schedule(guided, 100) \
            private(ampli, source, ells, es1, es2, chi2i) \
            reduction(+: lh0_add, chi2_add) num_threads(num_threads)
        for (i = 0; i < narclet; i++)
        {
            // Get amplification matrix:
            ampli = e_unmag_gal(&arclet[i]);

//           WARNING! if you uncomment something here put
//           "local variables" into private section of #pragma omp
//            a = ampli.a; b = ampli.b;
//            e = (a * a - b * b) / (a * a + b * b);
//            g1 = e * cos(2. * ampli.theta);
//            g2 = e * sin(2. * ampli.theta);
//
//            a = arclet[i].E.a; b = arclet[i].E.b;
//            e = (a * a - b * b) / (a * a + b * b);
//            e1 = e * cos(2. * arclet[i].E.theta);
//            e2 = e * sin(2. * arclet[i].E.theta);
//
//            es1 = e1 - g1;
//            es2 = e2 - g2;

            // Apply amplification matrix to image shape to get source shape:
            isoima(&arclet[i].E, &ampli, &source);
            // In weak lensing regime, have that e = e_s + g, or that e_s = e - g
            // i.e. in this limit source contains the residual that goes into chi-sq
            ells = (source.a * source.a - source.b * source.b) / (source.a * source.a + source.b * source.b);
            es1 = ells * cos(2 * source.theta);
            es2 = ells * sin(2 * source.theta);

            //WARNING! if you uncomment something here then put
            //"local variables" into private section of #pragma omp
            // Ellipticity compnents are the objects that are Gaussian-distributed:
            //   cos2th = cos(2.0*source.theta);
            //   sin2th = sin(2.0*source.theta);
            //   ells1 = ells*cos2th;
            //   ells2 = ells*sin2th;
            // Product of Gaussian distributions centred on zero:
            //
            // WARNING! if you uncoment it then put chis to reduction statement
            //   chis += ells1*ells1*wht_ell + ells2*ells2*wht_ell;

            chi2i = es1 * es1 / (I.sig2ell + arclet[i].var1 );
            chi2i += es2 * es2 / (I.sig2ell + arclet[i].var2 );

            chi2_add += chi2i;
            //if ( I.statmode == 1 || I.dsigell != -1 )
            //{
            lh0_add += log( 2 * M_PI * (I.sig2ell + arclet[i].var1) ) +
                       log( 2 * M_PI * (I.sig2ell + arclet[i].var2) );

            //}

            // TODO: test this on simulated data...
#ifdef CHIRES
#pragma omp critical
            {
                ga = (ampli.a - ampli.b) / 2.;       //   gamma_i
                k = (ampli.a + ampli.b) / 2.;      // divided by 1-k_i
                g = ga / k;
                NPRINTF(stdout, "INFO: compute chi2 for arclet %ld/%ld\r", i + 1, narclet);
                fprintf(OUT, "%6ld   %6s    %.3lf    %6.4lf    %6.3le    %6.3le    %6.3le    %6.3le    %6.3le\n", i, arclet[i].n, arclet[i].z, chi2i, ga, k, g, es1, es2 );
            }
#endif

        }

//       printf( "%lf\n", chi2_add);
        *chi2 += chi2_add;
        *lh0  += lh0_add;
    }
#ifdef CHIRES
    else if (I.stat == 7 || I.stat == 8) // optimisation ellipticite superior a la Marshall
#else
    else if (I.stat == 7) // optimisation ellipticite superior a la Marshall
#endif
    {
        // TODO include shape estimation error in quadrature?
//       wht_ell = 1. / I.sig2ell;

        // Compute the normalisation factor
//        *lh0 += narclet * log( 2 * M_PI * I.sig2ell );
        double lh0_add  = 0;
        double chi2_add = 0;
        int num_threads = 1;
#ifdef _OPENMP
        num_threads = omp_get_max_threads();
#endif
        if (num_threads > narclet / 100)
            num_threads = narclet / 100;

#ifdef CHIRES
        num_threads = 1;
#endif

#pragma omp parallel for if (num_threads > 1) schedule(guided, 100) \
            private(grad2, k, g1, g2, g, e, e1, e2, x, es1, es2, chi2i) \
            reduction(+: lh0_add, chi2_add) num_threads(num_threads)
        for (i = 0; i < narclet; i++)
        {
            grad2 = e_grad2_gal(&arclet[i], NULL);

            // Get amplification matrix:
            grad2.a /= arclet[i].dos;
            grad2.b /= arclet[i].dos;
            grad2.c /= arclet[i].dos;

//            a = ampli.a; b = ampli.b;
//          // k = 0.5 * (grad2.a + grad2.c) <--> k = 1 - 0.5 * (ampli.a + ampli.b)
//          // gam = 0.5 * (ampli.a - ampli.c) <--> gam = sqrt(g1^2 + g2^2)
//          // gam * cos(2*ampli.theta) <-> g1 = 0.5 * (grad2.c - grad2.a)
//          // gam * sin(2*ampli.theta) <-> g2 = - grad2.b
            k = 0.5 * (grad2.a + grad2.c);
            g1 = 0.5 * (grad2.c - grad2.a) / (1. - k); // reduced shear
            g2 = - grad2.b / (1. - k);
            g = g1 * g1 + g2 * g2;
//            e = (grad2.a - grad2.c) / (a + b);  // g = gamma / (1-kappa)
//            g1 = e * cos(2. * ampli.theta);
//            g2 = e * sin(2. * ampli.theta);

//            a = arclet[i].E.a; b = arclet[i].E.b;
//            e = (a - b) / (a + b);
//            e1 = e * cos(2. * arclet[i].E.theta);
//            e2 = e * sin(2. * arclet[i].E.theta);

            // a & b are gamma1 and gamma2 (calibrated ellipticities)
            // In COSMOS : e = (a*a-b*b)/(a*a+b*b) (cos(2theta) + i sin(2theta))
            // (See Bartelmann & Schneider 2001 : Eq. 4.6 pg 49)

            e = (arclet[i].E.a * arclet[i].E.a - arclet[i].E.b * arclet[i].E.b) / (arclet[i].E.a * arclet[i].E.a + arclet[i].E.b * arclet[i].E.b);
            e1 = e * cos(2. * arclet[i].E.theta);
            e2 = e * sin(2. * arclet[i].E.theta);
            g1 *= -1; g2 *= -1; /// to match Bartelmann definition in Eq 4.6
            x = 1. + g - 2. * (g1 * e1 + g2 * e2);
            es1 = (e1 - 2. * g1 + e1 * (g1 * g1 - g2 * g2) + 2. * g1 * g2 * e2) / x;
            es2 = (e2 - 2. * g2 - e2 * (g1 * g1 - g2 * g2) + 2. * g1 * g2 * e1) / x;

            // Apply amplification matrix to image shape to get source shape:
            //isoima(&arclet[i].E, &ampli, &source);
            // In weak lensing regime, have that e = e_s + g, or that e_s = e - g
            // i.e. in this limit source contains the residual that goes into chisq
            //ells = (source.a - source.b) / (source.a + source.b);
            //ells should be Rayleigh-distributed:
            //*chi2 += ells * ells * wht_ell; // - 2.0 * log(ells);

            chi2i = es1 * es1 / (I.sig2ell + arclet[i].var1 );
            chi2i += es2 * es2 / (I.sig2ell + arclet[i].var2 );

            chi2_add += chi2i;
            lh0_add += log( 2 * M_PI * (I.sig2ell + arclet[i].var1) ) +
                       log( 2 * M_PI * (I.sig2ell + arclet[i].var2) );

#ifdef CHIRES
#pragma omp critical
            {
            ampli = e_unmag_gal(&arclet[i]);
            ga = (ampli.a - ampli.b) / 2.;       //   gamma_i
            k = (ampli.a + ampli.b) / 2.;      // divided by 1-k_i
            g = ga / k;
            NPRINTF(stdout, "INFO: compute chi2 for arclet %ld/%ld\r", i + 1, narclet);
            fprintf(OUT, "%6ld   %6s    %.3lf    %6.4lf    %6.3le    %6.3le    %6.3le    %6.3le    %6.3le\n", i, arclet[i].n, arclet[i].z, chi2i, ga, k, g, es1, es2 );
            }
#endif

        }

        *chi2 += chi2_add;
        *lh0 += lh0_add;

    }
//    else
//    {
//        fprintf(stderr, "ERROR in o_chi/arclet: I.stat = %d\n", I.stat);
//        exit(-1);
//    }
#ifdef CHIRES
    NPRINTF( stdout, "\n" );
    fprintf(OUT, "chis               %.2lf\n", *chi2);
    fprintf(OUT, "log(Likelihood)    %.2lf\n", -0.5*(*chi2 + (*lh0)));
#endif
}

/* Return the chi2 for a single image.
 * Parameters :
 * - pima : pointer to an observed single image
 * - ps : pointer to the corresponding predicted source position
 */
static double chi2SglImage( struct galaxie *pima, struct point *ps )
{
    const extern struct g_image I;
    struct bitriplet Tsol[NIMAX];   // list of triangles containing predicted
    //arclets for a single obs image
    struct point Bs; /*barycenter of a familly of I.mult[i] sources*/
    double chi2, I2x, I2y, dx, dy;
    int j, nimages;

#ifdef CHIRES
    double rmsi_tot, rmsi, chi22;

    rmsi_tot = 0.;
#endif

    chi2 = 0;

    I2x = I2y = I.sig2pos[0][0];

    //just an over check
    if (pima->gsource != NULL && pima->grid_dr < 0)
    {
        fprintf(stderr, "You should initializa galaxie::gsource as NULL\n");
        exit(EXIT_FAILURE);
    }

    //if it is first run we should initialize pima->gsource
    if (pima->gsource == NULL)
    {
        pima->gsource = (struct point (*)[NGGMAX][NGGMAX])
                        malloc(sizeof(struct point[NGGMAX][NGGMAX]));
        const extern struct g_grille  G;
//  fprintf(stderr,"%d %d\n", NGGMAX, G.ngrid);
        pima->grid_dr = -1;
    }

    // Check if the grid is initialize to the image redshift
    // But we should reconstruct grid in any case
//    if ( pima->grid_dr != pima->dr )
    {
        pima->grid_dr = pima->dr;
        e_unlensgrid( *(pima->gsource), pima->grid_dr);
    }

    //raise chip to a high value if we have more than 1 image
    nimages = inverse( (const struct point (*)[NGGMAX]) * (pima->gsource),  ps, Tsol);

    if ( nimages > 1 )
    {
        if ( fabs(I.forme) == 10 )
            I2x = I2y = pima->E.a * pima->E.a;
        else if ( fabs(I.forme) == 11 )
        {
            I2x = pima->E.a * cos(pima->E.theta);
            I2y = pima->E.b * cos(pima->E.theta);
            I2x = I2x * I2x;
            I2y = I2y * I2y;
        }

        // compute chip as the total distance between the predicted images
        // and the observed image
        for ( j = 0; j < nimages; j++ )
        {
            Bs = barycentre(&Tsol[j].i);
            dx = Bs.x - pima->C.x;
            dy = Bs.y - pima->C.y;
            chi2 += dx * dx / I2x + dy * dy / I2y;
#ifdef CHIRES
            chi22 = dx * dx / I2x + dy * dy / I2y;
            rmsi = dx * dx + dy * dy;
            rmsi_tot += rmsi;
            fprintf(OUT, " 0 %6s %.3lf   %d   %7.2lf  %6.2lf  %6.2lf  %6.2lf  %.2lf  %.2lf    %d\n", pima->n, pima->z, 1, chi22, dx*dx / I2x, dy*dy / I2y,
                    0., 0., sqrt(rmsi), 0);
#endif
        }
    }

#ifdef CHIRES
    if ( nimages != 0 )
        rmsi_tot = sqrt(rmsi_tot / nimages);

    fprintf(OUT, " 0 %6s %.3lf   %d   %7.2lf  %6.2lf  %6.2lf  %6.2lf  %.2lf  %.2lf    %d\n", pima->n, pima->z, nimages, chi2, 0., 0., 0., 0., rmsi_tot, 0);
#endif

    return(chi2);
}

/*******************************************************
 * Optimization with the arclet shape in the SOURCE PLANE.
 *******************************************************/
static int chi2_src( double *chi2, double *lh0, double *np_b0 )
{
    const extern struct g_grille  G;
    const extern struct g_image   I;
    extern double chip, chia, chix, chiy;
    extern struct galaxie   multi[NFMAX][NIMAX];
    extern struct galaxie source[NFMAX];
    const extern int optim_z;
    extern struct matrix amplifi_mat[NFMAX][NIMAX], amplifi_matinv[NIMAX][NIMAX];
    extern struct sigposStr sigposAs;

    struct point Ps[NIMAX];
    struct point Bs; //barycenter of a familly of I.mult[i] sources
    struct galaxie *pima; // pointer to an arclet

    double MA, MB, MC, MD;
    double sigx2, sigy2, da, fluxS, sig2flux;// ,siga2;

    int ntot; // total number of multiple images
    double wtot; // sum of images weight (multi[i][j].var1, default:1)
    double Ix, Iy, I2x, I2y, Ixy; //sigma and sigma^2 in the source plane
    double w2x, w2y, wxy; //1/sigma^2 in the image & source plane
    double lhood0;  // local lh0
    double dx, dy; //error in position in the source plane
    double chipij;
    int    optEllip;  // test accelerator (boolean)
#define NDXDY  300
    double tmp[NDXDY];
    int    dcnt, ndxdy;   // image counter
    register int i, j, k;

#ifdef CHIRES
    struct point multib[NIMAX];
    int    nimages;  // number of valid images in multib
    extern int nwarn; /*counts each time a barycenter is not found in a
                             small source triangle. see e_im_prec() function*/
    char   fname[15];  // root name of a multiple images system

    double chipi, chixi, chiyi, chiai; // per family chi2
    double chixij, chiyij, chiaij; // per image chi2
    double rmss, rmsi;  // per family rms in source and image plane
    double rmssj, rmsij;    // per image "rms" in source and image plane
    double rmss_tot, rmsi_tot;  // image rms tot in source and image plane
    int    nwarni;  // +1 if error in predicted image

    fprintf(OUT, "chi multiples\n");
    fprintf(OUT, " N    ID    z   Narcs    chip    chix    chiy    chia   rmss     rmsi    dx      dy    nwarn\n");

    rmss_tot = rmsi_tot = 0;
    nwarn = 0;
#endif
    //
    double d[40];
//  double sigposMean;
    chip = chix = chiy = chia = 0.;
    lhood0 = 0.;
    ntot = wtot = 0.;

    // Just to prevent warning messages during compilation
    // I2x and I2y are properly computed later
    I2x = I2y = sigposAs.min * sigposAs.min;

    // Regularisation of sig2pos
//  sigposMean=0.; ntot=0;
//  for( i=0 ; i< I.n_mult; i++)
//      for( j = 0; j < I.mult[i]; j++)
//      {
//          sigposMean+=sqrt(sig2pos[i][j]);
//          ntot++;
//      }
//  sigposMean/=ntot;
//  for( i=0 ; i< I.n_mult; i++)
//      for( j = 0; j < I.mult[i]; j++)
//      {
//          dx=sqrt(sig2pos[i][j])-sigposMean;
//          chip+=dx*dx/I2x;
//      }
//  *lh0 += log( 2.*M_PI*I2x );
//
    // It's not very necessary for I.forme==5 but it can be worth that
    // in the rest of the code all the images have their amplification
    // computed at some point, even if it's not necessary today (EJ: 07/07).
//    if ( I.forme <= 7 || I.forme == 10 )
//        amplif(np_b0, amplifi);   //amplification of each image of all families
//        EJ(22-05-2011): amplifi is removed and a call to e_amp_gal is made in each function which requires it

    if ( I.forme == 8 || I.forme == 12 || I.forme == 13 )
        amplif_mat();

    if ( I.forme == 9 )
        amplif_matinv();

    // compute the total number of images, and initialise tmp[]
    if ( I.forme == 12 )
    {
        ndxdy = 0;
        dcnt = 0;
        for ( i = 0; i < I.n_mult && ndxdy < NDXDY; i++ )
            for ( j = 0; j < I.mult[i] && ndxdy < NDXDY; j++ )
            {
                tmp[ndxdy] = tmp[ndxdy+1] = 0.;
                ndxdy += 2;
            }

        if ( ndxdy >= NDXDY )
        {
            fprintf(stderr, "ERROR: Too many images. Increase NDXDY in o_chi.c\n");
            exit(0);
        }
    }
 
    // Define the number of threads we want to use
    int num_threads = 1;
#ifdef _OPENMP
    num_threads = omp_get_max_threads();
#endif

#ifdef CHIRES
    num_threads = 1;
#endif
    if ( I.forme == 12 )
        num_threads = 1;

    /* for each familly of arclets*/
#pragma omp parallel for if (num_threads > 1) \
            private(optEllip,fluxS,i,I2x,I2y,Ix,Iy,Ps,pima,Bs, \
            j,chipij,dx,dy,MA,MB,MC,MD, \
            w2x,w2y,wxy,k,tmp,d,sigx2,sigy2,da,sig2flux) \
            reduction(+: lhood0,chip,chix,chiy,chia,wtot,ntot) num_threads(num_threads)
    for (i = 0; i < I.n_mult; i++)
    {
#ifdef CHIRES
        chipi = chixi = chiyi = chiai = 0.;
        rmss = rmsi = 0.;
        nwarni = 0;
#endif
        //int n_famille = i;

        // Initialise shortcut tests
        // for optimization with the ellipticity of the arclets
        optEllip = 0;
        if ( ( I.forme == 1 || I.forme == 3 ) &&
             ! ( multi[i][0].E.a == multi[i][0].E.b && multi[i][0].E.theta == 0. )
           )
        {
            optEllip = 1;
            o_mag_m(I.mult[i], multi[i]);
        }

        /*optimization with the flux of the arclets*/
        if ( I.forme == 2 || I.forme == 13 )
            o_flux(I.mult[i], &fluxS, i, np_b0);

        if ( I.forme <= 3 || I.forme == 10 )
            // etendue method. Source plane error averaged
            // over all the images of the system.
            I2x = I2y = sig2posS(I.sig2pos[i][0], I.mult[i], i, np_b0);

        if ( I.forme == 10 )
            Ix = I2x / I.sig2pos[i][0]; // Ix: normalised etendue source plane

        /*optim_z can be not null only with the o_runz*() functions*/
        if (optim_z == 0)
            /*Compute the sources positions in Ps */
            o_dpl(I.mult[i], multi[i], Ps, np_b0);
        else
            /*Compute the source position for each arclet*/
        {
            for (j = 0; j < I.mult[i]; j++)
            {
                pima = &multi[i][j];
                Ps[j].x = pima->C.x - pima->Grad.x * multi[i][0].dr;
                Ps[j].y = pima->C.y - pima->Grad.y * multi[i][0].dr;
            }

            Ps[I.mult[i]] = Ps[0];

            /* NPRINTF(stderr," dr=%.3lf Ps=%.3lf\n",multi[i][0].dr,Ps[0].x); */
        }

        /*Find the barycenter position of the computed sources*/
        if (I.forme == 4 || I.forme == 6)
            Bs = weight_baryc(Ps, multi[i], I.mult[i], i);
        else
            Bs = bcentlist(Ps, I.mult[i]);  /* barycentre */

        if( G.nlens != G.nmsgrid )
        {
            Bs.x = source[i].C.x;
            Bs.y = source[i].C.y;
        }

#ifdef CHIRES
        /* Compute the rms in the image plane.
           multib[] contains the predicted image positions
           in source plane.
           multib[] has the same size as multi[]
           In case of lost images, the predicted image position is
           the observed image position.
         */
        int det_stop = 0;
        nimages = unlens_bc(Ps, Bs, multi[i], multib, I.mult[i], i, &det_stop);


#endif

        // ****************************************************
        // If we contraint with a SINGLE IMAGE...
        // Exception!!! -->constraint in the image plane
        //
        if ( I.mult[i] == 1 )
        {
            chip += chi2SglImage( &multi[i][0], &Ps[0] ) * multi[i][0].var1;
            wtot += multi[i][0].var1;
            ntot ++;
        }
        else
            /******************************************************
             * constraint with MULTIPLE IMAGES (SOURCE PLANE)
             * in each familly i, for each arclet j, compute chip
             **/
            for (j = 0; j < I.mult[i]; j++)
            {
                // Initialisation of per image variables
                chipij = 0.;
                ntot ++;
#ifdef CHIRES
                chiaij = chixij = chiyij = 0.;
#endif

                /* distance in arcsec in the source plane between barycenter
                 * and computed source position for each arclet*/
                dx = Bs.x - Ps[j].x;
                dy = Bs.y - Ps[j].y;

                /* Modify the error in position according to the shear and
                 * convergence parameters in the source plane.
                * (EJ 09/01/09: This is the proper way of computing the chi2,
                * \theta_pred - \theta_obs = M^-1 ( \hat{\beta} - \beta_obs )
                */
                if (I.forme == 8 || I.forme == 13)
                {
                    double cost, sint, a2, b2;
                    double dx_i, dy_i;
                    double wi2x, wi2y, wixy;
                    MA = amplifi_mat[i][j].a;
                    MB = amplifi_mat[i][j].b;
                    MC = amplifi_mat[i][j].c;
                    MD = amplifi_mat[i][j].d;
                    dx_i = dx * MA + dy * MD;
                    dy_i = dx * MB + dy * MC;

                    // covariance matrix in amplification axis
                    cost = cos(multi[i][j].E.theta);
                    sint = sin(multi[i][j].E.theta);
             //       cost = 1.;
             //       sint = 0.;
                    a2 = multi[i][j].E.a * multi[i][j].E.a;
                    b2 = multi[i][j].E.b * multi[i][j].E.b;
                    wi2x = cost * cost / a2 + sint * sint / b2;
                    wixy = cost * sint / a2 - cost * sint / b2;
                    wi2y = cost * cost / b2 + sint * sint / a2;
                    // Not used in chi2 because of outflow errors in the denominator chi2
                    w2x = wi2x * MA * MA + 2.*wixy * MA * MD + wi2y * MD * MD;
                    w2y = wi2y * MC * MC + 2.*wixy * MB * MC + wi2x * MB * MB;
                    wxy = wi2x * MA * MB + wixy * (MA * MC + MB * MD) + wi2y * MC * MD;

//                    chipij = dx*dx*wi2x*MA*MA + dx*dx*wi2y*MD*MD + 2.*dx*dx*wixy*MA*MD
//                           + dy*dy*wi2y*MC*MC + dy*dy*wi2x*MB*MB + 2.*dy*dy*wixy*MB*MC
//                       + 2.*(dx*dy*wi2x*MA*MB + dx*dy*wi2y*MC*MD + dx*dy*wixy*MA*MC + dx*dy*wixy*MB*MD);

//                    chipij = dx * dx * w2x + dy * dy * w2y + 2. * dx * dy * wxy;
                    chipij = dx_i * dx_i * wi2x + dy_i * dy_i * wi2y + 2. * dx_i * dy_i * wixy;
                    // update normalisation factor for image ij
//                    lhood0 += log( 4.*M_PI * M_PI / (w2x * w2y - wxy * wxy));
                    lhood0 += log( 4.*M_PI * M_PI / (wi2x * wi2y - wixy * wixy));

                    // Detect infinite amplification issues
                    //if ( chipij < 0. ) return(1); // error
                }
                else if ( I.forme == 12 )
                    // compute chi2 with covariance matrix
                {
                    double dx_i, dy_i;
                    MA = amplifi_mat[i][j].a;
                    MB = amplifi_mat[i][j].b;
                    MC = amplifi_mat[i][j].c;
                    MD = amplifi_mat[i][j].d;
                    dx_i = dx * MA + dy * MD;
                    dy_i = dx * MB + dy * MC;

                    for ( k = dcnt + 1; k < ndxdy; k++ )
                        tmp[k] += 2. * dx_i * I.weight[dcnt][k];

                    chipij = dx_i * dx_i * I.weight[dcnt][dcnt] + tmp[dcnt] * dx_i;
                    d[dcnt] = dx_i;
                    dcnt++;

                    for ( k = dcnt + 1; k < ndxdy; k++ )
                        tmp[k] += 2. * dy_i * I.weight[dcnt][k];

                    chipij += dy_i * dy_i * I.weight[dcnt][dcnt] + tmp[dcnt] * dy_i;
                    // lhood0 is updated once by detCov at the end
                    d[dcnt] = dy_i;
                    dcnt++;
                }
                else if( G.nlens != G.nmsgrid ) 
                {
                        Ix = multi[i][j].E.a * multi[i][j].mu;
                        Iy = multi[i][j].E.b * multi[i][j].mu;
                        I2x = Ix * Ix;
                        I2y = Iy * Iy;

                        chipij = dx * dx / I2x + dy * dy / I2y;
                        lhood0 += log( 4.*M_PI * M_PI * I2x * I2y);
                }
                else
                {
                    if ( I.forme == 4  || I.forme == 7 )
                        // sqrt(etendue) method. Image plane error = seeing
                        I2x = I2y = sig2posSj(I.sig2pos[i][j], multi[i], I.mult[i], j, i);
                    else if ( I.forme == 5 || I.forme == 6 )
                        // brute force method with 4 points
                        I2x = I2y = sig2posS4j(I.sig2pos[i][j], multi[i], Ps[j], j);
                    else if ( I.forme == 10 )
                        // etendue method. Different image plane errors
                        I2x = I2y = Ix * multi[i][j].E.a * multi[i][j].E.b;
                    else if ( I.forme == 11 )
                        // etendue method. Elliptical error in the source plane
                        // considering the image size as 1sigma error
                        sig2posSe(&multi[i][j], &I2x, &I2y);

                    w2x = 1. / I2x;
                    w2y = 1. / I2y;
                    wxy = 0.;

                    /*modify the sigma according to the shear and convergence
                     * parameters in the source plane*/
                    if (I.forme == 9)
                    {
                        MA = amplifi_matinv[i][j].a;
                        MB = amplifi_matinv[i][j].b;
                        MC = amplifi_matinv[i][j].c;
                        MD = amplifi_matinv[i][j].d;
                        Ix = (MA + MD) * sqrt(I.sig2pos[i][j]);
                        Iy = (MB + MC) * sqrt(I.sig2pos[i][j]);
                        I2x = Ix * Ix;
                        I2y = Iy * Iy;
                        w2x = 1. / I2x;
                        w2y = 1. / I2y;
                    }

                    chipij = dx * dx * w2x + dy * dy * w2y + 2. * dx * dy * wxy;
                    // update normalisation factor for image ij
                    lhood0 += log( 4.*M_PI * M_PI / (w2x * w2y - wxy * wxy));
                }

                // sum total chip
                chip += chipij * multi[i][j].var1;
                wtot += multi[i][j].var1;

#ifdef CHIRES
                // chi2 and rms for familly i
                chipi += chipij;
                rmssj = dx * dx + dy * dy;
                rmss += rmssj;
#endif

                //optimization with the ellipticity of the arclets...
                if ( optEllip )
                {
                    o_shape(j, multi[i], &dx, &sigx2, &dy, &sigy2, &da);
                    chix += dx * dx / sigx2;
                    chiy += dy * dy / sigy2;
#ifdef CHIRES
                    chixij = dx * dx / sigx2;
                    chixi += chixij;
                    chiyij = dy * dy / sigy2;
                    chiyi += chiyij;
#endif
                    // ...  and the flux
                    if (I.forme == 3)
                    {
                        chia += da * da / I.sig2amp;
#ifdef CHIRES
                        chiaij = da * da / I.sig2amp;
                        chiai += chiaij;
#endif
                    }
                }

                //optimization with the flux of the arclets only
                if ( (I.forme == 2 || I.forme == 13) && multi[i][j].mag != 0. )
                {
                    o_chi_flux(&multi[i][j], fluxS, &da, &sig2flux);
                    chia += da * da / sig2flux;
                    lhood0 += log( 2.*M_PI * sig2flux );
#ifdef CHIRES
                    chiaij = da * da / sig2flux;
                    chiai += chiaij;
#endif
                }

#ifdef CHIRES
                // compute the error in image plane from the multib[]
                // computed earlier
                dx = multib[j].x - multi[i][j].C.x;
                dy = multib[j].y - multi[i][j].C.y;
                rmsij = dx * dx + dy * dy;
                int warnij = 0;
                if ( rmsij > 0 )
                    rmsi += rmsij;
                else
                {
                    warnij = 1;
                    nwarni += 1;
                }

                // Summarize all these calculations in one print
                fprintf(OUT, "%2d %6s %.3lf   1   %7.2lf  %6.2lf  %6.2lf  %6.2lf  %6.3lf  %6.2lf  %6.2lf  %6.2lf  %d\n",
                        i, multi[i][j].n, multi[i][j].z, chipij, chixij, chiyij, chiaij, sqrt(rmssj), sqrt(rmsij), -dx, dy, warnij);
#endif

            } /*end for each image*/

#ifndef CHIRES
    } /*end for each familly i ( in o_chi() )*/
#else

        // Finalise the statistics
        if ( nimages - nwarni > 0 )
            rmsi /= nimages - nwarni;
        else
            rmsi = 0;

        rmss /= I.mult[i];
        rmss_tot += rmss;
        rmsi_tot += rmsi;
        nwarn += nwarni;

        strcpy(fname, multi[i][0].n);
        fname[strlen(multi[i][0].n) - 1] = 0;
        fprintf(OUT, "%2d %6s %.3lf   %d   %7.2lf  %6.2lf  %6.2lf  %6.2lf  %6.3lf  %6.2lf    N/A     N/A   %d\n",
                i, fname, multi[i][0].z, I.mult[i], chipi, chixi, chiyi, chiai, sqrt(rmss), sqrt(rmsi), nwarni);

    } /*end for each familly i ( in chires() ) */
#endif
 
    // normalize chip, but we still want the full chi2, not the reduced chi2
    chip = chip * ntot / wtot;

#ifdef CHIRES
    rmss_tot /= I.n_mult;
    rmsi_tot /= I.n_mult;
    fprintf(OUT, "chimul                %7.2lf  %6.2lf  %6.2lf  %6.2lf  %6.3lf  %6.2lf    N/A     N/A   %d\n", chip, chix, chiy, chia, sqrt(rmss_tot), sqrt(rmsi_tot), nwarn);
    fprintf(OUT, "log(Likelihood)       %7.2lf\n", -0.5*(chip + chia + chix + chiy + lhood0));
#endif

    // update the total statistics
    if ( I.forme == 12 )
        *lh0 += ndxdy * log(2. * M_PI) + log(I.detCov);
    else
        *lh0 += lhood0;

    *chi2 = chip + chia + chix + chiy;
    return(0); // OK

}

static int chi2_img( double *chi2, double *lh0 )
{
#ifndef CHIRES
    extern struct g_mode    M;
#endif
    const extern struct g_image   I;
    extern double chip, chia, chix, chiy;
    extern struct galaxie   multi[NFMAX][NIMAX];
    extern int nwarn; /*counts each time a barycenter is not found in a
                             small source triangle. see e_im_prec() function*/
    const extern int optim_z;

#ifdef CHIRES
    double rmss_tot, rmsi_tot;  // image rms tot in source and image plane
    int ntot; // total number of multiple images

    fprintf(OUT, "chi multiples\n");
    fprintf(OUT, " N    ID    z   Narcs    chip    chix    chiy    chia   rmss     rmsi    dx      dy    nwarn\n");

    ntot = rmss_tot = rmsi_tot = 0;
#endif
    chip = chia = chix = chiy = 0.;
    nwarn = 0;

    //I. We calculate numbers of all images in all families.
    //It is equal to the number of times we should run unlens_bc_single
    //or chi2SglImage

    int task_num = 0; //number of time we should call unlens_bc_single
    int task_idx;     //index of task
    int i, j;
    for (i = 0; i < I.n_mult; i++)
        task_num  += I.mult[i];

    struct point *vPs;    //array of source position for all images from all familes
    struct point *vmultib;//array of multib for all images from all families
    struct point *vBs; /*array of barycenters of a familly of I.mult[i] sources*/
    int* task_imap;   //task_imap[task_idx] = n_familly
    int* task_jmap;   //task_imap[task_idx] = j

    vPs     = (struct point*) malloc(task_num * sizeof(struct point));
    vmultib = (struct point*) malloc(task_num * sizeof(struct point));
    vBs     = (struct point*) malloc(I.n_mult * sizeof(struct point));
    task_imap = (int*) malloc(task_num * sizeof(int));
    task_jmap = (int*) malloc(task_num * sizeof(int));

#ifdef CHIRES
    int *task_ok = (int*) malloc(task_num * sizeof(int)); //find or not image
#endif

    //II. We create task_imap and task_jmap to be able to make single for-loop
    //over all images from all families
    task_idx = 0;
    for (i = 0; i < I.n_mult; i++)
    {
        for (j = 0 ; j < I.mult[i] ; j++, task_idx++)
        {
            task_imap[task_idx] = i;
            task_jmap[task_idx] = j;
        }
    }

    //III. We calculate vPs[task_idx] and vBs[i]
    task_idx = 0;
    for (i = 0; i < I.n_mult; i++)
    {
        struct point Ps[NIMAX + 1];
        /* optim_z can be different from 0 with the o_runz*() functions*/
        if (optim_z == 0)
            /* compute the sources positions in Ps*/
            o_dpl(I.mult[i], multi[i], Ps, NULL);
        else
        {
            for (j = 0; j < I.mult[i]; j++)
            {
                struct galaxie *pima = &multi[i][j];
                Ps[j].x = pima->C.x - pima->Grad.x * multi[i][0].dr;
                Ps[j].y = pima->C.y - pima->Grad.y * multi[i][0].dr;
            }
            Ps[I.mult[i]] = Ps[0];
        }

        /*vBs[i] contains the barycenter of the I.mult[i] sources positions Ps*/
        vBs[i] = bcentlist(Ps, I.mult[i]);
        for (j = 0 ; j < I.mult[i]; j++, task_idx++)
        {
            vPs[task_idx] = Ps[j];
        }
    }

    //IV. We run unlens_bc_single for all images from all families in
    //the single for cycle. For single image families we run chi2SglImage.
    //This part is the most time consuming part so we can parallelise it.

    int det_stop = 0; // if det_stop then we should return -1

    // Define the number of threads we want to use
    int num_threads = 1;
#ifdef _OPENMP
    num_threads = omp_get_max_threads();
#endif
    if (num_threads > task_num / 2)
        num_threads = task_num / 2;

#ifdef CHIRES
    num_threads = 1;
#endif

#pragma omp parallel for schedule(dynamic,1) private(i,j) num_threads(num_threads)
    for (task_idx = 0 ; task_idx < task_num ; task_idx++)
    {
        if (det_stop)
            continue;
        i = task_imap[task_idx]; //number of family
        j = task_jmap[task_idx]; //number of image in family i
        if ( I.mult[i] == 1 ) // SINGLE IMAGE system
        {
            double chip_tmp = chi2SglImage( &multi[i][0], &vPs[task_idx] );

#pragma omp atomic
            chip += chip_tmp;
        }
        else //MULTIPLE IMAGES system
        {
            /*return in vmultib[task_idx] an image for current arclet */
            int ok = unlens_bc_single(vPs[task_idx], vBs[i], &multi[i][j],
                                      &vmultib[task_idx], i);
            //if !ok then we cannot find this image

#ifdef CHIRES
            task_ok[task_idx] = ok;
#endif
#ifndef CHIRES
            if ( !ok && M.inverse > 2 )
            {
#pragma omp atomic
                det_stop++;     //we "return -1"
                //flush(det_stop) surprisingly reduce performance,
                //so we don't do it
            }
            if (det_stop)
                continue;
#endif
        }
    }

    //V. Now we can caluculate chi2 using constructed vmultib
    task_idx = 0;
    for (i = 0; (i < I.n_mult) && (!det_stop); i++)
    {
        if (I.mult[i] == 1) //SINGLE IMAGE system
        {
            //we already add results to chip
            task_idx++;
        }
        else  //MULTIPLE IMAGES systems
        {
            double I2x, I2y; //sigma and sigma^2
            double dLhood0;  // lh0 for a given familly of multiple images (I.forme=10 || sigposAs.bk !=0)
            double dx, dy; //error in position in the source plane
            int j;

#ifdef CHIRES
            int nimages = 0;
            double chipi, chixi, chiyi, chiai;
            double rmss, rmsi;   // image rms in source and image plane
            chipi = chixi = chiyi = chiai = 0.;
            rmss = rmsi = 0.;
#endif
            dLhood0 = 0.;

            // In all cases
            I2x = I2y = I.sig2pos[i][0];


            // For each image in familly i
            for (j = 0; j < I.mult[i];  j++, task_idx++ )
            {
                //current multib ---> vmultib[task_idx]
                //current Bs     ---> Bs[i]
                //current task_ok---> task_ok[task_idx]

                struct galaxie *pima = &multi[i][j];

                if ( I.forme == -1 )
                {
                    I2x = I2y = I.sig2pos[i][j];
                }
                else if ( I.forme == -10 )
                {
                    I2x = I2y = pima->E.a * pima->E.b;
                }
                else if ( I.forme == -11 )
                {
                    I2x = pima->E.a; //* cos(pima->E.theta);
                    I2y = pima->E.b; //* cos(pima->E.theta);
                    I2x = I2x * I2x;
                    I2y = I2y * I2y;
                }

                // update normalisation factor for image ij
                dLhood0 += log( 4.*M_PI * M_PI * (I2x * I2y) );

                dx = vmultib[task_idx].x - multi[i][j].C.x;
                dy = vmultib[task_idx].y - multi[i][j].C.y;

                chip += dx * dx / I2x + dy * dy / I2y;
#ifdef CHIRES
                if (task_ok[task_idx])
                    nimages++;
                double rmsij, rmssj, chipij;
                chipi += dx * dx / I2x + dy * dy / I2y;
                chipij = dx * dx / I2x + dy * dy / I2y;
                rmsi  += dx * dx + dy * dy;
                rmsij  = dx * dx + dy * dy;
                // source plane rmss
                double dxs, dys;
                dxs = vBs[i].x - vPs[task_idx].x;
                dys = vBs[i].y - vPs[task_idx].y;
                rmss += dxs * dxs + dys * dys;
                rmssj = dxs * dxs + dys * dys;

                //TODO: do we need warnj variable here?
                int warnj = 0;

                fprintf(OUT, "%2d %6s %.3lf   1   %7.2lf  %6.2lf  %6.2lf  %6.2lf  %6.3lf  %6.2lf  %6.2lf  %6.2lf  %d\n", i, multi[i][j].n, multi[i][j].z, chipij, 0., 0., 0., sqrt(rmssj), sqrt(rmsij), -dx, dy, warnj);
#endif
            } //end "for j" cycle

            // update the total statistics
            *lh0 += dLhood0;

            // TODO: Implement chix, chiy and chia calculation in image plane
            /*
             for (j=0;j<I.mult[i];j++)
             {
             da=diff_mag(multi[i][j],multib[j]);
             chia += da*da/I.sig2amp;
             };*/
#ifdef CHIRES
            rmss_tot += rmss;
            rmsi_tot += rmsi;
            ntot += nimages;
            nwarn += I.mult[i] - nimages;

            if ( nimages != 0 )
            {
                rmss = rmss / nimages;
                rmsi = rmsi / nimages;
            }

            fprintf(OUT, "%2d %6s %.3lf   %d   %7.2lf  %6.2lf  %6.2lf  %6.2lf  %6.3lf  %6.2lf    N/A     N/A   %d\n",
                    i, multi[i][0].n, multi[i][0].z, I.mult[i], chipi, chixi, chiyi, chiai, sqrt(rmss), sqrt(rmsi), I.mult[i] - nimages);
#endif
        } // end of case with I.mult[i] > 1
    } /*end for each familly*/

    free(vPs);
    free(vmultib);
    free(vBs);
    free(task_imap);
    free(task_jmap);
#ifdef CHIRES
    free(task_ok);
    // Total statistics
    rmss_tot = sqrt(rmss_tot / ntot);
    rmsi_tot = sqrt(rmsi_tot / ntot);
    fprintf(OUT, "chimul                %7.2lf  %6.2lf  %6.2lf  %6.2lf  %6.3lf  %6.2lf    N/A     N/A   %d\n", chip, chix, chiy, chia, rmss_tot, rmsi_tot, nwarn);
    fprintf(OUT, "log(Likelihood)       %7.2lf\n", -0.5*(chip + (*lh0)));
#endif
    if (det_stop)
        return -1;
    *chi2 = chip + chia + chix + chiy;
#ifdef DEBUG
    printf( "All images found!\n");
#endif
    return 0; // no error, no warning
}


/*****************************************************************
 * Send the points in srcfit[] to source plane and fit a line
 * to the source plane points.
 */
static void srcfitLinFit(double *chi2, double *lh0)
{
    const extern struct g_image I;
    extern struct galaxie *srcfit;

    struct point P[NSRCFIT];
    double wt[NSRCFIT], t;
    double ss, sx, sy, st2, a, b, tmp;
    int i;

    // Send points to source plane and compute source plane weights
    o_dpl(I.nsrcfit, srcfit, P, NULL);

    // Algorithm from IDL linfit.pro
    ss = sx = sy = 0.;
    for ( i = 0; i < I.nsrcfit; i++ )
    {
        wt[i] = srcfit[i].E.a * srcfit[i].E.b ; //* fabs(e_amp_gal(&srcfit[i]));
        wt[i] = 1. / wt[i];
        ss += wt[i] * wt[i];
        sx += P[i].x * wt[i] * wt[i];
        sy += P[i].y * wt[i] * wt[i];
    }

    b = 0.;
    st2 = 0.;
    for ( i = 0; i < I.nsrcfit; i++ )
    {
        t = (P[i].x - sx / ss) * wt[i];
        b += t * P[i].y * wt[i];
        st2 += t * t;
    }

    b /= st2;
    a = (sy - sx * b) / ss;

#ifdef CHIRES
    fprintf(OUT, "\nsource plane fitting y = %lf + %lf * x\n", a, b);
    fprintf(OUT, " N    ID      wt     chi2\n");
#endif

    // Compute chi2
    for ( i = 0; i < I.nsrcfit; i++ )
    {
        tmp = (P[i].y - b * P[i].x - a) * wt[i];
#ifdef CHIRES
        fprintf(OUT, "%2d  %4s    %.3lf    %.2lf\n", i, srcfit[i].n, wt[i], tmp*tmp ); // (1. + a * a));
#endif
        *chi2 += tmp * tmp ;   /// ( 1.  + a * a);
        *lh0 += log( 2.*M_PI / wt[i] / wt[i] );
    }

#ifdef CHIRES
    fprintf(OUT, "chisrcfit            %.2lf\n", *chi2);
    fprintf(OUT, "log(Likelihood)      %7.2lf\n", -0.5*(*chi2 + (*lh0)));
#endif

}
