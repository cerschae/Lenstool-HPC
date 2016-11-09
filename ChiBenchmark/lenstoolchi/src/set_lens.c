#include<stdio.h>
#include<stdlib.h>
#include<float.h>
#include<string.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_cblas.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"

/****************************************************************/
/*      nom:        set_lens_par            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

static void createInvMat();
//void rhos2b0(double *np_b0, double *array);

/* Initialize the lens parameters after they have been read in the
 * .par file and in the potfile.
 *
 * Called after reading the .par file in init_grille()
 * Called iteratively by runpot1() and runpot2().
 *
 * For this iterations, the number of initializations for each lens
 * has to be minimum.
 */
void    set_lens()
{
    extern struct g_mode    M;
    extern struct g_grille  G;
    extern struct g_source  S;
    extern struct g_pot     P[NPOTFILE];
    extern struct pot       lens[];
    extern struct pot       lmin[], lmax[], prec[];
    extern int block[][NPAMAX];
    extern struct g_cosmo   C;

    double GG = 10.867;

    long int i;

    double d1;

    // Set scaling relations
    for( i = 0; i < G.npot; i++ )
        if ( P[i].ftype != 0 )
            setScalingRelations(&P[i]);
    
    //************************************************************
    // SET THE CLUMPS PARAMETERS
    //************************************************************
    for ( i = 0 ; i < G.nlens ; i++ )
    {
        // Converting distance in kpc to arcsec.
        d1 = d0 / C.h * distcosmo1( lens[i].z );

        // Compute the DLS/DS ratio for all clumps.
        lens[i].dlsds = dratio( lens[i].z, S.zs );

        // Set rcore value in kpc or in arcsec.
        if ( lens[i].rckpc != 0. )
            lens[i].rc = lens[i].rckpc / d1;
        else
            lens[i].rckpc = lens[i].rc * d1;

        // Set rcut value in kpc or in arcsec.
        if ( lens[i].rcutkpc != DBL_MAX )
            lens[i].rcut = lens[i].rcutkpc / d1;
        else if ( lens[i].rcut != DBL_MAX )
            lens[i].rcutkpc = lens[i].rcut * d1;

        // Set the Mass to Light (mtol) property.
        /* the 1.5 factor come from the fact that we used 6pi\sigma^2 instead of
         4pi\sigma^2 - still don't understand why we need the 2.5 factor */
        if ( lens[i].lum != 0 )
            lens[i].mtol = 2.5 * 1.5 * PI / GG *
                           (lens[i].sigma / 1000) * (lens[i].sigma / 1000) *
                           lens[i].rcutkpc / lens[i].lum;

        // elliptical parameters epot and q
        if ( lens[i].type == 0 || lens[i].type == 2 )
        {
            lens[i].emass = 0.;
            lens[i].epot = 0.;
            lens[i].theta = 0.;
            lens[i].type++;
        }
        else if ( lens[i].type ==  8  ||
                  lens[i].type == -2  ||
                  ( lens[i].type > 80 && lens[i].type < 90 ) )
        {
            if ( lens[i].emass == 0. && lens[i].epot != 0. )
                // emass is (a2-b2)/(a2+b2)
                lens[i].emass = 2.*lens[i].epot / (1. + lens[i].epot * lens[i].epot);
            else if ( lens[i].emass == 0. && lens[i].epot == 0. )
                lens[i].epot = 0.00001;
            else
                // epot is (a-b)/(a+b)
                lens[i].epot = (1. - sqrt(1 - lens[i].emass * lens[i].emass)) / lens[i].emass;
        }
        else
        {
            if ( lens[i].emass == 0. && lens[i].epot != 0. )
                lens[i].emass = 3.*lens[i].epot;
            else
                lens[i].epot = lens[i].emass / 3.;

        }

        // Dynamical parameters
        if ( lens[i].type == 12 && lens[i].alpha != 0 )
            read_lenstable();
	if ( lens[i].type == 15) //einasto
		{read_table_einasto();
		}
        // TODO
        //if ( lens[i].type == 12 )
        //    lens[i].cr = rho_crit(lens[i].z);

        set_dynamics(i);
    }

    //************************************************************
    // GRID POTENTIALS : CREATE G.INVMAT & CONVERT RHOS TO B0
    //************************************************************
    if ( G.nmsgrid < G.nlens )
    {
        for ( i = G.nmsgrid; i < G.nlens; i++ )
        {
            if (lens[i].pmass != 0. && lens[i].sigma == 0.) 
                lens[i].b0 = lens[i].pmass;
        }

        // Create G.invmat
        //createInvMat();
    
//        int check = 0;
        // check if potentials are defined with rhos or v_disp. v_disp is dominant
//        for ( i = G.nmsgrid; i < G.nlens; i++ )
//            if (lens[i].pmass != 0. && lens[i].sigma == 0.)      check++;
    
        // if defined with rhos, then convert rhos to b0
//        if (check > 0)
//        {
//            double *atmp = (double *) malloc((G.nlens - G.nmsgrid) * sizeof(double));
//    
//            for (i = G.nmsgrid; i < G.nlens; i++)
//                atmp[i] = lens[i].pmass;
//    
//            setgrid_rhos2b0(atmp);
//    
//            for (i = G.nmsgrid; i < G.nlens; i++)
//                lens[i].b0 = atmp[i - G.nmsgrid];
//    
//            free(atmp);
//        }
    }

    //************************************************************
    // SET THE LIMITS PARAMETERS
    //************************************************************
    for ( i = 0 ; i < G.no_lens ; i++ )
    {
        lmin[i].type = lmax[i].type = prec[i].type = lens[i].type;

        // Converting distance in kpc to arcsec
        d1 = d0 / C.h * distcosmo1(lens[i].z);

        // set RCORE value, lmin, lmax and prec in kpc or in arcsec
        if ( lmin[i].rckpc != 0. )
            lmin[i].rc = lmin[i].rckpc / d1;
        else
            lmin[i].rckpc = lmin[i].rc * d1;

        if ( lmax[i].rckpc != 0. )
            lmax[i].rc = lmax[i].rckpc / d1;
        else
            lmax[i].rckpc = lmax[i].rc * d1;

        if ( prec[i].rckpc != 0. )
            prec[i].rc = prec[i].rckpc / d1;
        else
            prec[i].rckpc = prec[i].rc * d1;

        // set RCUT value, lmin, lmax and prec in kpc or in arcsec
        if ( lmin[i].rcutkpc != 0. )
            lmin[i].rcut = lmin[i].rcutkpc / d1;
        else
            lmin[i].rcutkpc = lmin[i].rcut * d1;

        if ( lmax[i].rcutkpc != 0. )
            lmax[i].rcut = lmax[i].rcutkpc / d1;
        else
            lmax[i].rcutkpc = lmax[i].rcut * d1;

        if ( prec[i].rcutkpc != 0. )
            prec[i].rcut = prec[i].rcutkpc / d1;
        else
            prec[i].rcutkpc = prec[i].rcut * d1;

        // Elliptical parameters
        if ( lens[i].type == 0 || lens[i].type == 2 || lens[i].type == 15/*einasto*/ )
        {
            // Circular SIS
            block[i][EPOT] = 0;
            block[i][EMASS] = 0;
        }

        // Dynamical parameters
        if (lens[i].type == 1 || lens[i].type == -1)
        {
            lmin[i].b0 = 4.*pia_c2 * lmin[i].sigma * lmin[i].sigma;

            if ( M.inverse >= 3 && block[i][B0] == 3 )
                lmax[i].sigma = 8.*pia_c2 * lmin[i].sigma * lmax[i].sigma;
            else
                lmax[i].b0 = 4.*pia_c2 * lmax[i].sigma * lmax[i].sigma;

            prec[i].b0 = 4.*pia_c2 * prec[i].sigma * prec[i].sigma;
        }
        else if (lens[i].type == 7)
        {
            lmin[i].b0 = 4.*RTA * GM_c2 *
                         lmin[i].masse / (D0 / C.h * distcosmo1(lens[i].z));
            lmax[i].b0 = 4.*RTA * GM_c2 * lmax[i].masse / (D0 / C.h * distcosmo1(lens[i].z));
            prec[i].b0 = 4.*RTA * GM_c2 *
                         prec[i].masse / (D0 / C.h * distcosmo1(lens[i].z));
        }
        else if (lens[i].type == 9)
        {
            lmin[i].b0 = lmin[i].pmass / cH2piG / C.h * distcosmo1(lens[i].z);
            lmax[i].b0 = lmax[i].pmass / cH2piG / C.h * distcosmo1(lens[i].z);
            prec[i].b0 = prec[i].pmass / cH2piG / C.h * distcosmo1(lens[i].z);
        }
        else if (lens[i].type == 5)
        {
            lmin[i].b0 = 18.*RTA * lmin[i].sigma * lmin[i].sigma / vol / vol;
            lmax[i].b0 = 18.*RTA * lmax[i].sigma * lmax[i].sigma / vol / vol;
            prec[i].b0 = 18.*RTA * prec[i].sigma * prec[i].sigma / vol / vol;
        }
        else if ( lens[i].type == 13 ) //Sersic : Keep the limits as they are
        {
            d1 = cH0_4piG * C.h / distcosmo1(lens[i].z);    // in 10^12 Msol/kpc^2
            lmin[i].b0 = lmin[i].sigma * 1e-12 / d1;
            lmax[i].b0 = lmax[i].sigma * 1e-12 / d1;
            prec[i].b0 = prec[i].sigma * 1e-12 / d1;
        }
	else if( lens[i].type==15) //einasto
	{
		d1=cH0_4piG*C.h/distcosmo1(lens[i].z);
		lmin[i].b0=lmin[i].pmass*1e-12/d1;
		lmax[i].b0=lmax[i].pmass*1e-12/d1;
		prec[i].b0=prec[i].pmass*1e-12/d1;
	}
        else
        {
            lmin[i].b0 = 6.*pia_c2 * lmin[i].sigma * lmin[i].sigma;
            lmax[i].b0 = 6.*pia_c2 * lmax[i].sigma * lmax[i].sigma;
            prec[i].b0 = 6.*pia_c2 * prec[i].sigma * prec[i].sigma;
        }
    }

    // Assign limits to msgrid potentials
//    if ( G.nmsgrid < G.nlens )
//    {
//        double *atmp = (double *) malloc((G.nlens - G.nmsgrid) * sizeof(double));
//        double *ab0 = (double *) malloc((G.nlens - G.nmsgrid) * sizeof(double));
//        for( i = G.nmsgrid; i < G.nlens; i++ )
//            atmp[i - G.nmsgrid] = lmax[i].pmass;
//        
//        rhos2b0(ab0, atmp);
//    
//        for( i = G.nmsgrid; i < G.nlens; i++ )
//            lmax[i].b0 = ab0[i - G.nmsgrid];
//
//        free(atmp);
//        free(ab0);
//    }

}

void set_dynamics(long int i)
{
    extern struct g_mode   M;
    extern struct g_grille G;
    extern struct pot      lens[];
    extern struct g_cosmo  C;
    
    register int ii, jj;
    double test;

    extern double *v_xx;
    extern double *v_yy;
    extern double **map_p;
    extern double **tmp_p;
    extern double **map_axx;
    extern double **map_ayy;

    double **tmpf, lc200, lm200;
    char    mode[20], nature[20], type[20], comment[1024];
    struct pot *ilens;

    ilens = &lens[i];

    switch (ilens->type)
    {
        case(1): // SIS
            ilens->b0 = 4.*pia_c2 * ilens->sigma * ilens->sigma;
            ilens->ct = ilens->b0 * ilens->dlsds;
            ilens->cr = 0.;
            break;
        case(12): // NFW
            // rescale from c,mvir,rhos --> sigmas (km/s), rs (arcsec)
            if ( ilens->beta != 0 && ilens->masse != 0 )
                // NFW defined by concentration and m200
            {
                e_nfw_cm200_sigrs( ilens->beta, ilens->masse,
                                   &ilens->sigma, &ilens->rckpc, ilens->z );
                ilens->rc = ilens->rckpc / (d0 / C.h * distcosmo1(ilens->z));
            }
            else if ( ilens->beta != 0 && ilens->rcut != DBL_MAX )
                // NFW defined by concentration and r200
            {
                ilens->rcutkpc = ilens->rcut * (d0 / C.h * distcosmo1(ilens->z));
                e_nfw_cr200_sigrs( ilens->beta, ilens->rcutkpc,
                                   &ilens->sigma, & ilens->rckpc, ilens->z );
                ilens->rc = ilens->rckpc / (d0 / C.h * distcosmo1(ilens->z));
            }
            else if ( ilens->beta != 0 && ilens->rc != 0 )
                // NFW defined by concentration and scale_radius
            {
                ilens->rckpc = ilens->rc * (d0 / C.h * distcosmo1(ilens->z));
                e_nfw_crs2sig( ilens->beta, ilens->rckpc, &ilens->sigma, ilens->z );
            }
            else if ( ilens->rc != 0 && ilens->masse != 0 )
                // NFW defined by scale_radius and m200
            {
                ilens->rckpc = ilens->rc * (d0 / C.h * distcosmo1(ilens->z));
                e_nfw_rsm200_sigrs( ilens->rckpc, ilens->masse, &ilens->sigma, ilens->z);
            }
            else if ( ilens->rc != 0 && ilens->rcut != DBL_MAX )
                // NFW defined by scale_radius and r200
            {
                ilens->rckpc = ilens->rc * (d0 / C.h * distcosmo1(ilens->z));
                ilens->rcutkpc = ilens->rcut * (d0 / C.h * distcosmo1(ilens->z));
                e_nfw_rsr200_sigrs( ilens->rckpc, ilens->rcutkpc, &ilens->sigma, ilens->z);
            }
            else if ( ilens->masse != 0 && ilens->beta == 0 )
                // NFW scaling relation from Maccio 2008
            {
                lm200 = log( ilens->masse ) / LOG10;
                lc200 = -0.098  * (lm200 - 12 ) + 0.830;
                ilens->beta = exp( lc200 * LOG10 );
                e_nfw_cm200_sigrs( ilens->beta, ilens->masse, &ilens->sigma, &ilens->rckpc, ilens->z );
                ilens->rc = ilens->rckpc / (d0 / C.h * distcosmo1(ilens->z));
            }
//            else if ( ilens->sigma == 0. )
//            {
//                fprintf( stderr, "ERROR: NFW potential with ID %s, badly defined\n", ilens->n );
//                exit(-1);
//            }

            ilens->b0 = 6.*pia_c2 * ilens->sigma * ilens->sigma;

            break;
		case(121): // Triaxial NFW
			// rescale from c,mvir,rhos --> sigmas (km/s), rs (arcsec)
			if( ilens->beta != 0 && ilens->masse != 0 )
			{
                fprintf(stderr, "ERROR: c and m200 to sigma, rs not yet implemented\n");
                exit(1);
			} 
			else if( ilens->beta != 0 && ilens->rcut != DBL_MAX )
			{
                fprintf(stderr, "ERROR: c and r200 to sigma, rs not yet implemented\n");
                exit(1);
			}
			else if( ilens->beta != 0 && ilens->rc != 0 )
			{
                double c2D = 0.;
                e_nfw_c3D2c2D(ilens->beta, &c2D);
				ilens->rckpc = ilens->rc * (d0/C.h*distcosmo1(ilens->z));
				e_nfw_crs2sig( c2D, ilens->rckpc, &ilens->sigma, ilens->z );
			}
			else if( ilens->rc != 0 && ilens->masse != 0 )
			{
                fprintf(stderr, "ERROR: rs and m200 to sigma, rs not yet implemented\n");
                exit(1);
			}
			else if( ilens->rc != 0 && ilens->rcut != DBL_MAX )
			{
                fprintf(stderr, "ERROR: rs and r200 to sigma, rs not yet implemented\n");
                exit(1);
			}

			ilens->b0=6.*pia_c2*ilens->sigma*ilens->sigma;
			
			break;
        case(13): // Sersic
            // here b0 is kappa(Re) = Sigma_e/Sigma_crit // cH0_4piG in 10^12 Msol/kpc^2
            ilens->b0 = ilens->sigma * 1e-12 / (cH0_4piG * C.h / distcosmo1(ilens->z));
            break;
	case(15) : //Einasto
	    //ici b0 représente Kappa_critique=sigma(x)/sigma_critique avec CH0_4piG en 10^12 Msol/kpac^2
	    ilens->b0=ilens->pmass*1e-12/(cH0_4piG*C.h/distcosmo1(ilens->z));
	    break;
        case(16): // Hernquist model
            // Default expression for testing. Then, need to match MOKA definition
            ilens->b0 = 6 * pia_c2 * ilens->sigma * ilens->sigma;
	    break;
        case(14): // External shear
            break;
        case(-1): //SIE
            ilens->b0 = 4.*pia_c2 * ilens->sigma * ilens->sigma;
            ilens->ct = ilens->b0 * ilens->dlsds;
            ilens->cr = 0.;
            break;
        case(-2):
//          NPRINTF(stderr,"Clump %d: True Elliptical BBS model\n",i);
//          NPRINTF(OUT,"-------- Clump %d: True Elliptical BBS model \n",i);
            ilens->b0 = 4.*pia_c2 * ilens->sigma * ilens->sigma;
            ilens->ct = ilens->b0 * ilens->dlsds;
            ilens->cr = 0.;
            break;
        case(7):
//          NPRINTF(stderr,"Clump %d: Point Masse \n",i);
//          NPRINTF(OUT,"Clump %d: Point Masse \n",i);
            ilens->b0 = 4.*RTA * GM_c2 *
                        ilens->masse / (D0 / C.h * distcosmo1(ilens->z));
            ilens->ct = sqrt(ilens->b0 * ilens->dlsds);
            ilens->cr = 0.;
            break;
        case(9):
//          NPRINTF(stderr,"Clump %d: Plan Masse \n",i);
//          NPRINTF(OUT,"Clump %d: Plan Masse \n",i);
            ilens->b0 = ilens->pmass / cH2piG / C.h * distcosmo1(ilens->z);
            ilens->ct = 0.;
            ilens->cr = 0.;
            break;
        case(5):
//          NPRINTF(stderr,"Clump %d: Hubble Modified Law \n",i);
//          NPRINTF(OUT,"Clump %d: Hubble Modified Law \n",i);
            ilens->b0 = 18.*RTA * ilens->sigma * ilens->sigma / vol / vol;
            /*
            ilens->phi0=ilens->b*ilens->rc;
            phi0=ilens->phi0; coeur=ilens->rc;
            ilens->ct=iter(ilens->phi0,ilens->rc,ilens->rc);
            ilens->cr=zero(0.,ilens->ct,fz_mhl_cr);
            */
            break;
        case(8):
//          NPRINTF(stderr,"Clump %d: PIEMD Kovner\n",i);
//          NPRINTF(OUT,"Clump %d: PIEMD Kovner\n",i);
            ilens->b0 = 6.*pia_c2 * ilens->sigma * ilens->sigma;
            break;
        case(81):
            /*          NPRINTF(stderr,"Clump %d: trunc. PIEMD Kovner:",i);
                        NPRINTF(OUT,"Clump %d: PIEMD Kovner, truncated\n",i);
                        NPRINTF(OUT," Total mass: %lf(cor) %lf (10^12 M_sol)\n",
                            4*M_PI/3*M_PI/GG*(ilens->sigma/1000)*(ilens->sigma/1000)*
                            ilens->rcut*(d0/C.h*distcosmo1(ilens->z)),
                            1.5*M_PI/GG*(ilens->sigma/1000)*(ilens->sigma/1000)*
                        ilens->rcut*(d0/C.h*distcosmo1(ilens->z)) );
                        NPRINTF(OUT," rcut:%.2lf(\") %.2lf(kpc)\n",
                            ilens->rcut,ilens->rcut*(d0/C.h*distcosmo1(ilens->z)) );
            */
            ilens->b0 = 6.*pia_c2 * ilens->sigma * ilens->sigma;
            break;

        case(82):
            /*          NPRINTF(stderr,"Clump %d: PIEMD Kovner, shallow center\n",i);
                        NPRINTF(OUT,"Clump %d: PIEMD Kovner, shallow center\n",i);
                        NPRINTF(OUT," Steep radius:%.2lf(\") %.2lf(kpc)\n",
                            ilens->rcut,ilens->rcut*(d0/C.h*distcosmo1(ilens->z)) );
            */
            ilens->b0 = 6.*pia_c2 * ilens->sigma * ilens->sigma;

            break;
        case(83):
            /*          NPRINTF(stderr,"Clump %d: EMD Kovner, 3/2\n",i);
                        NPRINTF(OUT,"Clump %d: PIEMD Kovner, shallow center\n",i);
                        NPRINTF(OUT," Steep radius:%.2lf(\") %.2lf(kpc)\n",
                            ilens->rcut,ilens->rcut*(d0/C.h*distcosmo1(ilens->z)) );
            */
            ilens->b0 = 6.*pia_c2 * ilens->sigma * ilens->sigma;
            break;
        case(84):
            /*          NPRINTF(stderr,"Clump %d: EMD Kovner, 0.5a-0.5s+1.5s\n",i);
                        NPRINTF(OUT,"Clump %d: PIEMD Kovner, 0.5a-0.5s+1.5s\n",i);
                        NPRINTF(OUT," Steep radius:%.2lf(\") %.2lf(kpc)\n",
                            ilens->rcut,ilens->rcut*(d0/C.h*distcosmo1(ilens->z)) );
            */
            ilens->b0 = 6.*pia_c2 * ilens->sigma * ilens->sigma;
            break;
        case(85):
            /*          NPRINTF(stderr,"Clump %d: EMD Kovner, 1\n",i);
                        NPRINTF(OUT,"Clump %d: PIEMD Kovner, 1a\n",i);
                        NPRINTF(OUT," Steep radius:%.2lf(\") %.2lf(kpc)\n",
                            ilens->rcut,ilens->rcut*(d0/C.h*distcosmo1(ilens->z)) );
            */
            ilens->b0 = 6.*pia_c2 * ilens->sigma * ilens->sigma;
            break;
        case(86):
            /*          NPRINTF(stderr,"Clump %d: EMD Kovner, 1a-1s\n",i);
                        NPRINTF(OUT,"Clump %d: PIEMD Kovner, 1a-1s\n",i);
                        NPRINTF(OUT," Steep radius:%.2lf(\") %.2lf(kpc)\n",
                            ilens->rcut,ilens->rcut*(d0/C.h*distcosmo1(ilens->z)) );
            */
            ilens->b0 = 6.*pia_c2 * ilens->sigma * ilens->sigma;
            break;
        case(87):
            /*          NPRINTF(stderr,"Clump %d: EMD Kovner, 1a-1s+0.5a-0.5s\n",i);
                        NPRINTF(OUT,"Clump %d: PIEMD Kovner, 1a-1s+0.5a-0.5s\n",i);
                        NPRINTF(OUT," Steep radius:%.2lf(\") %.2lf(kpc)\n",
                            ilens->rcut,ilens->rcut*(d0/C.h*distcosmo1(ilens->z)) );
            */
            ilens->b0 = 6.*pia_c2 * ilens->sigma * ilens->sigma;
            break;
        case(88):
            /*          NPRINTF(stderr,"Clump %d: EMD Kovner, 1a-1s+1.5s\n",i);
                        NPRINTF(OUT,"Clump %d: PIEMD Kovner, 1a-1s+1.5s\n",i);
                        NPRINTF(OUT," Steep radius:%.2lf(\") %.2lf(kpc)\n",
                            ilens->rcut,ilens->rcut*(d0/C.h*distcosmo1(ilens->z)) );
            */
            ilens->b0 = 6.*pia_c2 * ilens->sigma * ilens->sigma;
            break;
        case(89):
            /*          NPRINTF(stderr,"Clump %d: EMD Kovner, 1a-1s+0.5a-0.5s\n",i);
                        NPRINTF(OUT,"Clump %d: PIEMD Kovner, 1a-1s+0.5a-0.5s\n",i);
                        NPRINTF(OUT," Steep radius:%.2lf(\") %.2lf(kpc)\n",
                            ilens->rc*ilens->beta,ilens->rc*ilens->beta*(d0/C.h*distcosmo1(ilens->z)) );
            */
            ilens->b0 = 6.*pia_c2 * ilens->sigma * ilens->sigma;
            break;
        case(10):
//          NPRINTF(stderr,"Clump %d: Spline Potential\n",i);
//          NPRINTF(OUT,"Clump %d: Spline Potential\n",i);

            /* get the potential map */

            tmpf = (double **)rdf_ipx(G.splinefile, &G.ny, &G.nx,
                                      type, mode, nature, comment,
                                      &G.xmin, &G.xmax, &G.ymin, &G.ymax);

            map_p = (double **) alloc_square_double(G.nx, G.ny);
            tmp_p = (double **) alloc_square_double(G.nx, G.ny);

            for (ii = 0; ii < G.nx; ii++)
                for (jj = 0; jj < G.ny; jj++)
                    map_p[ii][jj] = tmpf[jj][ii];

            free_square_double(tmpf, G.nx);

            G.dx = (G.xmax - G.xmin) / (G.nx - 1);
            G.dy = (G.ymax - G.ymin) / (G.ny - 1);

            NPRINTF(stderr, "COMP: vecteurs x & y\n");
            v_xx = (double *) alloc_vector_double(G.nx);
            v_yy = (double *) alloc_vector_double(G.ny);
            for (ii = 0; ii < G.nx; ii++)
                v_xx[ii] = G.xmin + ii * G.dx;

            for (ii = 0; ii < G.ny; ii++)
                v_yy[ii] = G.ymin + ii * G.dy;

            /* compute the derivatives map */

            map_axx = (double **) alloc_square_double(G.nx, G.ny);
            map_ayy = (double **) alloc_square_double(G.nx, G.ny);

            sp_set(map_p, G.nx, G.ny, map_axx, map_ayy);

            NPRINTF(stderr, "Clump %ld: done\n", i);
            break;
        default:
//          NPRINTF(stderr,"Clump %d: Pseudo-Elliptical Potential with Core Radius\n",i);
//          NPRINTF(OUT,"Clump %d: Pseudo-Elliptical Potential with Core Radius\n",i);
            ilens->b0 = 6.*pia_c2 * ilens->sigma * ilens->sigma;
            test = ilens->dlsds * ilens->dlsds * ilens->b0 * ilens->b0 - ilens->rc * ilens->rc;

            if (test > 0.)
            {
                ilens->ct = sqrt(test);
                ilens->cr = sqrt(pow(ilens->b0 * ilens->dlsds, 2. / 3.) *
                                 pow(ilens->rc, 4. / 3.) - ilens->rc * ilens->rc);
            }
            else
                ilens->ct = ilens->cr = 0.;

            if (ilens->type > 20)
                updatecut(i);

            break;
    }
}



/* Convert the list of rhos values contained in np_b0 to b0 using matrix G.invmat
 * array must be the size G.nlens - G.nmsgrid. No check is performed.
 *
 * Global variables used: G.nmsgrid, G.nlens, lens[], np_b0
 */
void rhos2b0()
{
    extern struct g_grille  G;
    extern struct pot       lens[];
    extern double *np_b0;  // contains rhos values
    long int  i, j;
    double tmp;

    for( i = 0; i < G.nlens - G.nmsgrid; i++ )
    {
        tmp = 0.;  // reinitialization
        for( j  = 0; j < G.nlens - G.nmsgrid; j++ )
            tmp += G.invmat[i][j] * np_b0[j];

        lens[G.nmsgrid + i].b0 = tmp;
    }
    
}


//    gsl_vector *vsigma, *vb0;
//
//    // Allocate memory for vectors
//    vb0 = gsl_vector_alloc(G.nlens - G.nmsgrid);
//    vsigma = gsl_vector_alloc(G.nlens - G.nmsgrid);
//
//    // Inverse rhos to b0 with matrix G.invmat
//    memcpy(vsigma->data, array, (G.nlens - G.nmsgrid) * sizeof(double));
//    gsl_blas_dgemv(CblasNoTrans, 1., G.invmat, vsigma, 0., vb0);
//    memcpy( array, vb0->data, (G.nlens - G.nmsgrid) * sizeof(double));
//
//    // Clean vectors
//    gsl_vector_free(vsigma);
//    gsl_vector_free(vb0);
//  }    

/* Convert an array of rhos to b0 and assign the values to lens[] list
 * Argument array is modified into a list of b0
 */
//void setgrid_pmass2b0()
//{
//    extern struct g_grille G;
//    extern struct pot lens[];
//    extern double *np_b0;
//    long int ilens;
//
//    // Update lens[] with input rhos
//    for ( ilens = G.nmsgrid ; ilens < G.nlens ; ilens++ )
//        lens[ilens].pmass = array[ilens - G.nmsgrid];
//    
//    // Convert array of rhos to array of b0
//    rhos2b0(np_b0, array);
//
//    // Store data in lens[] array
//    for ( ilens = G.nmsgrid ; ilens < G.nlens ; ilens++ )
//        lens[ilens].b0 = np_b0[ilens - G.nmsgrid];
//}
//
/* Create the inverse matrix used to convert surface 
 * density Sigma vector to b0 vector
 *     [b0] = M^{-1} . [Sigma]
 *
 * This matrix is used in readBayesModel.c:setBayesModel()
 */
static void createInvMat()
{
    extern struct g_grille G;
    extern struct g_cosmo  C;
    extern struct pot      lens[];
    extern double *np_b0;

    struct matrix grad2;
    struct ellipse ampli;
    double m_ij, dl;
    int i, j, signum;
    gsl_matrix *mat, *invmat;
    gsl_permutation *p;
 
    // Create a global array for the lens.b0
    np_b0 = (double *) calloc(G.nlens - G.nmsgrid, sizeof(double));

    // set all lens[*].b0 to 1
    for ( i = G.nmsgrid; i < G.nlens; i++ )
    {
        np_b0[i - G.nmsgrid] = lens[i].b0;
        lens[i].b0 = 1.;
    }
    
    // Create the inverse matrix used to convert surface density Sigma vector to b0 vector
    invmat = gsl_matrix_alloc(G.nlens - G.nmsgrid, G.nlens - G.nmsgrid);
    // allocate matrix M and permutation vector p for LU decomposition
    mat = gsl_matrix_alloc(G.nlens - G.nmsgrid, G.nlens - G.nmsgrid);
    p = gsl_permutation_calloc(G.nlens - G.nmsgrid);

    dl = distcosmo1(lens[0].z);
    for( i = G.nmsgrid; i < G.nlens; i++ )
        for( j = G.nmsgrid; j < G.nlens; j++ )
        {
            grad2 = e_grad2_pot(&lens[i].C, j);
            ampli = formeli(1. - grad2.a, -grad2.b, 1. - grad2.c);
            m_ij = MCRIT12 / C.h * dl * (1. - (ampli.a + ampli.b) / 2.);  // 1e12 M/arcsec2
            gsl_matrix_set(mat, i - G.nmsgrid, j - G.nmsgrid, m_ij);
        }

    // Compute inverse of mat
    gsl_linalg_LU_decomp(mat, p, &signum);
    if( gsl_linalg_LU_invert(mat, p, invmat) )
    {
        fprintf(stderr, "ERROR: Singular matrix inversion in %s:%d\n", __FILE__, __LINE__);
        exit(1);
    }

    // Convert invmat to double[][]
    G.invmat = alloc_square_double(G.nlens - G.nmsgrid, G.nlens - G.nmsgrid);
    for( i = 0; i < G.nlens - G.nmsgrid; i++ )
        for( j = 0; j < G.nlens - G.nmsgrid; j++ )
            G.invmat[i][j] = gsl_matrix_get(invmat, i, j);

    // restore lens[i].b0 values
    for ( i = G.nmsgrid; i < G.nlens; i++ )
        lens[i].b0 = np_b0[i-G.nmsgrid];

    // Free matrix
    gsl_matrix_free(invmat);
    gsl_matrix_free(mat);
    gsl_permutation_free(p);
}



