#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>
#include<errors.h>
#include<gsl/gsl_rng.h>

#ifdef _OPENMP
#include "omp.h"
#endif

static void readWeightMat(double *** weight, double *detCov);

struct shear   shm[NAMAX];
struct galaxie multi[NFMAX][NIMAX];
struct z_lim   zlim[NFMAX];
struct z_lim   zalim;
struct galaxie *srcfit;  // points for source plane fitting
struct z_lim   zlim_s[NFMAX];
struct pot     lmin_s[NLMAX], lmax_s[NLMAX], prec_s[NLMAX];
struct sigposStr sigposAs_s, sigposAs;
struct matrix  amplifi_mat[NFMAX][NIMAX], amplifi_matinv[NFMAX][NIMAX];
struct pixlist plo[100000]; /*list of points that compose the polygon around the image to inverse*/

double z_dlsds;
double chi_im, chip, chix, chiy, chia, chil, chis, chi_vel,chi_mass;    //I added chi_vel and chi_mass TV
//double amplifi[NFMAX][NIMAX];
double **imo, **wo, **soo, **ero;
double ***cubeo,***ercubeo;
double **fluxmap;
double drima;   /* distance ratio between the image to inverse and the main lens*/
double **sm_p;      /* optimization spline map : initial save */

int    nshmap;
int    optim_z;
int    block_s[NLMAX][NPAMAX];
int  **imuo;
/*number of points of the polygon around the image to inverse*/

double distmin[NFMAX]; /* List of minimal euclidean distances between
                          the images of a same familly i*/

int    nwarn;   // warning emitted during image plane chi2

double *np_b0;  // non parametric accelerator
int    nplo;    // number of points for criticinv.c

/****************************************************************/
/*      nom:        global              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void  o_global()
{

    /*************  declaration de common et locale ****************/

    extern struct g_mode    M;
    extern struct g_grille  G;
    extern struct g_source  S;
    extern struct g_image   I;
    extern struct pot       lens[];


    double evidence;    // in the bayesian optimisation
    int  constraints, parameters; // for the computation of dof


    /**************************************************************/

    NPRINTF(stderr, "------- GLOBAL - Lens Optimisation - JPK, March 1994 ------\n");

    /* preparation de l'optimisation */
    readConstraints();

    /*-------------------------------------------------------*/
    /* optimisation */
    NPRINTF(stderr, "OPTIMISATION");
    if ( G.nmsgrid != G.nlens )
        NPRINTF(stderr, "(NON PARAM)");

    NPRINTF(stderr, " :\n");

    // set the evidence and compute the number of dof.
    evidence = 0;
    if ( M.inverse >= 3 ) NPRINTF( stderr, "Bayesys Rate : %lf \n", M.rate);

    constraints = getNConstraints();
    parameters = getNParameters();
    NPRINTF(stderr, "Number of contraints : %d\n", constraints);
    NPRINTF(stderr, "Number of free parameters : %d\n", parameters);


//    if ( parameters > constraints )
//    {
//        NPRINTF( stderr, "WARN : Too many free parameters.\n");
//      exit(-1);
//      return;
//    }

    // Start the optimisation engines
    if ( M.inverse == 3 || M.inverse == 4 )
        evidence = o_run_bayes();
    else if ( M.inverse == 5 )
     {
//        evidence = o_run_nest();
     }
    else
    {
        printf("ERROR : inverse modes 1 & 2 are not supported anymore\n");
        exit(-1);

//        /* sauvegarde des parametres d'optimisations */
//
//        /* is the potentiel a spline? */
//        if (lens[0].type != 10)
//        {
//            for (i = 0; i < G.no_lens; i++)
//            {
//                lmin_s[i] = lmin[i];
//                lmax_s[i] = lmax[i];
//                prec_s[i] = prec[i];
//                for (j = 0; j < NPAMAX; j++)
//                    block_s[i][j] = block[i][j];
//            }
//        }
//        else
//        {
//            sm_p = (double **) alloc_square_double (G.nx, G.ny);
//            sm_p = map_p;
//        }
//
//        for (i = 0; i < I.nzlim; i++)
//            zlim_s[i] = zlim[i];
//
//        sigposAs_s = sigposAs;
//
//        /* selection/identification des parametres pour la carte de chi2 */
//        // if any block[ils][ipx] is < 0, set the ip.map global variable
//        o_set_map();
//        if ( ip.map > 0 ) NPRINTF( stderr, "map dimensions : %d \n", ip.map);
//
//        if (M.inverse == 1) /* direct fitting */
//        {
//            if (mc.optMC == 1)
//                o_run_mc();
//            else if (mc.optMC == 2)
//                o_run_mc0();
//            else if (zlim[0].opt == 2)
//                o_runz_mc();
//            else if (zlim[0].opt == 3)
//                o_runz_mc0();
//            else if (ip.map == 0)
//                o_run();
//            else if (ip.map == 1)
//                o_run1();
//            else if (ip.map == 2)
//                o_run2();
//        }
//        else if (M.inverse == 2)  /* potfile fitting */
//        {
//            if ((P.ircut > 1) && (P.isigma > 1))
//                o_runpot2();
//            else if (P.ircut > 1)
//                o_runpot1(1);
//            else if (P.isigma > 1)
//                o_runpot1(2);
//            else
//                o_run();
//        }
    }

    /*-------------------------------------------------------
     * Reset the lens parameters and print the optimisation results */

    /*Fill the chires.dat file*/
    o_chires("chires.dat");

    /*Reset the lens parameters*/
    set_res_par();

    /*Print the results of optimisation*/
    o_print_res(o_chi(), evidence);

//  NPRINTF(stderr,"chi: %.1lf (%.1lf,%.1lf)\n",chi0,chip,chil);

    /*-------------------------------------------------------*/
    /*  resultats  */

    if (lens[0].type != 10)
    {
        int    i;
        for (i = 0; i < G.no_lens; i++)
        {
            if (lens[i].type > 1)
                lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
            else
                lens[i].sigma = sqrt(lens[i].b0 / 4. / pia_c2);

            lens[i].dlsds = dratio(lens[i].z, S.zs);
        }
    }

    if (I.stat != 0)
    {
        extern struct galaxie   arclet[];
        extern long int    narclet;
        o_stat(narclet, arclet);
    }

    o_global_free();
    NPRINTF(stderr, "------- Lens Optimisation Done -----------\n");
}

/*
 * Free all memory in o_global.readConstraints()
 */
void o_global_free()
{
    extern struct g_image   I;
    extern struct g_grille  G;
    extern struct g_mode    M;
    
    if ( I.srcfit != 0 )
        free(srcfit);

    if ( G.nlens != G.nmsgrid )
    {
        extern struct galaxie   arclet[];

        free(np_b0);
        if( I.n_mult > 0 )
        {
            free(multi[0][0].np_grad);  // root of the malloc
            free(multi[0][0].np_grad2a);
            free(multi[0][0].np_grad2b);
            free(multi[0][0].np_grad2c);
        }
        else
        {
            free(arclet[0].np_grad);
            free(arclet[0].np_grad2a);
            free(arclet[0].np_grad2b);
            free(arclet[0].np_grad2c);
        }
    }

    if (M.iclean != 0)
    {
        extern struct g_source  S;
        const extern struct g_pixel   ps;
        extern struct g_pixel   imFrame;
        extern struct g_pixel   wFrame;
        extern struct g_cube    cubeFrame;
        extern double **map_axx, **map_ayy;
        int i;

        //  free the square maps
        if (M.iclean == 2 )
        {
            free_square_double(ero, imFrame.ny);
        }
        else if (M.iclean == 3 )
        {
            free_cubic_double(cubeo, cubeFrame.ny);
            free_cubic_double(ercubeo, cubeFrame.ny, cubeFrame.nx);
        }
        else
        {
            free_square_double(soo, ps.ny);
            free_square_int(imuo, ps.ny);
        }

        free_square_double(imo, imFrame.ny);
        free_square_double(wo, wFrame.ny);
        free_square_double(fluxmap, imFrame.ny*imFrame.nx);

        int pot_nopt = 0;
        extern struct g_pot P[NPOTFILE];
        for (i = 0; i < G.npot; i++ )
        {
            pot_nopt += P[i].ircut;
            pot_nopt += P[i].isigma;
            pot_nopt += P[i].islope;
            pot_nopt += P[i].ivdslope;
            pot_nopt += P[i].ivdscat;
            pot_nopt += P[i].ircutscat;
            pot_nopt += P[i].ia ;
            pot_nopt += P[i].ib;
        }

        if (G.no_lens == 0 && pot_nopt == 0)
        {
            free_square_double(map_axx, imFrame.ny);
            free_square_double(map_ayy, imFrame.ny);
        } 

        // Reset source positions to absolute positions
        if( I.n_mult > 0 )
        {
            extern struct galaxie source[NFMAX];
            struct point Bs;  // list of source barycenters  
            struct point Ps[NIMAX];
            char str[IDSIZE], *pch;  // get image identifier
            int    j;

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
                    Bs = bcentlist(Ps, I.mult[j]);
                    source[i].C.x += Bs.x;
                    source[i].C.y += Bs.y;
                }
            }
        }

        // Reset source counter to zero
        S.ns = 0;
    }
}

/* Read the different constraints to prepare the optimisation
 * according to the global variables set in the .par file.
 */
void readConstraints()
{
    extern struct g_mode    M;
    extern struct g_image   I;
    extern struct g_source  S;
    extern struct g_grille  G;
    extern struct pot       lens[];
    extern struct cline     cl[];
    extern struct sigposStr sigposAs;
    extern struct galaxie   multi[NFMAX][NIMAX];

    extern struct galaxie   arclet[];
    extern long int    narclet;

    int  np;
    long int ntmult;
    int    i, j;
    FILE    *IN;

    if (M.iclean != 0)
    {
//        imo = (double **) o_invim(M.zclean, &drima, plo, &nplo);
        extern struct g_pixel imFrame;
        extern struct g_pixel   wFrame;
        extern struct g_cube    cubeFrame;
        extern struct galaxie source[NFMAX];
        const extern struct g_pixel   ps;
        extern double **map_axx, **map_ayy;

        // Read input image from cleanlens section
        if(M.iclean<3)
        {
		    imo = (double **) readimage(&imFrame);
        }
        if(M.iclean==3)
        {
            cubeo = (double ***) readcube(&cubeFrame);
        }
//        for( i = 0; i < imFrame.nx; i++ )
//            for( j = 0; j < imFrame.ny; j++ )
//                if( imo[j][i] == 0. )
//                {
//                    fprintf(stderr, "ERROR: O value found in %s. Cannot compute chi2.\n", imFrame.pixfile);
//                    exit(-1);
//                }


        if(wFrame.format != 0)
        {
            wo = (double **) readimage(&wFrame);
            
            
            // Test of positivity of wo 
            for( i = 0; i < wFrame.nx; i++ )
                for( j = 0; j < wFrame.ny; j++ )
                    if( wo[j][i] < 0. )
                    {
                        fprintf(stderr, "ERROR: negative value found in %s.\n", wFrame.pixfile);
                        exit(-1);
                    }

            // Test of dimension matching
            if ( wFrame.nx != imFrame.nx || wFrame.ny != imFrame.ny )
                {
                    fprintf(stderr, "ERROR:  dimensions mismatch in %s and %s", wFrame.pixfile, imFrame.pixfile);        
                    exit(-1);
                }
        }

        if ( M.iclean == 2 || M.iclean ==3 )
        {

            if(M.iclean==2)
            {
              // allocate temporary array to compute chi2 in o_chi()
              ero = (double **) alloc_square_double(imFrame.ny, imFrame.nx);
            }
            if(M.iclean==3)
            {
              // allocate temporary array to compute chi2 in o_chi()
              ercubeo = (double ***) alloc_cubic_double(cubeFrame.ny, cubeFrame.nx, cubeFrame.nz);
            }

            int pot_nopt = 0;
            extern struct g_pot P[NPOTFILE];
            for (i = 0; i < G.npot; i++ )
            {
                pot_nopt += P[i].ircut;
                pot_nopt += P[i].isigma;
                pot_nopt += P[i].islope;
                pot_nopt += P[i].ivdslope;
                pot_nopt += P[i].ivdscat;
                pot_nopt += P[i].ircutscat;
                pot_nopt += P[i].ia ;
                pot_nopt += P[i].ib;
            }

            if( G.no_lens == 0 && pot_nopt == 0 )
            {
                map_axx = (double **) alloc_square_double(imFrame.ny, imFrame.nx);
                map_ayy = (double **) alloc_square_double(imFrame.ny, imFrame.nx);
                struct point pi, ps;
                for (i = 0; i < imFrame.ny; i++)
                {
                    pi.y = imFrame.ymin + i * imFrame.pixelx;
                    for (j = 0; j < imFrame.nx; j++)
                    {
                        pi.x = imFrame.xmin + j * imFrame.pixelx;
                        e_dpl(&pi, 1., &ps);
                        map_axx[i][j] = pi.x - ps.x;
                        map_ayy[i][j] = pi.y - ps.y;
                    }
                }
            }

            // read and prepare sources for o_chi()
            if( M.source > 0 )
            {
                long int istart = S.ns;
                f_shape(&S.ns, source, M.sourfile, M.source);  // with var1 & var2 parameters
                if( M.source == 2 )
                    for( i = istart; i < S.ns; i++ )
                        source[i].type = source[i].var1;

                pro_arclet(S.ns, source);  // compute eps = (a-b)/(a+b)
            }

            // prepare sources to have relative positions wrt barycenter
//            if( I.n_mult > 0 )
//            {
//                struct point Ps[NIMAX];
//                char str[IDSIZE], *pch;  // get image identifier
//
//                for( i = 0; i < S.ns; i++ )
//                {
//                    // check if source is attached to an image
//                    strcpy(str, source[i].n);
//                    pch = strtok(str, "_");
//                    pch = strtok(NULL, "_");
//                    if( pch == NULL )    break;  // no attachment to an image
//
//                    // look for the corresponding image family
//                    j = 0; 
//                    while ( indexCmp( pch, multi[j][0].n ) && j < I.n_mult ) j++;
//
//                    // set source position into relative coordinates
//                    if( j < I.n_mult )
//                        source[i].C.x = source[i].C.y = 0.;
//                }
//            }

            // prepare source redshifts
            for ( i = 0; i < S.ns; i++ )
            {
                if( source[i].z == 0 )
                {
                    source[i].z = M.zclean;
                    NPRINTF(stderr, "WARN: source %s redshift set to cleanset z=%lf\n", source[i].n, M.zclean);
                }
                dratio_gal(&source[i], lens[0].z);
            }
        }
        else
        {
            //  allocate square map
            soo = (double **) alloc_square_double(ps.ny, ps.nx);
            ero = (double **) alloc_square_double(ps.ny, ps.nx);
            imuo = (int **) alloc_square_int(ps.ny, ps.nx);
        }

    }

    if (I.shmap != 0)
    {
        nshmap = 0;
        f_shmap(&nshmap, shm, I.shfile);
        if (nshmap < 1)
            I.shmap = 0;

        I.dl0ssh = distcosmo2( lens[0].z , I.zsh );
        I.dossh = distcosmo1( I.zsh );
        I.drsh = I.dl0ssh / I.dossh;
    }

    if (I.stat != 0)
    {
        narclet = 0;

        // Test if the number of arclets is lower than NAMAX
        IN = fopen(I.arclet, "r");
        if ( IN == NULL )
        {
            fprintf( stderr, "ERROR: File %s not found.\n", I.arclet );
            exit(-1);
        }

        if ( wc(IN) > NAMAX )
        {
            fprintf( stderr, "ERROR: Too many arclets in %s (%d). Max allowed %d\n",
                     I.arclet, (int)wc(IN), NAMAX);
            fclose(IN);
            exit(-1);
        }
        fclose(IN);

        f_shape(&narclet, arclet, I.arclet, I.statmode);

        if (M.seeing > 0.)
            cor_seeing(narclet, arclet, M.seeing);

        if (narclet < 1)
            I.stat = 0;
        else
        {
            sort(narclet, arclet, comparer_z);
            zalim.opt = 0;
            gsl_rng *seed = gsl_rng_alloc(gsl_rng_taus);
            gsl_rng_set(seed, S.rand);
            gsl_ran_discrete_t *gsmail;
            if( S.distz == 2 && S.par1 != 0 )
                gsmail = smailpreproc();

            while( arclet[zalim.opt].z == 0 && zalim.opt < narclet ) 
            {
                if( I.zarclet > 0 )
                    arclet[zalim.opt].z = I.zarclet;
                else if( zalim.min > 0 )
                    arclet[zalim.opt].z = zalim.min;
                else if( S.distz == 0 )
                    arclet[zalim.opt].z = S.ns;
                else if( S.distz == 1 )
                    arclet[zalim.opt].z = S.zsmin;
                else if( S.distz == 2 && S.par1 != 0 )
                    arclet[zalim.opt].z = d_rndzsmail(seed, gsmail);

                if( arclet[zalim.opt].z > 0 )
                    zalim.opt++;
                else
                {
                    fprintf(stderr, "ERROR: Arclets with unknown redshift, and no redshift set in image section\n");
                    exit(E_ZALIM_MISSING);
                }
            }

            if( zalim.opt > 0 )
                FPRINTF(stdout, "INFO: Set %d arclets with unknown redshift\n", zalim.opt);
  
            for( i = 0; i < narclet; i++ )
                dratio_gal(&arclet[i], lens[0].z);

            pro_arclet(narclet, arclet);

            gsl_rng_free(seed);
            if( S.distz == 2 && S.par1 != 0 )
                gsl_ran_discrete_free(gsmail);
        }

        // Intrinsic ellipticity distribution width:
        I.sig2ell = I.sigell * I.sigell;

    }

    if (I.n_mult != 0)
    {
        ntmult = 0;
        struct galaxie *mult = (struct galaxie *) malloc((unsigned long int) NFMAX * NIMAX * sizeof(struct galaxie));
        if ( mult == NULL )
        {
            fprintf( stderr, "ERROR: mult[NFMAX*NIMAX] memory allocation failed.\n");
            exit(-1);
        }

        if (I.n_mult == 2 )
            f_shape(&ntmult, mult, I.multfile, 2);
        else 
        {
            f_shape(&ntmult, mult, I.multfile, 1);
            for ( i = 0; i < ntmult; i++ )
                mult[i].var1 = mult[i].var2 = 1.;  //weight = 1
        }

        if (M.seeing > 0.)
            cor_seeing(ntmult, mult, M.seeing);

        if (ntmult < 1)
            I.n_mult = 0;
        else
        {
            pro_arclet(ntmult, mult);
            o_prep_mult(ntmult, mult);
            dist_min();

            // Set the image position error in image plane
            for ( i = 0; i < I.n_mult; i++ )
                for ( j = 0; j < I.mult[i]; j++ )
                    I.sig2pos[i][j] = sigposAs.min * sigposAs.min;

            if( I.forme == 12)
                readWeightMat(&I.weight, &I.detCov);
        }
        free(mult);
    }

    // Prepare Accelerated potentials
    if ( G.nmsgrid != G.nlens )
    {
        prep_non_param();
 
        // Error message to prevent deletion of imported data
        if( M.source != 0 && I.n_mult != 0 )
        {
            fprintf(stderr, "ERROR: You cannot import a catalog of sources while optimizing with the grid.\n");
            exit(1);
        }
    }


    // Source plane fit of a brightness profile
    if (I.srcfit != 0 )
    {
        srcfit = (struct galaxie *)malloc((unsigned int)NSRCFIT * sizeof(struct galaxie));
        if ( srcfit == NULL )
        {
            fprintf(stderr, "ERROR: srcfit memory allocation failed.\n");
            exit(1);
        }

        f_shape(&I.nsrcfit, srcfit, I.srcfitFile, 0);
        for (i = 0; i < I.nsrcfit; i++)
            dratio_gal(&srcfit[i], lens[0].z);
    }

    if (I.npcl != 0)
    {
        /* preparation de l'optimisation dans le cas de contraintes de ligne */
        /* critiques (cassures d'arcs) */
        np = 0;
        for (i = 0; i < I.npcl; i++)
        {
            if (cl[i].n != 0)
            {
                if (lens[0].z < cl[i].z)
                    cl[i].dl0s = distcosmo2(lens[0].z, cl[i].z);
                else
                    cl[i].dl0s = 0.;

                cl[i].dos = distcosmo1(cl[i].z);
                cl[i].dlsds = cl[i].dl0s / cl[i].dos;

                cl[np] = cl[i]; // stack all the CL contraints at the beginning of <cl> list
                np++;
            };
        }
        I.npcl = np;  // new number of CL constraints
        NPRINTF(stderr, "INF: critical points: %d\n", I.npcl);
    }
}


// Compute the number of constraints with multiple images
int getNConstraints()
{
    extern struct g_mode  M;
    extern struct g_image I;
    extern struct g_dyn   Dy;   //TV Oct 2011
    extern long int narclet;
    int constraints, i;
    int opi;        // number of observable per images

    opi = 2;     // default only consider x,y
    if ( I.forme == 1 ) opi = 4;  // +shape
    if ( I.forme == 2 ) opi = 3;  // +flux only
    if ( I.forme == 3 ) opi = 5;  // +flux + shape

    constraints = 0;
    // Multiple images constraints
    for ( i = 0; i < I.n_mult; i ++)
        if ( I.mult[i] > 1 )
            constraints += opi * (I.mult[i] - 1);
        else  // for singly imaged galaxy
            constraints += opi; 

    // Sums up the critical lines
    constraints += 2 * I.npcl;

    if ( I.srcfit != 0 )
        constraints += 2 * I.nsrcfit;

    // arclets
    if ( I.stat > 0 )
        constraints += 2 * narclet;

    //dynamics                     // TV
    if(Dy.dynnumber == 1)
        constraints = constraints+1;
    
    if(Dy.dynnumber == 2)
        constraints = constraints+2;
    
    if(Dy.dynnumber == 3)
        constraints = constraints+1;
    
    if(Dy.dynnumber == 4)
        constraints = constraints+1;
	

    // FITS image reconstruction
    if ( M.iclean == 2 )
    {
        extern struct g_pixel imFrame;
        constraints += imFrame.nx * imFrame.ny;
    }

    return(constraints);
}


// Compute the number of free parameters
int getNParameters()
{
    extern  struct  g_grille    G;
    extern  struct  g_pot       P[NPOTFILE];
    extern  struct  g_image     I;
    extern  struct  g_source    S;
    extern  int     block[][NPAMAX];
    extern  int     cblock[], sblock[NFMAX][NPAMAX];
    extern  struct  sigposStr sigposAs;

    long int parameters, i, j;
    int ipx;

    parameters = 0;
    for ( i = 0; i < G.no_lens; i ++ )
        for ( ipx = CX; ipx <= PMASS; ipx ++ )
            if ( block[i][ipx] != 0 )
                parameters++;

    // multiscale grid clumps
    parameters += G.nlens - G.nmsgrid;

    // source parameters
    for ( i = 0; i < S.ns; i++ )
        for (ipx = SCX; ipx <= SFLUX; ipx++ )
            if ( sblock[i][ipx] != 0 )
                parameters++;

    // cosmological parameters
    for ( ipx = OMEGAM; ipx <= WA; ipx++ )
        if ( cblock[ipx] ) parameters++;

    // redshift optimization
    for ( i = 0; i < I.nzlim; i++ )
        if( zlim[i].bk > 0 )
            parameters ++;

    // source fitting
    if ( I.srcfit != 0 )
    {
        if ( !strcmp(I.srcfitMethod, "LINFIT") )
            parameters += 2;
    }

    for ( i = 0; i < G.npot; i++ )
    {
        struct g_pot *pot = &P[i];

        if ( pot->ftype && pot->ircut != 0 ) parameters++;
        if ( pot->ftype && pot->isigma != 0 ) parameters++;
        if ( pot->ftype && pot->islope != 0 ) parameters++;
        if ( pot->ftype && pot->ivdslope != 0 ) parameters++;
        if ( pot->ftype && pot->ivdscat != 0 ) parameters++;
        if ( pot->ftype && pot->ircutscat != 0 ) parameters++;
        if ( pot->ftype && pot->ia != 0 ) parameters++;
        if ( pot->ftype && pot->ib != 0 ) parameters++;
    }

    if ( sigposAs.bk != 0 )
    {
        for ( i = 0; i < I.n_mult; i++ )
            for ( j = 0; j < I.mult[i]; j++ )
                parameters++;
    }
    if ( I.dsigell != -1. ) parameters++;

    return(parameters);
}

/* Initialise the np_grad and np_grad2 global variables
 * used after in e_grad() and e_grad2().
 *
 * These variables are used to speed up the computation
 * of the gradients for the many identical potentials
 * involved in the non parametric model.
 *
 */
void prep_non_param()
{
    extern struct g_image  I;
    extern struct g_grille G;
    extern struct pot      lens[];
    extern long int        narclet;
    extern struct galaxie  arclet[];
    extern struct pot      lmax[];

    struct galaxie *image;
    long int i, j, k, l, n, nimages;
    struct point *pGrad;
    double  *pGrad2a, *pGrad2b, *pGrad2c, *tmp;
    double dls, oldz;
    struct matrix grad2;
    long int gcount, lcount;  // global and local counters
    struct ellipse ampli;

    // Number of clumps
    n = G.nlens - G.nmsgrid;

    // Create a global array for the lens.b0
    np_b0 = (double *) calloc((unsigned) n, sizeof(double));

    // set all lens[*].b0 so that all lens[*].pmass=1
    for ( i = 0; i < G.nlens - G.nmsgrid; i++ )
    {
        np_b0[i] = lens[i + G.nmsgrid].b0;
        lens[i + G.nmsgrid].b0 = 1.;
//        lens[i + G.nmsgrid].b0 = 0.;
//        for( j = 0; j < G.nlens - G.nmsgrid; j++ )
//            lens[i + G.nmsgrid].b0 += G.invmat[i][j];
    }

    // Initialise np_grad and np_grad2 arrays
    nimages = narclet; 
    for ( i = 0; i < I.n_mult; i++ )
        nimages += I.mult[i];

    pGrad = (struct point *) malloc((unsigned int) n * nimages * sizeof(struct point));
    pGrad2a = (double *) malloc((unsigned int) n * nimages * sizeof(double));
    pGrad2b = (double *) malloc((unsigned int) n * nimages * sizeof(double));
    pGrad2c = (double *) malloc((unsigned int) n * nimages * sizeof(double));
    tmp = (double *) malloc((unsigned int) n * nimages * sizeof(double));

    // Check memory allocation
    if( pGrad == NULL || pGrad2a == NULL || pGrad2b == NULL || pGrad2c == NULL || tmp == NULL )
    {
        fprintf(stderr, "ERROR in prep_non_param() during memory allocation\n");
        exit(1);
    } 

    nimages = 0;
    for ( i = 0; i < I.n_mult; i++ )
        for ( j = 0; j < I.mult[i]; j++ )
        {
            image = &multi[i][j];

            image->np_grad = pGrad + nimages * n;
            image->np_grad2a = pGrad2a + nimages * n;
            image->np_grad2b = pGrad2b + nimages * n;
            image->np_grad2c = pGrad2c + nimages * n;
            nimages++;

            oldz = lens[0].z; dls = image->dl0s;
            for ( k = G.nmsgrid; k < G.nlens; k++ )
            {
                if( lens[k].z != oldz )
                {
                    dls  = distcosmo2(lens[k].z, image->z);
                    oldz = lens[k].z;
                }

                // dls/ds multiplication is done in o_dpl.c for e_grad_pot()
                image->np_grad[k-G.nmsgrid] = e_grad_pot(&image->C, k);

                grad2 = e_grad2_pot(&image->C, k); 
                image->np_grad2a[k - G.nmsgrid] = grad2.a * dls;
                image->np_grad2b[k - G.nmsgrid] = grad2.b * dls;
                image->np_grad2c[k - G.nmsgrid] = grad2.c * dls;
            }
        }

    // For SL: Compute amplification due to the fixed potentials
    if( I.n_mult > 0 )
    {
        extern struct g_pot P[NPOTFILE];
        struct matrix *amatinv; // correction factor for the error, which is not exactly in the image plane if there are fixed potentials
        int nimage = 0;
        long int ilens;
        for ( i = 0; i < I.n_mult; i++ )
            nimage += I.mult[i];
    
        amatinv = (struct matrix *)calloc((unsigned int)nimage, sizeof(struct matrix));
    
        // Contribution from the fixed potentials
        for( ilens = G.no_lens; ilens < G.nplens[0]; ilens++ )
        {
            nimage = 0;
            for ( i = 0; i < I.n_mult; i++ )
                for ( j = 0; j < I.mult[i]; j++ )
                {
                    grad2 = e_grad2_pot(&multi[i][j].C, ilens);
                    amatinv[nimage].a += grad2.a * multi[i][j].dr;
                    amatinv[nimage].b += grad2.b * multi[i][j].dr;
                    amatinv[nimage].c += grad2.c * multi[i][j].dr;
                    nimage++;
                }
        }
        
        // Contribution from the potfile, if it is fixed 
        for( k = 0; k < G.npot; k++ )
            if( P[k].isigma == 0 && P[k].ircut == 0 )
                for( ilens = G.nplens[k]; ilens < G.nplens[k+1]; ilens++ )
                {
                    nimage = 0;
                    for ( i = 0; i < I.n_mult; i++ )
                        for ( j = 0; j < I.mult[i]; j++ )
                        {
                            grad2 = e_grad2_pot(&multi[i][j].C, ilens);
                            amatinv[nimage].a += grad2.a * multi[i][j].dr;
                            amatinv[nimage].b += grad2.b * multi[i][j].dr;
                            amatinv[nimage].c += grad2.c * multi[i][j].dr;
                            nimage++;
                        }
                }
        
        
        // Compute resulting amplification
        nimage = 0;
        for ( i = 0; i < I.n_mult; i++ )
            for ( j = 0; j < I.mult[i]; j++ )
            {
                multi[i][j].mu = sqrt(fabs((1. - amatinv[nimage].a) * (1. - amatinv[nimage].c) - amatinv[nimage].b * amatinv[nimage].b));
                nimage++;
            }
        
        free(amatinv);

    }

    //
    // Process arclets with OPENMP
    //
    gcount = lcount = 0; // count number of processed arclets (big/small loop)

//#pragma omp parallel default(shared) firstprivate(pGrad,pGrad2a,pGrad2b,pGrad2c) private(i,image,k,lcount,grad2)
    {
//#pragma omp for nowait
        for ( i = 0; i < narclet; i++ )
        {
            lcount++;
            if ( lcount == 200 )
            {
//#pragma omp critical
                {
                    gcount = gcount + lcount;
                    printf("INFO: prepare lens profiles for arclet %ld/%ld\r", gcount, narclet);
                }

                lcount = 0;
            }

            image = &arclet[i];

            image->np_grad = pGrad + nimages * n;
            image->np_grad2a = pGrad2a + nimages * n;
            image->np_grad2b = pGrad2b + nimages * n;
            image->np_grad2c = pGrad2c + nimages * n;
            nimages++;

            oldz = lens[0].z; dls = image->dl0s;
            for ( k = G.nmsgrid; k < G.nlens; k++ )
            {
                if( lens[k].z != oldz )
                {
                    dls  = distcosmo2(lens[k].z, image->z);
                    oldz = lens[k].z;
                }

                image->np_grad[k-G.nmsgrid] = e_grad_pot(&image->C, k);
                image->np_grad[k-G.nmsgrid].x *= image->dr;
                image->np_grad[k-G.nmsgrid].y *= image->dr;
                grad2 = e_grad2_pot(&image->C, k);
                image->np_grad2a[k-G.nmsgrid] = grad2.a * dls;
                image->np_grad2b[k-G.nmsgrid] = grad2.b * dls;
                image->np_grad2c[k-G.nmsgrid] = grad2.c * dls;
            }
        }
    }
    if( narclet > 0 )
        printf("INFO: prepare lens profiles for arclet %ld/%ld\n", narclet, narclet);

    // restore lens[i].b0 values
    for ( i = G.nmsgrid; i < G.nlens; i++ )
        lens[i].b0 = np_b0[i-G.nmsgrid];

    // clean allocated data not used afterwards
    free(tmp);
}

// Read a FITS file containing the weight matrix
// Return the determinant of the Covariance matrix, for Likelihood normalization
static void readWeightMat(double *** weight_ext, double *detCov)
{
    int i, j, nx, ny, nullity;
    char *header;
    double *p, result, **tmp, **weight;

    weight = rdf_fits("weight.fits", &nx, &ny, &header);
    p = (double *) calloc(nx, sizeof(double));
    tmp = (double **) malloc(nx * ny * sizeof(double));
    memcpy( tmp, weight, nx * ny * sizeof(double));

    cholesky(tmp, nx, p);

    result = 1.;
    for( i = 0; i < nx; i++ )
        result *= p[i] * p[i];

    *detCov = 1. / result; 
    *weight_ext = weight;

    free(tmp);
    free(p);
    free(header);
}
