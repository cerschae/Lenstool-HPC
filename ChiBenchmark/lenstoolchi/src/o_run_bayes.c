#include<stdio.h>
#include<signal.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "bayesys3.h"
#include "userstr.h"
#include "lt.h"

/****************************************************************/
/*      nom:        o_run_bayes         */
/*      auteur:     Eric Jullo          */
/*      date:       01/06               */
/*      place:      ESO, Chile          */
/****************************************************************
 * Optimise the potential parameters using the bayesys3 algorithm.
 *
 */
typedef void (*sighandler_t)(int);
static void signalReset(int);
int optInterrupt;   // Global variable read in bayesapp.c/UserMonitor()
double  *tmpCube;      // temporary Cube to be used in UserBuild()
double  **zbitstmp1;
int    *ibitstmp;

// Mutex functions implemented in bayesapp.c
void init_mutex();
void destroy_mutex();

double  o_run_bayes()
{
    extern struct g_mode    M;
    extern struct g_image   I;
    extern struct g_grille  G;
// Addressing
    CommonStr     Common[1];
    ObjectStr     *Objects;
    UserCommonStr UserCommon[1];
    int           nDim;
    long int      i, j;
    FILE *burnin;


    // Get the number of free parameters, and write header of bayes.dat
    nDim = bayesHeader();

#ifdef DEBUG
    burnin = fopen( "burnin.dbg.dat", "w" );
#else
    burnin = fopen( "burnin.dat", "w" );
#endif
    fclose(burnin);
 

    //Common Bayesys parameters
    Common->Ndata = 0;
    Common->ENSEMBLE    = 10;         // # objects in ensemble
    if ( M.minchi0 > 10 )
        Common->ENSEMBLE    = (int)ceil(M.minchi0);
    Common->Method      = -1;         // Algorithm method
    Common->Rate        = M.rate;     // Speed of calculation (dimensionless)
    Common->Iseed       = -4321;       // Random seed, -ve is time seed
#if defined(DEBUG) || defined(__CONSTANT_SEED__)
    Common->Iseed     = 4321;     // Reproductible samples
#endif
    Objects = (ObjectStr *)malloc((unsigned int) Common->ENSEMBLE * sizeof(ObjectStr));


    // Initialize Bayesys for MassInf. WL only with a Grid
    if( G.nmsgrid != G.nlens )
    {
        if( I.stat == 8 || I.n_mult != 0 )
        {
            extern long int narclet;
            extern struct g_pot P[NPOTFILE];
            extern struct galaxie   arclet[];
            extern struct galaxie  multi[NFMAX][NIMAX];
            extern struct g_frame F;
            double gam1, gam2;
            struct matrix grad2;
            struct point grad;
            long int ilens;
            double e;
            int k;
    
            int Ndata = narclet;
            for ( i = 0; i < I.n_mult; i++ )
                Ndata += I.mult[i];
            int nmult = Ndata - narclet;

            Ndata = 2 * Ndata;
            double *Data = (double *)malloc((unsigned int)(Ndata * sizeof(double)));
            double *Acc = (double *)malloc((unsigned int)(Ndata * sizeof(double)));
            double *zbitstmp0 = (double *)calloc(Ndata, sizeof(double));
            zbitstmp1 = alloc_square_double(G.nlens - G.nmsgrid + 2 * I.n_mult, Ndata);
            ibitstmp = (int *)calloc(Ndata, sizeof(int));

            // Contribution from the grid potentials
            for( ilens = 0; ilens < G.nlens - G.nmsgrid; ilens++ )
            {
                int nimage = 0;
                for ( i = 0; i < I.n_mult; i++ )
                    for ( j = 0; j < I.mult[i]; j++ )
                    {
                        zbitstmp1[ilens][nimage] = multi[i][j].np_grad[ilens].x * multi[i][j].dr;
                        zbitstmp1[ilens][nmult + nimage] = multi[i][j].np_grad[ilens].y * multi[i][j].dr;
                        nimage++;
                    }

                for( i = 0; i < narclet; i++ )
                {
                    gam1 = 0.5 * (arclet[i].np_grad2a[ilens] - arclet[i].np_grad2c[ilens]) / arclet[i].dos; 
                    gam2 = arclet[i].np_grad2b[ilens] / arclet[i].dos; 
                    zbitstmp1[ilens][2 * nmult + i] = 2. * gam1;   // e1 defined as (a^2 - b^2) / (a^2 + b^2) and kap <<1
                    zbitstmp1[ilens][2 * nmult + narclet + i] = 2. * gam2; // g2 (see B&S01 Eq 4.6)
                }
            }
     
            // Contribution from the fixed potentials
            int nimage = 0;
            for( ilens = G.no_lens; ilens < G.nplens[0]; ilens++ )
            {
                nimage = 0;
                for ( i = 0; i < I.n_mult; i++ )
                    for ( j = 0; j < I.mult[i]; j++ )
                    {
                        grad = e_grad_pot(&multi[i][j].C, ilens);
                        zbitstmp0[nimage] += grad.x * multi[i][j].dr;
                        zbitstmp0[nmult + nimage] += grad.y * multi[i][j].dr;
                        nimage++;
                    }

                for( i = 0; i < narclet; i++ )
                {
                    grad2 = e_grad2_pot(&arclet[i].C, ilens);
                    zbitstmp0[2 * nmult + i] += (grad2.a - grad2.c) * arclet[i].dr;
                    zbitstmp0[2 * nmult + narclet + i] += 2. * grad2.b * arclet[i].dr;
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
                                grad = e_grad_pot(&multi[i][j].C, ilens);
                                zbitstmp0[nimage] += grad.x * multi[i][j].dr;
                                zbitstmp0[nmult + nimage] += grad.y * multi[i][j].dr;
                                nimage++;
                            }
    
                        for( i = 0; i < narclet; i++ )
                        {
                            grad2 = e_grad2_pot(&arclet[i].C, ilens);
                            zbitstmp0[2 * nmult + i] += (grad2.a - grad2.c) * arclet[i].dr;
                            zbitstmp0[2 * nmult + narclet + i] += 2. * grad2.b * arclet[i].dr;
                        }
                    }

            // Shift source positions to (F.xmin, F.ymin)
//            nimage = 0;
//            double xmin = 1000.;
//            double ymin = 1000.;
//            for( i = 0; i < I.n_mult; i++ )
//                for( j = 0; j < I.mult[i]; j++ )
//                {
//                    if( multi[i][j].C.x < xmin )     xmin = multi[i][j].C.x; 
//                    if( multi[i][j].C.y < ymin )     ymin = multi[i][j].C.y;
//                }
//
//            for( i = 0; i < I.n_mult; i++ )
//                for( j = 0; j < I.mult[i]; j++ )
//                {
//                    zbitstmp0[nimage] += xmin;
//                    zbitstmp0[nmult + nimage] += ymin;
//                    nimage++;
//                }
 
            // Indexes for the source positions
            nimage = 0;
            for ( i = 0; i < I.n_mult; i++ )
                for ( j = 0; j < I.mult[i]; j++ )
                {
                    zbitstmp1[G.nlens - G.nmsgrid + i][nimage] = 1.;
                    zbitstmp1[G.nlens - G.nmsgrid + I.n_mult + i][nmult + nimage] = 1.;
                    nimage++;
                }
    
            // Data
            nimage = 0;
            for ( i = 0; i < I.n_mult; i++ )
                for ( j = 0; j < I.mult[i]; j++ )
                {
                    Data[nimage] = multi[i][j].C.x - zbitstmp0[nimage];
                    Data[nmult + nimage] = multi[i][j].C.y - zbitstmp0[nmult + nimage];
                    Acc[nimage] = 1. / multi[i][j].E.a / multi[i][j].mu;  // cos(multi[i][j].E.theta);
                    Acc[nmult + nimage] = 1. / multi[i][j].E.b / multi[i][j].mu; // cos(multi[i][j].E.theta);;

                    ibitstmp[nimage] = nimage;
                    ibitstmp[nmult + nimage] = nmult + nimage;
                    nimage++;
                }

            for( i = 0; i < narclet; i++ )
            {
                e = (arclet[i].E.a * arclet[i].E.a - arclet[i].E.b * arclet[i].E.b) / (arclet[i].E.a * arclet[i].E.a + arclet[i].E.b * arclet[i].E.b);
                Data[2 * nmult + i] = e * cos(2. * arclet[i].E.theta) - zbitstmp0[2 * nmult + i];  // gamma1
                Data[2 * nmult + narclet + i] = e * sin(2. * arclet[i].E.theta) - zbitstmp0[2 * nmult + narclet + i]; // gamma2
                Acc[2 * nmult + i] = 1. / sqrt(arclet[i].var1 + I.sig2ell);
                Acc[2 * nmult + narclet + i] = 1. / sqrt(arclet[i].var2 + I.sig2ell);
    
                ibitstmp[2 * nmult + i] = 2 * nmult + i;
                ibitstmp[2 * nmult + narclet + i] = 2 * nmult + narclet + i;
            }
    
            free(zbitstmp0);

            Common->Ndata       = Ndata;
            Common->Data        = Data;
            Common->Acc         = Acc;
           
            // UserBuild uses mock data
            for( i = 0; i < Common->ENSEMBLE; i++ ) 
                Objects[i].Mock = (double*)malloc(Ndata * sizeof(double));
    
        }
        
        // Number of Atoms free
        Common->MinAtoms    = 1;          // >= 1
        Common->MaxAtoms    = 0;          // >= MinAtoms, or 0 = infinity(not implemented)
        Common->Alpha       = -0.02 * (G.nlens - G.nmsgrid) - 2 * I.n_mult - 1;       // +ve for Poisson, -ve for geometric

        if( nDim < G.nlens )
        {
            // MassInf parameters
            Common->Ndim      = 1;     // Coordinate in the grid
            Common->Valency   = 1;        // On/off switch for MassInf/Grid optim (Default: off)
            Common->MassInf     = 1;    // Positive prior on the flux
            if( I.n_mult > 0 && G.no_lens > 0 )
            {
                Common->Valency = 3;    // Add 2 valencies for the SL source positions X and Y 
                Common->Method      = -3;         // Remove lifeStory2 algorithm for conflict of atom valencies
                Common->MassInf     = 2;    // Positive/Negative prior on the flux
            }
            Common->ProbON      = 0.7;
            Common->FluxUnit0   = -10.;  // 1 for pmass, 10 for b0
            nDim = G.nlens - G.nmsgrid;  // To initialize statistics arrays err, tmpCube and avg
        }
        else
        {
            // Normal optimization but with Natoms (see UserBuild)
            Common->Ndim      = nDim;     // Number of parameters of the optimized clumps
            Common->Valency   = 0;        // On/off switch for MassInf/Grid optim (Default: off)
        }
    }
    else
    { 
        // Number of Atoms = 1
        Common->MinAtoms    = 1;          // >= 1
        Common->MaxAtoms    = 1;          // >= MinAtoms, or 0 = infinity(not implemented)
        Common->Alpha       = -1.;       // +ve for Poisson, -ve for geometric

        Common->Ndim      = nDim;     // Number of parameters of the optimized clumps
        Common->Valency   = 0;        // On/off switch for MassInf/Grid optim (Default: off)
    }

    
    // Initialize some statistics
    Common->UserCommon = (void*)UserCommon;
    UserCommon->Nsample = 0;
    UserCommon->atoms   = 0.0;
    UserCommon->Nchi2 = 0;
    UserCommon->sum = 0.;

    UserCommon->Ndim = nDim;
    tmpCube = (double *)malloc((unsigned int) nDim * sizeof(double));
    UserCommon->err = (double *)malloc((unsigned int) nDim * sizeof(double));
    UserCommon->avg = (double *)malloc((unsigned int) nDim * sizeof(double));

    for ( i = 0; i < nDim; i++ )
    {
        UserCommon->err[i] = 0.;
        UserCommon->avg[i] = 0.;
    }

    // Handle CTRL-C to interrupt the optimisation but not lenstool
    signal(SIGINT, signalReset);
    optInterrupt = 0;

    // Initialize mutex for nchi2 counts
    init_mutex();

    //Run bayesys
    int code = BayeSys3(Common, Objects);
    printf("\nBayeSys3 return code %d\n", code);

    // Destroy mutex
    destroy_mutex();

    // Reset the SIGINT signal to its default behavior
    signal(SIGINT, SIG_DFL);

    fprintf(stdout, "\n");

    /*Set the lens parameters that give the best model
     * and the Gaussian prior at 3sigma limits*/
    o_set_lens_bayes(-4, 3, 3.);

    NPRINTF( stderr, "log(Evidence)          = %lf\n", Common->Evidence );
    NPRINTF( stderr, "Information = -Entropy = %lf\n", Common->Information );
    NPRINTF( stderr, "Calls to chi2() function = %ld\n", UserCommon->Nchi2 );

    /*Free the structures*/
    if( I.stat == 8 && G.nlens != G.nmsgrid )
    {
        free_square_double(zbitstmp1, G.nlens-G.nmsgrid);
        free( ibitstmp );
        free( Common->Data );
        free( Common->Acc );
        free( Objects->Mock );
    }
    free( Objects );
    free( UserCommon->err );
    free( UserCommon->avg );
    free( tmpCube );

    return( Common->Evidence );
}

/* Write the header of the bayes.dat file.
 * Return the number of free parameters/dimensions.
 */
int bayesHeader()
{
    extern struct g_mode    M;
    extern struct g_grille  G;
    extern struct g_image   I;
    extern struct g_pot     P[NPOTFILE];
    extern struct pot       lens[];
    extern struct z_lim     zlim[];
    extern struct z_lim     zalim;
    extern struct galaxie   multi[NFMAX][NIMAX];
    extern int block[][NPAMAX];
    extern int cblock[NPAMAX];
    extern int vfblock[NPAMAX];
    extern struct sigposStr sigposAs;

    int nDim = 0; //# of parameters
    int ipx;      //index of the current parameter in the Param
    int i;        //index of the clump in the lens[] global variable
    int k;        //index of the z_m_limit images to optimise
    int nimages;
    char limages[ZMBOUND][IDSIZE];

    FILE *bayes;
    char name[60];


// Clean the results files
#ifdef DEBUG
    bayes = fopen( "bayes.dbg.dat" , "w" );
#else
    bayes = fopen( "bayes.dat" , "w" );
#endif

    fprintf( bayes, "#Nsample\n");
    //if ( sigposAs.bk != 0 || I.dsigell != -1 ||
    //     (I.stat == 6 && I.statmode == 1) )
        fprintf( bayes, "#ln(Lhood)\n");
    //else
    //    fprintf( bayes, "#Chi2\n");

//  fprintf( bayes, "#Evidence\n");

// Number of lens parameters in the atom
    for ( i = 0; i < G.no_lens; i++ )
    {
        int nDim_old = nDim;
        for ( ipx = CX; ipx <= PMASS; ipx++ )
            if ( block[i][ipx] != 0 )
            {
                nDim++;
                fprintf( bayes, "#%s : %s\n", lens[i].n, getParName(ipx, name, lens[i].type) );
            }
 
        // check against optimized potential without free parameters
        if( nDim_old == nDim )
        {
            fprintf(stderr, "ERROR: Individually optimized potential %s without free parameter.\n", lens[i].n);
            exit(1);
        }
    }

    // multiscale grid clumps
    for ( i = G.nmsgrid; i < G.nlens; i++ )
    {
        fprintf( bayes, "#%s : %s\n", lens[i].n, getParName(B0, name, lens[i].type) );
        if( block[i][B0] != 0 )
            nDim++;
    }

    // source positions of multiple images with grid optimization
    if( G.nmsgrid != G.nlens && nDim < G.nlens - G.nmsgrid )
    {
        for ( i = 0; i < I.n_mult; i++ )
            fprintf( bayes, "#X src %s\n", multi[i][0].n );
        for ( i = 0; i < I.n_mult; i++ )
            fprintf( bayes, "#Y src %s\n", multi[i][0].n );
    }

    // source positions for image reconstruction
    if( M.iclean == 2 )
    {
        extern int sblock[NFMAX][NPAMAX];
        extern struct galaxie source[NFMAX];
        extern struct g_source  S;

        for ( i = 0; i < S.ns; i++ )
            for ( ipx = SCX; ipx <= SFLUX; ipx++ )
                if ( sblock[i][ipx] != 0 )
                {
                    nDim++;
                    fprintf( bayes, "#src %s : %s\n", source[i].n, getParName(ipx, name, 0));
                }
    }

// Number of cosmological parameters in the atom
    for ( ipx = OMEGAM; ipx <= WA; ipx++ )
        if ( cblock[ipx] != 0 )
        {
            nDim++;
            fprintf( bayes, "#%s\n", getParName(ipx, name, 0) );
        }

// Number of redshifts to optimise
    for ( k = 0; k < I.nzlim; k++ )
        if ( zlim[k].bk != 0 )
        {
            nimages = splitzmlimit( zlim[k].n, limages );
            i = 0;
            while ( indexCmp( multi[i][0].n, limages[0] ) ) i++;
            fprintf( bayes, "#Redshift of %s\n", multi[i][0].n );
            nDim++;
        }

// redshift of the arclets
    if ( zalim.bk != 0 )
    {
        fprintf( bayes, "#Redshift of arclets\n");
        nDim++;
    }

// Velocity field parameters to optimise
    for ( ipx = VFCX; ipx <= VFSIGMA; ipx++ )
        if ( vfblock[ipx] != 0 )
        {
            nDim++;
            fprintf( bayes, "#%s\n", getParName(ipx, name, 0) );
        }


// Potfile parameters to optimise
    for ( i = 0; i < G.npot; i++ )
        if ( P[i].ftype != 0 )
        {
            if ( P[i].ircut != 0 )
            {
                nDim++;
                fprintf( bayes, "#Pot%d rcut (arcsec)\n", i);
            }
            if ( P[i].isigma != 0 )
            {
                nDim++;
                fprintf( bayes, "#Pot%d sigma (km/s)\n", i);
            }
            if ( P[i].islope != 0 )
            {
                nDim++;
                if ( P[i].ftype == 62 )
                    fprintf( bayes, "#Pot%d m200slope\n", i);
                else
                    fprintf( bayes, "#Pot%d rcut_slope\n", i);
            }
            if ( P[i].ivdslope != 0 )
            {
                nDim++;
                if ( P[i].ftype == 62 )
                    fprintf( bayes, "#Pot%d c200slope\n", i);
                else
                    fprintf( bayes, "#Pot%d sigma_slope\n", i);
            }
            if ( P[i].ivdscat != 0 )
            {
                nDim++;
                fprintf( bayes, "#Pot%d vdscatter\n", i);
            }
            if ( P[i].ircutscat != 0 )
            {
                nDim++;
                fprintf( bayes, "#Pot%d rcutscatter\n", i);
            }
            if ( P[i].ia != 0 )
            {
                nDim++;
                if ( P[i].ftype == 62 )
                    fprintf( bayes, "#Pot%d m200\n", i);
                else
                    fprintf( bayes, "#Pot%d a\n", i);
            }
            if ( P[i].ib != 0 )
            {
                nDim++;
                if ( P[i].ftype == 62 )
                    fprintf( bayes, "#Pot%d c200\n", i);
                else
                    fprintf( bayes, "#Pot%d b\n", i);
            }
        }

// The noise
    if ( sigposAs.bk != 0 )
    {
        for ( i = 0; i < I.n_mult; i++ )
            for ( k = 0; k < I.mult[i]; k++)
            {
                fprintf( bayes, "#SigposArsec %s (arcsec)\n", multi[i][k].n);
                nDim++;
            }
    }
    if ( I.dsigell != -1. )
    {
        nDim++;
        fprintf( bayes, "#Sigell (arcsec)\n");
    }

    //add Chi2 as last column
    fprintf( bayes, "#Chi2\n");
    
    fclose(bayes);
    return(nDim);
}



static void signalReset(int param)
{
    signal(SIGINT, SIG_DFL);
    optInterrupt = 1;
    fprintf( stderr, "INFO: Optimisation interrupted by CTRL-C\r");
}
