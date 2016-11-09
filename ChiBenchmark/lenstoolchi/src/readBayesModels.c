#include<stdlib.h>
#include<float.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include"lt.h"

/****************************************************************/
/*      nom:        readBayesModels */
/*      auteur:     Eric Jullo          */
/*      date:       07/08/06            */
/*      place:      ESO, Chile          *
 *
 * Read the bayes.dat file and return an array[nParam][val].
 * where nParam = [chi2 -> Last param]
 *
 * It also analyse the header of bayes.dat and set
 * the <G.no_lens>, <G.nlens>, <block>, <cblock> and <P.ircut>,
 * <P.isigma> global variables to their corresponding values.
 *
 * Return:
 * - the bayes.dat array and nParam and nVal the size of the array.
 * - 0 for nVal and nParam if bayes.dat cannot be read.
 *
 * Warning :
 *  - Theta in <bayes.dat> is in degree but is converted to radian in <array>
 *  - sigma in <bayes.dat> is in km/s but is converted according to the lens type
 *    in <array>.
 *
 * Then, the user can generate an error mass map even if he doesn't perform
 * an optimisation.
 */

static double lhood_mode; // If 1, bayes.dat contains #ln(Lhood), otherwise 0

static int analyseHeader( char *line );
static double bayes2lt(long int ilens, int ipx, double val);
static void convertArray(long int nVal, double **array );
static double method(long int n, int param, double **array);
static double  **readBayesASCII(int *nParam, long int *nVal);
static double  **readBayesFITS(int *nParam, long int *nVal);
static void reset_grad_multi();
static void reset_grad_arclet();

double  **readBayesModels(int *nParam, long int *nVal)
{
    extern struct g_grille  G;
    extern struct g_image   I;
    extern struct g_pot     P[NPOTFILE];
    extern struct g_source  S;

    extern int block[][NPAMAX];
    extern int cblock[NPAMAX];
    extern int sblock[NFMAX][NPAMAX];
    extern struct sigposStr sigposAs;

    double **array;
    char   line[LINESIZE];
    long int    nLines;
    long int ilens;
    int    i, ipx;

    // Initialize ouput variables
    *nParam = 0;
    *nVal = 0;

    //
    // REINITIALIZE ALL VARIABLES
    //

    // Get G.no_lens from bayes.dat
    // Reset G and block global variables
    G.no_lens = 0;
    for ( ilens = 0; ilens < G.nlens; ilens++)
        for ( ipx = CX; ipx <= PMASS; ipx ++)
            block[ilens][ipx] = 0;

    // Source parameters
//    for ( i = 0; i < S.ns; i++ )
//        for( ipx = SCX; ipx <= SFLUX; ipx++ )
//            sblock[i][ipx] = 0;

    // Reset the cblock global variables
    for ( ipx = OMEGAM; ipx <= WA; ipx++)
        cblock[ipx] = 0;
 
    // Reset the redshift I.nzlim
//    I.nzlim = 0;

    // Reset the potfile
    for( i = 0; i < G.npot; i++ )
    {
        P[i].ircut = 0;
        P[i].isigma = 0;
        P[i].islope = 0;
        P[i].ivdslope = 0;
        P[i].ivdscat = 0;
        P[i].ircutscat = 0;
    }

    // Reset the nuisance
    sigposAs.bk = 0;
    I.dsigell = -1.;

    // Read bayes.dat
    lhood_mode = 0;
    array = readBayesASCII(nParam, nVal);

    if ( array == NULL )
        return 0;

    // Convert <array> to values directly usable in lenstool
    convertArray(*nVal, array);

    return array;
}

static int imethod;     // use a specified method as parameter estimator
// -1 : mean value
// -2 : median value
// -3 : mode value
// -4 : best chi2 value
// -5 : to be used by o_rescale()

// Reset the parameters of a potfile
// Return 1 if the potfile parameters have been updated
int setpotBayes(struct g_pot *pot, int *iParam_in, long int nVal, double **array)
{
    int iParam = *iParam_in;
    int computePot = 0;

    if ( pot->ircut != 0 )
    {
        pot->cut = method(nVal, iParam++, array);
        computePot = 1;
    }
    if ( pot->isigma != 0 )
    {
        pot->sigma = method(nVal, iParam++, array);
        computePot = 1;
    }
    if ( pot->islope != 0 )
    {
        pot->slope = method(nVal, iParam++, array);
        computePot = 1;
    }
    if ( pot->ivdslope != 0 )
    {
        pot->vdslope = method(nVal, iParam++, array);
        computePot = 1;
    }
    if ( pot->ivdscat != 0 )
    {
        pot->vdscat = method(nVal, iParam++, array);
        computePot = 1;
    }
    if ( pot->ircutscat != 0 )
    {
        pot->rcutscat = method(nVal, iParam++, array);
        computePot = 1;
    }
    if ( pot->ia != 0 )
    {
        pot->a = method(nVal, iParam++, array);
        computePot = 1;
    }
    if ( pot->ib != 0 )
    {
        pot->b = method(nVal, iParam++, array);
        computePot = 1;
    }

    *iParam_in = iParam;
    return computePot;
}

/* Set all the lenstool parameters to the values contained in the
 * bayes.dat file at line iVal.
 * If iVal < 0 then iVal specifies the method corresponding to the imethod index.
 */
void    setBayesModel( long int iVal, long int nVal, double **array)
{
    extern struct g_mode    M;
    extern struct g_grille  G;
    extern struct g_cosmo   C;
    extern struct g_image   I;
    extern struct g_pot P[NPOTFILE];

    extern struct galaxie   multi[NFMAX][NIMAX];
    extern struct galaxie   arclet[NAMAX];
    extern struct z_lim zlim[];
    extern struct z_lim zalim;
    extern struct cline     cl[];
    extern struct pot       lens[];
    extern int  block[][NPAMAX];
    extern int  cblock[NPAMAX];
    extern int  vfblock[NPAMAX];
    extern struct sigposStr sigposAs;
    extern double *np_b0;
    extern long int narclet;

    int     computedr_m, computedr_a, computePot[NPOTFILE];  // flags
    int     nimages;
    char    limages[ZMBOUND][IDSIZE];
    double  tmp;
    double  dl0s, dos;
    long int ilens, i;
    int     iParam, ipx;
    int     j, k, l;

    imethod = iVal;
    if ( imethod == -5 )
        iParam = 0;   // assign cube[0][iParam] to the model parameters
    else
        iParam = IDPARAM1;

    computedr_m = computedr_a = 0;    // do not recompute dlsds ratio for multiple images
    for( i = 0; i < G.npot; i++ )
        computePot[i] = 0;   // potfile parameters are not modified

    for ( ilens = 0; ilens < G.no_lens; ilens++ )
    {
        for ( ipx = CX; ipx <= PMASS; ipx++ )
            if ( block[ilens][ipx] != 0 )
                o_set_lens( ilens, ipx, method(nVal, iParam++, array) );

        // TODO in o_set_lens() assign sigma directly to lens.sigma rather than lens.b0
        if ( block[ilens][B0] != 0 )
            lens[ilens].sigma = lens[ilens].b0;

        if ( block[ilens][ZLENS] != 0 )
            computedr_m = computedr_a = 1;

        set_dynamics(ilens);
    }

    // Compute b0 from msgrid potentials either from sigma or rhos random values
//    if ( G.nmsgrid < G.nlens )
//    {
//        if ( block[G.nmsgrid][B0] != 0 )
//            for ( ilens = G.nmsgrid; ilens < G.nlens; ilens++ )
//            {
//                tmp = method(nVal, iParam++, array);
//                tmp = 6.*pia_c2 * tmp * tmp;
//                np_b0[ilens-G.nmsgrid] = tmp;
//                lens[ilens].b0 = tmp;   
//            }
//        else
//        {
//            double *atmp = (double *) malloc((G.nlens - G.nmsgrid ) * sizeof(double));
//            for ( ilens = G.nmsgrid; ilens < G.nlens; ilens++ )
//                atmp[ilens - G.nmsgrid] = method(nVal, iParam++, array);
//        
//            setgrid_rhos2b0(atmp);
//            free(atmp);
//        }
//    }
    for ( ilens = G.nmsgrid; ilens < G.nlens; ilens++ )
        np_b0[ilens - G.nmsgrid] = lens[ilens].b0 =  method(nVal, iParam++, array);

    // Set Multiscale grid SL source positions
    if ( G.nmsgrid != G.nlens && I.n_mult != 0)
    {
        extern struct galaxie source[NFMAX];
        for ( i = 0; i < I.n_mult; i++ )
            source[i].C.x = method(nVal, iParam++, array);

        for ( i = 0; i < I.n_mult; i++ )
            source[i].C.y = method(nVal, iParam++, array);
    }
 
// set the source positions
    if( M.iclean == 2 )
    {
        extern struct g_source  S;
        extern struct galaxie source[NFMAX];
        extern int  sblock[NFMAX][NPAMAX];
        char str[IDSIZE], *pch;

        for ( i = 0; i < S.ns; i++ )
        {
            // if source is attached to an image, reset its coordinates
            strcpy(str, source[i].n);
            pch = strtok(str, "_");
            pch = strtok(NULL, "_");
            if( pch != NULL )
                source[i].C.x = source[i].C.y = 0.;

            for( ipx = SCX ; ipx <= SFLUX; ipx++ )
                if( sblock[i][ipx] != 0 )
                    o_set_source( &source[i], ipx, method(nVal, iParam++, array));

            if( sblock[i][SEPS] != 0 )
                source[i].E.b = source[i].E.a * (1. - source[i].eps) / (1. + source[i].eps);
       }
    }

    // Set the cosmological parameters
    for ( ipx = OMEGAM; ipx <= WA; ipx++ )
        if ( cblock[ipx] != 0 )
        {
            o_set_lens( 0, ipx,  method(nVal, iParam++, array) );
            computedr_m = computedr_a = 1;
        }

    // Set the optimized redshifts to the corresponding images
    for ( i = 0 ; i < I.nzlim; i++ )
        if ( zlim[i].bk != 0 )
        {
            tmp = method(nVal, iParam++, array);
            nimages = splitzmlimit(zlim[i].n, limages);
            // for all the families in zlim[ipx].n
            for ( j = 0; j < nimages; j++ )
            {
                if( ! strncmp(limages[j], "cl", 2) )
                {
                    // loop over the critical lines
                    int zmid;
                    sscanf(limages[j], "cl%d", &zmid);
                    cl[zmid].z = tmp;
                }
                else
                {
                    // loop over the images families
                    k = 0;
                    while ( indexCmp( multi[k][0].n, limages[j] ) && k < I.n_mult ) k++;
                    // set the reshift of all the images of the family
                    for ( l = 0; l < I.mult[k]; l++ )
                        multi[k][l].z = tmp;
                }
            }

            computedr_m = 1;
        }
 
    // Set the optimized redshifts of the arclets with unknown redshift
    if ( zalim.bk != 0 )
    {
        I.zarclet = method(nVal, iParam++, array);
        for( j = 0; j < zalim.opt; j++ )
            arclet[j].z = I.zarclet;

        computedr_a = 1;
    }

    // Set the velocity field parameters
    for ( ipx = VFCX; ipx <= VFSIGMA; ipx++ )
        if ( vfblock[ipx] != 0 )
        {
            o_set_lens( 0, ipx,  method(nVal, iParam++, array) );
        }

    //rescale the pot parameters
    for( i = 0; i < G.npot; i++ )
        if (P[i].ftype != 0 )
            computePot[i] = setpotBayes(&P[i], &iParam, nVal, array);

    //set the nuisance parameters
    if ( sigposAs.bk != 0. )
    {
        for ( i = 0; i < I.n_mult; i++ )
            for ( j = 0; j < I.mult[i]; j++ )
            {
                tmp = method(nVal, iParam++, array);
                I.sig2pos[i][j] = tmp * tmp;
            }
    }

    if ( I.dsigell != -1. )
    {
        I.sigell = method(nVal, iParam++, array);
        I.sig2ell = I.sigell * I.sigell;
    }

    // recompute the dlsds ratios if necessary
    if ( computedr_a && computedr_m )
    {
        if ( C.kcourb == 0. )
            C.omegaX = 1 - C.omegaM;
        else
            C.kcourb = C.omegaM + C.omegaX - 1;

        // ...same for the critical line dlsds
        for ( i = 0 ; i < I.npcl; i++ )
        {
            if ( lens[0].z < cl[i].z )
                cl[i].dl0s = distcosmo2(lens[0].z, cl[i].z);
            else
                cl[i].dl0s = 0.;
            cl[i].dos = distcosmo1(cl[i].z);
            cl[i].dlsds = cl[i].dl0s / cl[i].dos;
        }
    }

    if ( computedr_m )
    {
        // set the DLS/DS ratio for all the images according
        // to the new cosmological parameters
        for ( i = 0 ; i < I.n_mult; i++ )
            for ( j = 0; j < I.mult[i]; j++ )
            {
                dos = distcosmo1(multi[i][0].z);

                if ( lens[0].z < multi[i][0].z )
                    dl0s = distcosmo2(lens[0].z, multi[i][0].z);
                else
                    dl0s = 0.;

                for ( j = 0 ; j < I.mult[i]; j++ )
                {
                    multi[i][j].dos = dos;
                    multi[i][j].dl0s = dl0s;
                    multi[i][j].dr = dl0s / dos;
                }
            }
    }

    //...the same for the arclets
    if( computedr_a )
    {
        //... arclets with optimized redshifts
        if( zalim.bk != 0 )
        {
            dos = distcosmo1(I.zarclet);
            if ( lens[0].z < I.zarclet )
                dl0s = distcosmo2(lens[0].z, I.zarclet);
            else
                dl0s = 0.;
    
            for ( i = 0 ; i < zalim.opt; i++ )
            {
                arclet[i].dos = dos;
                arclet[i].dl0s = dl0s;
                arclet[i].dr = dl0s / dos;
            }
        }

        // ... the rest of the arclets
        if( computedr_m )
            for( i = zalim.opt; i < narclet; i++ )
                dratio_gal(&arclet[i], lens[0].z);
    }

// recompute the potfile if necessary
    for( i = 0; i < G.npot; i++ )
        if ( computePot[i] )
        {
            scale_pot(&P[i]);
            computedr_a = computedr_m = 1;
        }

// Reset the image->grad and image->grad2 accelerators
    if( computedr_m )
        reset_grad_multi();

    if( computedr_a )
        reset_grad_arclet();
}

// Read the header of bayes.dat and call analyseHeader()
// to process each line
static double ** readBayesASCII(int *nParam, long int *nVal)
{
    extern struct g_grille  G;
    long int    nLines;
    FILE     *bayes;
    int    error;
    char   line[LINESIZE];
    char   *pch;    // token to split a line
    double **array;
    int j;

    // Read bayes.dat
#ifdef DEBUG
    bayes = fopen( "bayes.dbg.dat", "r");
#else
    bayes = fopen( "bayes.dat" , "r" );
#endif

    if ( bayes == NULL )
    {
        //fprintf(stderr, "WARNING: Reading bayes.dat. No such file found.\n");
        return 0;
    }

    // Read header
    error = 0;
    nLines = wc(bayes);
    fgets( line, LINESIZE, bayes);
    line[strlen(line)-1] = '\0';  // replace ending \n by \0
    while ( !feof(bayes) && !ferror(bayes) && line[0] == '#' && !error )
    {
        (*nParam)++;
        error = analyseHeader( line );
        fgets( line, LINESIZE, bayes);
        line[strlen(line)-1] = '\0'; 
    }
    
    // Deal with error code
    if ( error )
    {
        fclose(bayes);
        return 0;
    }
 
    // Reinitialise G.no_lens in case part of the clumps are accelerated
    if ( G.no_lens > G.nmsgrid )
        G.no_lens = G.nmsgrid;

    // Initialize array
    *nParam = (*nParam) - 1; //because Nsample is excluded from array
    array = alloc_square_double( *nParam , nLines - (*nParam) - 1);

    // Read the data (first line comes from header analysis)
    nLines = *nParam + 1; 
    while ( !feof(bayes) && !ferror(bayes) )
    {
        // Test for line too long
        if ( strlen(line) + 1 == LINESIZE )
        {
            fprintf(stderr, "ERROR: while reading bayes.dat line %ld. Increase LINESIZE(current value is=%d) to match the line size in bayes.dat.\n", nLines, LINESIZE);
            exit(1);
        }

        // Skip the fist column and read the 2nd one
        pch = strtok( line, " " );

        // Fill array[][*nVal]
        j = 0; 
        while ( j < *nParam && pch != NULL )
        {
            pch = strtok (NULL, " ");
            if ( pch != NULL )
            {
                int lnval;
                double val;
                if( sscanf( pch, "%dx%lf", &lnval, &val ) == 2 )
                {
                    long int i;
                    for ( i = 0; i < lnval; i++ )
                        array[j++][*nVal] = val;
                }
                else if( sscanf( pch, "%lf", &array[j++][*nVal] ) != 1 )
                {
                    fprintf( stderr, "ERROR: while reading bayes.dat line %ld, column %d. pch = %s\n", *nVal, j + 1, pch);
                    pch = NULL;
                    *nVal -= 1;
                }
            }
        }
        // normal exit: (j == *nParam && pch != NULL) || pch == NULL (it points on the last element)

        if( j < *nParam )  // if true pch != NULL
        {
            fprintf( stderr, "ERROR: while reading bayes.dat line %ld, not enough parameters (%d/%d).\n", nLines, j+1, *nParam+1);
            *nVal -= 1;
        }

        fgets(line, LINESIZE, bayes);
        line[strlen(line)-1] = '\0';

        nLines ++;

        (*nVal) ++;
    }
    fprintf( stderr, "INFO: %ld good lines found in bayes.dat\n", *nVal);

    // Catch errors
    if ( ferror(bayes) || ( feof(bayes) && *nVal == 0 ) )
    {
//      fprintf( stderr, "WARNING : reading bayes.dat\n" );
        *nVal = *nParam = 0;
        fclose(bayes);
        return array;
    }
    fclose(bayes);

    return array;

}

// Read the header of bayes.fits and call analyseHeader()
// to process each line
static double ** readBayesFITS(int *nParam, long int *nVal)
{
}

/* Analyse a line of the bayes.dat header and eventually
 * In the case of reading a best.par file, reset the
 * <G.no_lens>, <block>, <cblock>, <zlim> and <P.ircut>, <P.isigma>
 * global variables as if they were optimized
 */
static int analyseHeader( char *line )
{
    extern struct g_grille  G;
    extern struct g_image   I;

    extern struct pot       lens[];
    extern struct z_lim     zlim[];
    extern int block[][NPAMAX];
    extern int cblock[NPAMAX];
    extern struct sigposStr sigposAs;

    long int   i;
    int ipx;
    char  *pch, *pch2, name[50];
    int updated;    // flag to prevent the case when the clump is already set as optimised

    // skip the '#'
    line = line + 1;

    if ( !strcmp(line, "Nsample")  ||
         !strcmp(line, "Evidence") ||
         !strcmp(line, "Chi2") )
        return 0; // not an error

    if ( !strcmp(line, "ln(Lhood)") )
    {
        lhood_mode = 1;
        return 0; // not an error
    }


    // Get the first word
    pch = strtok( line, " " );

    // Skip the "O" if it exists in the <bayes.dat> file... <pch> in bayes.dat
    // "O" means a lens[] element
    if ( pch[0] == 'O' )
    {
        pch++;

        // Look for the right lens to change
        updated = 0;
        for ( i = 0; i < G.nlens && updated == 0 ; i++ )
        {
            // ... and in Lenstool... <pch2> in Lenstool
            pch2 = lens[i].n;
            if ( pch2[0] == 'O' )
                pch2++;

            if ( !strcmp( pch, pch2 ) )
            {
                pch = strtok( NULL, " " );  // skip the ':'
                pch = strtok( NULL, " " );
                // Look for the right parameter to change
                for ( ipx = CX; ipx <= PMASS; ipx++ )
                {
                    pch2 = strtok( getParName(ipx, name, lens[i].type), " ");   // remove the unit if any
                    if ( !strcmp( pch, pch2) && block[i][ipx] == 0  )
                        block[i][ipx] = 1;
                }
                // Prevent the case when the clump is already set as optimised
                for ( ipx = CX; ipx <= PMASS; ipx++ )
                    updated += block[i][ipx];

                // Increment G.no_lens only once per clump
                if ( updated == 1 )
                {
                    G.no_lens ++;
                }
                else if ( updated == 0 )
                {
                    fprintf( stderr, "ERROR: bayes.dat header format unrecognised\n");
                    exit(-1);
                }
            }
        }
    }

    // Look for a source parameter
    // TODO: Not finished. Find a way to check the number of sources
    if ( !strcmp(pch, "src") )
    {
        extern struct g_source  S;
        extern struct galaxie source[NFMAX];
        pch = strtok( NULL, " " ); 
        i = 0;
        while( strcmp(pch, source[i].n) && i < S.ns ) i++;
        
        if ( i == S.ns && strcmp(pch, source[i - 1].n) )
        {
            fprintf(stderr, "ERROR: number of sources in lenstool (%ld) and found in bayes.dat mismatch\n", S.ns);
            exit(-1);
        }
    }
    
    // Look for cosmological parameters
    if ( !strcmp(pch, "omegaM") )
        cblock[OMEGAM] = 1;
    if ( !strcmp(pch, "omegaX") )
        cblock[OMEGAX] = 1;
    if ( !strcmp(pch, "wX") )
        cblock[WX] = 1;
    if ( !strcmp(pch, "wa") )
        cblock[WA] = 1;

    // Look for the redshifts
    // THIS HAS TO BE DEFINED IN THE .PAR FILE
//    if ( !strcmp(pch, "Redshift") )
//    {
//        pch = strtok( NULL, " " ); // pch contains 'of'
//        pch += 3;   // pch points on the last part of line
//        strcpy(zlim[I.nzlim].n, pch);
//        // Append an image to optimize
//        zlim[I.nzlim].opt = 1;
//        COMMENTED EJ20100306 : trust .par input file
//        I.nzlim++;
//    }

    // Look for the Potfile keyword "Pot" or "Pot1" or "Pot2" etc...
    int ipot = 0;
    if ( !strcmp(pch, "Pot") || sscanf(pch, "Pot%d", &ipot) == 1 )
    {
        extern struct g_pot     P[NPOTFILE];
        struct g_pot *pot = &P[ipot];  // default with P[0]

//        if ( pot->core == -1 )
//        {
//            fprintf( stderr, "WARNING : Potfile parameters in bayes.dat not found in parfile.\n");
//            return 1;
//        }
        pch = strtok( NULL, " " ); // pch contains rcut or sigma
        if ( !strcmp(pch, "rcut") )
            pot->ircut  = 1;
        if ( !strcmp(pch, "sigma") )
            pot->isigma  = 1;
        if ( !strcmp(pch, "rcut_slope") || !strcmp(pch, "m200slope") )
            pot->islope = 1;
        if ( !strcmp(pch, "sigma_slope") || !strcmp(pch, "c200slope") )
            pot->ivdslope = 1;
        if ( !strcmp(pch, "vdscatter") )
            pot->ivdscat = 1;
        if ( !strcmp(pch, "rcutscatter") )
            pot->ircutscat = 1;
        if ( !strcmp(pch, "a") || !strcmp(pch, "m200") )
            pot->ia = 1;
        if ( !strcmp(pch, "b") || !strcmp(pch, "c200") )
            pot->ib = 1;
    }

    // Look for the nuisance
    if ( !strcmp(pch, "SigposArcsec") )
        sigposAs.bk = 1;   // this is not the correct sigma but it doesn't matter
    // since the values are written in the bayes.dat file

    if ( !strcmp(pch, "Sigell") )
        I.dsigell = 1.;

    return 0;
}

/* Convert input <bayes.dat> values to the corresponding values readily usable
 * in lenstool.
 */
static double bayes2lt(long int i, int ipx, double val)
{
    extern struct g_cosmo C;
    extern struct pot       lens[];

    switch (ipx)
    {
        case(THETA):
            val *= DTR;
            break;
        case(B0):
            if ( lens[i].type == 13 )
                val *= 1e8;
            break;
            /* now it is in solar mass in bayes.dat also
                    case(PMASS):
                        if( lens[i].type == 12 )
                            val *= 1e14; // rhos in 10^14 Msol -> Msol
                        break;
                    case(MASSE):
                        if( lens[i].type == 12 )
                            val *= 1e14; // masse in 10^14 Msol -> Msol
                        break;
            */
        default:
            break;
    }

    return val;
}

static void convertArray(long int nVal, double **array )
{
    extern struct g_mode     M;
    extern struct g_grille   G;
    extern struct g_source   S;
    extern int sblock[NFMAX][NPAMAX];
    extern int block[][NPAMAX];

    long int   iVal, i;
    int   nParam, ipx;

    for ( iVal = 0; iVal < nVal; iVal++ )
    {
        nParam = IDPARAM1; // position of the first parameter in array
        // 1st column in array is chi2
        if ( lhood_mode )    array[0][iVal] *= -1.; // so that the best lhood has the lowest value
        // ie same behavior as chi2

        // Optimized clumps
        for ( i = 0; i < G.no_lens; i++ )
            for ( ipx = CX; ipx <= PMASS; ipx++ )
                if ( block[i][ipx] != 0 )
                {
                    array[nParam][iVal] = bayes2lt( i, ipx, array[nParam][iVal] );
                    nParam++;
                }
        
        // Grid clumps
        for ( i = G.nmsgrid; i < G.nlens; i++ )
            nParam++;

        // Sourcelimit parameters
        if( M.iclean == 2 )
            for ( i = 0; i < S.ns; i++ )
                for( ipx = SCX; ipx <= SFLUX; ipx ++ )
                    if( sblock[i][ipx] != 0 )
                    {
                        if( ipx == STHETA )
                            array[nParam][iVal] = array[nParam][iVal] * DTR;
    
                        nParam ++;
                    }
    }
}

/* Return the value in <list> that has produced the best chi2*/
static double bestchi2(long int n, double *chi2, double *list)
{
    double chi2min;
    long int i, imin;

    imin = 0;
    chi2min = DBL_MAX;

    for ( i = 0; i < n; i++)
        if ( chi2[i] < chi2min )
        {
            chi2min = chi2[i];
            imin = i;
        }

    return( list[imin] );
}

/* Return an estimated value of the values in <array> depending
 * on the value of the <imethod> value.
 * <imethod> may also be a line number in array[param].
 * n is the number of elements in array[param].
 */
static double method(long int n, int param, double **array)
{
    double val;

    switch ( imethod )
    {
        case(-1):
            val = mean(n, array[param]);
            break;
        case(-2):
            val = median(n, array[param]);
            break;
        case(-3):
            val = mode(n, array[param]);
            break;
        case(-4):
            val = bestchi2(n, array[0], array[param]);
            break;
        case(-5):   // o_rescale Atoms, with n parameters per atom 
            val = array[0][param];
            break;
        default:
            val = array[param][imethod];
            break;
    }

    return val;
}

/*
 * reset the multi.grad(2) parameters for the optimisation.
 */
static void reset_grad_multi()
{
    extern struct g_image   I;
    extern struct galaxie   multi[][NIMAX];

    long int i;
    int k;

    // reset the grad and grad2 arclet temporary structure for optimisation
    for ( i = 0; i < I.n_mult; i++ )
        for ( k = 0; k < I.mult[i]; k++ )
        {
            multi[i][k].grad2.a = 0;
            multi[i][k].grad2.c = 0;
            multi[i][k].grad.x = 0;
            multi[i][k].grad.y = 0;
        }
}

/*
 * reset the arclet.grad(2) parameters for the optimisation.
 */
static void reset_grad_arclet()
{
    extern struct galaxie   arclet[NAMAX];
    extern long int narclet;

    long int i;

    // reset the grad and grad2 arclet temporary structure for optimisation
    for ( i = 0 ; i < narclet ; i++ )
    {
        arclet[i].grad2.a = arclet[i].grad2.c = 0;
        arclet[i].grad.x = arclet[i].grad.y = 0;
    }
}
