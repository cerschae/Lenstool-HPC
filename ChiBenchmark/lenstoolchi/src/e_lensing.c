#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        e_lensing           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * For each source of the source list, find the list of corresponding
 * arclets and apply a special treatment to those with a tau parameter
 * larger than large_dist.
 *
 * The correspondance between the source and the image planes required
 * that the gimage grid is initialized.
 *
 * Parameters :
 * - source : a list of sources
 * - image[nb_familly][nb_arclets] : a list of arclets
 *
 * Calling sequence : e_lensing() -> e_lens() -> e-test() -> e_mag() -> diag()
 */

void    e_lensing( struct galaxie source[NFMAX],
                   struct galaxie image[NFMAX][NIMAX] )
{
    extern struct g_mode    M;
    const extern struct g_source  S;
    const extern struct g_large   L;
    extern struct point gsource_global[NGGMAX][NGGMAX];
   
    long int    i;

    double dlsds;
    int  ngiant;
    long int ni, nfound, nmissed;
    int    save;

    dlsds = 0;
    ngiant = ni = nfound = nmissed = 0;

    /* classement des sources par z croissant if the sources have different z*/
    i = 1;
    while ( i < S.ns && source[i].z == source[i-1].z ) i++;
    if ( i < S.ns )
        sort(S.ns, source, comparer_z);

    // Switch to silent mode if we create a weak-lensing catalog
    save = M.verbose;
    if ( M.image == 0 && M.source == 0 )
        M.verbose = -1;

    // Initialize the image array

    for ( i = 0 ; i < S.ns ; i++ )
    {
        /*if the source is behind the cluster */
        if ( source[i].dr > PREC_DLSDS )
        {
            /*initialize the gsource grid for the source i*/
            if ( source[i].dr != dlsds )
            {
                dlsds = source[i].dr;
                NPRINTF(stderr, "COMP: grid at z=%.3lf (%.3lf)\n", source[i].z, dlsds);
                e_unlensgrid(gsource_global, dlsds);
            }

            /*find the arclets corresponding to source[i] and fill the
             * image[familly][arclets] list*/
            if ( M.verbose == -1 )
                fprintf( stderr, "Compute multiple images for source %ld...", i);

            ni = e_lens( source[i] , image[i] );
            nfound += ni;
            if ( ni == 0 ) nmissed++;

            if ( M.verbose == -1 )
                fprintf( stderr, "images found %ld, not found %ld\r", nfound, nmissed);

            /*fill the tau parameter of the arclets*/
            e_tau(ni, image[i]);
            sort(ni, image[i], comparer_tau);

            if ( image[i][0].tau > L.dlarge || L.dlarge < 0.1 )
                e_giant( &ngiant, source[i], image[i] );
        }

    }

    if ( nmissed != 0 ) 
	    fprintf( stderr, "\nWARN: There were missed images probably due to the grid resolution.\n");

    M.verbose = save;

    if ( ngiant != 0 )
        w_gianti( ngiant , "gianti.dat" );
}
