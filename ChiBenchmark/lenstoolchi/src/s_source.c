#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>
#include <gsl/gsl_rng.h>

static void e_magshear(struct galaxie *arclet, double *ka, double *ga1, double *ga2);

/****************************************************************/
/*      nom:        s_source            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/* Create a list of sources from
 * - randomly with the tirer() function
 * - the source file M.sourfile
 * - the images file M.imafile
 *
 * Fill the global lists <source> and <arclet> and write the sources
 * in the file source.dat
 *
 ****************************************************************/

void    s_source()
{
    extern struct g_mode    M;
    extern struct g_image I;
    extern struct z_lim zlim[];
    int k, l, nimages;
    char limages[ZMBOUND][IDSIZE];
    extern struct g_source  S;
    extern struct galaxie   arclet[NAMAX], source[NFMAX];
    extern struct galaxie  multi[NFMAX][NIMAX];
    extern struct pot       lens[];
    long int   i, na;
    double  ka, ga1, ga2;
    gsl_ran_discrete_t *gsmail;
    gsl_rng *seed = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(seed, S.rand); 
 
    // in case Smail et al. distibution is used
    if( S.distz == 2 && S.par1 != 0 )
        gsmail = smailpreproc();

    NPRINTF(stderr, "SET: sources\n");
    if ( M.image == 0 && M.source == 0 )
        tirer(source);  /*optimisation dans le plan source*/
    else
    {
        if (S.ns > 0)
            NPRINTF(stderr, "[WARNING] Source index is already set to %ld\n", S.ns);

        if (M.source != 0)
        {
            if( M.source > 0 )  
                f_shape(&S.ns, source, M.sourfile, M.source);

            if( S.ns > NFMAX )
            {
                fprintf(stderr, "[ERROR] Too many sources in %s. NFMAX limit is %d\n", M.sourfile, NFMAX);
                exit(1);
            }

            for ( i = 0; i < S.ns; i++ )
            {
                if( source[i].z == 0 )
                {
                    if( S.distz == 0 )
                        source[i].z = S.zs;
                    else if( S.distz == 1 )
                        source[i].z = gsl_ran_flat(seed, S.zsmin, S.zsmax);
                    else if( S.distz == 2 && S.par1 != 0 )
                        source[i].z = d_rndzsmail(seed, gsmail);
//                    else
//                        source[i].z = d_rndz(S.zsmax, &seed); // TODO: update with GSL

                    NPRINTF(stderr, "WARN: source %s redshift set to z=%lf\n", source[i].n, source[i].z);
                }

                dratio_gal(&source[i], lens[0].z);
            }
        }

//      else if (M.source==2)
//          f_shape2(&S.ns,source,M.sourfile);

        if (M.image != 0)
        {
            na = 0;
            if( M.image > 0 )
                f_shape(&na, arclet, M.imafile, M.image);  

            if ( na > NFMAX )
            {
                printf("ERROR: Too many arclets in %s. NFMAX limit is %d\n", M.imafile, NFMAX );
                exit(1);
            }

            NPRINTF(stderr, "COMP: arclets shear components\n");

            pro_arclet(na, arclet);
            for (i = 0; i < na; i++)
            {
                // If images have no redshift, look up into multi[][]
                if( arclet[i].z == 0. )
                {
                    for( k = 0; k < I.n_mult && arclet[i].z == 0.; k++ )
                        if( !indexCmp(arclet[i].n, multi[k][0].n) )
                            arclet[i].z = multi[k][0].z;
                }

                // If still no redshift is found then look into z_m_limit statements
                if( arclet[i].z == 0. )
                {
                    for( k = 0; k < I.nzlim && arclet[i].z == 0.; k++ )
                    {
                        nimages = splitzmlimit(zlim[k].n, limages);
                        for( l = 0; l < nimages; l++ )
                            if( !indexCmp(limages[l], arclet[i].n) )
                            {
                                arclet[i].z = zlim[k].min;
                                NPRINTF(stderr, "INFO : Set redshift of %s to %lf\n",arclet[i].n, arclet[i].z);
                            }
                    }
                }

                // If still no redshift is found then assign one from source section
                if( arclet[i].z == 0 )
                {
                    if( S.distz == 0 )
                        arclet[i].z = S.zs;
                    else if( S.distz == 1 )
                        arclet[i].z = gsl_ran_flat(seed, S.zsmin, S.zsmax);
                    else if( S.distz == 2 && S.par1 != 0 )
                        arclet[i].z = d_rndzsmail(seed, gsmail);
//                    else
//                        arclet[i].z = d_rndz(S.zsmax, &seed);  //TODO: Update with GSL

                    NPRINTF(stderr, "WARN: arclet %s redshift set to z=%lf\n", arclet[i].n, arclet[i].z);
                }

                dratio_gal(&arclet[i], lens[0].z);
                e_magshear(&arclet[i], &ka, &ga1, &ga2);
                arclet[i].kappa = ka;
                arclet[i].gamma1 = ga1;
                arclet[i].gamma2 = ga2;
            }
            /*      ecrire_r(0,na,arclet,"aletshear.dat"); */

            e_unlens(na, arclet, &S.ns, source);
        }
    }

    pro_arclet(S.ns, source);

    ecrire_r(0, S.ns, source, "source.dat", 2); // relative coordinates

    // free gsl random variables
    gsl_rng_free(seed);
    if( S.distz == 2 && S.par1 != 0 )
        gsl_ran_discrete_free(gsmail);

}


static void e_magshear(struct galaxie *arclet, double *ka, double *ga1, double *ga2)
{
    struct matrix M;

    M = e_grad2_gal(arclet, NULL);
    M.a /= arclet->dos;
    M.b /= arclet->dos;
    M.c /= arclet->dos;
    M.d /= arclet->dos;
    *ka = (M.a + M.c) / 2.;
    *ga1 = (M.a - M.c) / 2.;
    *ga2 = M.b;
}
