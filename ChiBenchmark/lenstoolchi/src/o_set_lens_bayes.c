#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_set_lens_bayes    */
/*      auteur:     Eric Jullo          */
/*      date:       03/05/06            */
/*      place:      ESO, Chile          *
 *
 * Set the lens parameters that give the most probable model.
 *
 * Refer to setBayesModel() for the list of method possibities
 *
 * For each parameter of the <bayes.dat> assign the median value
 */
void o_set_lens_bayes(int method, int prior, double limit)
{


    double   **array;       // contains the bayes.dat data
    int      nParam;
    long int nVal;  //size of array

    // Read the bayes.dat file
    array = readBayesModels(&nParam, &nVal);

    if ( nVal == 0 )
    {
#ifdef DEBUG
        fprintf(stderr, "ERROR: bayes.dbg.dat file not found.\n");
#else
        fprintf(stderr, "ERROR: bayes.dat file not found.\n");
#endif
        exit(-1);
    }

    if ( method >= nVal )
    {
        fprintf(stderr, "ERROR: %d larger than maximum line index %ld\n", method, nVal - 1);
        exit(-1);
    }

    // set the model
    setBayesModel(method, nVal, array);

    // Set the lmin and lmax limits to Gaussian 3 sigma error
    o_set_limit_bayes(nParam, nVal, array, prior, limit);

    free(array);
}

// Set the lmin and lmax values to 1sigma error around the median value
void o_set_limit_bayes(int nParam, long int nVal, double **array,
                       int prior, double limit)
{
    extern struct g_mode   M;
    extern struct g_grille G;
    extern struct g_pot    P[NPOTFILE];
    extern struct g_cosmo  C;
    extern struct g_image  I;
    extern struct pot      lmin[];
    extern struct pot      lmax[];
    extern struct z_lim    zlim[];
    extern struct z_lim    zalim;
    extern int block[][NPAMAX];
    extern int cblock[NPAMAX];
    extern int vfblock[NPAMAX];
    extern struct sigposStr sigposAs;

    int      iParam, ipx;
    long int i;
    double   med, std, err, min, max;

    // Set the position of the first physical parameter in bayes.dat
    iParam = IDPARAM1;

    // Set the optimized clumps limits
    for ( i = 0; i < G.no_lens; i++ )
    {
        for ( ipx = CX; ipx <= PMASS; ipx++ )
            if ( block[i][ipx] != 0 )
            {
                block[i][ipx] = prior;
                med = median(nVal, array[iParam]);
                std = limit * stddev(nVal, array[iParam++]);
                std = std > fabs(med) ? fabs(med) : std;
                err = std / med;

                if ( prior == 3 )
                {
                    min = med;
                    max = std;
                }
                else
                {
                    min = med * ( 1 - err );
                    max = med * ( 1 + err );
                }

                o_set_lmin( i, ipx, min );
                o_set_lmax( i, ipx, max );

                if ( ipx == B0 )
                {
                    lmin[i].sigma = min;
                    lmax[i].sigma = max;
                }
            }
    }

    // for the grid clumps
    for ( i = G.nmsgrid; i < G.nlens; i++ )
    {
        med = median(nVal, array[iParam]);
        std = limit * stddev(nVal, array[iParam++]);
        //std = std > fabs(med) ? fabs(med) : std;
        err = std / med;

        // uniform prior only
        min = med * ( 1 - err );
        max = med * ( 1 + err );

        if( block[i][B0] != 0 )
        {
            lmin[i].sigma = min;
            lmax[i].sigma = max;
        }
        else
        {
            lmin[i].pmass = min;
            lmax[i].pmass = max;
        }
    }

    // Set the source parameters
    if( M.iclean == 2 )
    {
        extern int sblock[NFMAX][NPAMAX];
        extern struct g_source   S;

        for ( i = 0; i < S.ns; i++ )
            for ( ipx = SCX; ipx <= SFLUX; ipx++ )
                if ( sblock[i][ipx] != 0 )
                {
                    sblock[i][ipx] = prior;
                    med = median(nVal, array[iParam]);
                    std = limit * stddev(nVal, array[iParam++]);
                    std = std > fabs(med) ? fabs(med) : std;
                    err = std / med;
    
                    if ( prior == 3 )
                    {
                        min = med;
                        max = std;
                    }
                    else
                    {
                        min = med * ( 1 - err );
                        max = med * ( 1 + err );
                    }
    
                    o_set_lmin(i, ipx, min);
                    o_set_lmax(i, ipx, max);
                }
    }

    // Set the cosmo
    for ( ipx = OMEGAM; ipx <= WA; ipx++ )
        if ( cblock[ipx] != 0 )
        {
            cblock[ipx] = 1;
            med = median(nVal, array[iParam]);
            std = limit * stddev(nVal, array[iParam++]);
            std = std > fabs(med) ? fabs(med) : std;
            err = std / med;
            o_set_lmin(0, ipx, med * ( 1 - err ));
            o_set_lmax(0, ipx, med * ( 1 + err ));
        }

    // rescale the z_m_limit redshifts
    for ( i = 0; i < I.nzlim; i++ )
        if ( zlim[i].opt != 0 )
        {
            zlim[i].bk = prior;
            med = median(nVal, array[iParam]);
            std = limit * stddev(nVal, array[iParam++]);
            std = std > fabs(med) ? fabs(med) : std;
            err = std / med;

            if ( prior == 3 )
            {
                zlim[i].min = med;
                zlim[i].max = std;
            }
            else
            {
                zlim[i].min = med * ( 1 - err );
                zlim[i].max = med * ( 1 + err );
            }

        }

    // rescale the z_a_limit redshift
    if( zalim.bk != 0 )
    {
        zalim.bk = prior;
        med = median(nVal, array[iParam]);
        std = limit * stddev(nVal, array[iParam++]);
        std = std > fabs(med) ? fabs(med) : std;
        err = std / med;

        if ( prior == 3 )
        {
            zalim.min = med;
            zalim.max = std;
        }
        else
        {
            zalim.min = med * ( 1 - err );
            zalim.max = med * ( 1 + err );
        }
    }

    // Set the velocity field
    for ( ipx = VFCX; ipx <= VFSIGMA; ipx++ )
        if ( vfblock[ipx] != 0 )
        {
            vfblock[ipx] = 1;
            med = median(nVal, array[iParam]);
            std = limit * stddev(nVal, array[iParam++]);
            std = std > fabs(med) ? fabs(med) : std;
            err = std / med;
            o_set_lmin(0, ipx, med * ( 1 - err ));
            o_set_lmax(0, ipx, med * ( 1 + err ));
        }

    // Set the potfile limits
    for( i = 0; i < G.npot; i++ )
        if ( P[i].ftype != 0 )
        {
            if ( P[i].ircut != 0 )
            {
                P[i].ircut = prior;
                med = median(nVal, array[iParam]);
                std = limit * stddev(nVal, array[iParam++]);
                std = std > fabs(med) ? fabs(med) : std;
                err = std / med;
                if ( prior == 3 )
                {
                    min = med;
                    max = std;
                }
                else
                {
                    min = med * ( 1 - err );
                    if( min < 0 )   min = 0.;
                    max = med * ( 1 + err );
                }
                P[i].cut1 = min;
                P[i].cut2 = max;
                P[i].cutkpc1 = P[i].cut1 * (d0 / C.h * distcosmo1(P[i].zlens));
                P[i].cutkpc2 = P[i].cut2 * (d0 / C.h * distcosmo1(P[i].zlens));
            }
            if ( P[i].isigma != 0 )
            {
                P[i].isigma = prior;
                med = median(nVal, array[iParam]);
                std = limit * stddev(nVal, array[iParam++]);
                std = std > fabs(med) ? fabs(med) : std;
                err = std / med;
                if ( prior == 3 )
                {
                    min = med;
                    max = std;
                }
                else
                {
                    min = med * ( 1 - err );
                    if( min < 0 )   min = 0.;
                    max = med * ( 1 + err );
                }
                P[i].sigma1 = min;
                P[i].sigma2 = max;
            }
            if ( P[i].islope != 0 )
            {
                P[i].islope = 1;
                med = median(nVal, array[iParam]);
                std = limit * stddev(nVal, array[iParam++]);
                std = std > fabs(med) ? fabs(med) : std;
                err = std / err;
                P[i].slope1 = med * ( 1 - err );
                if( P[i].slope1 < 0 )   P[i].slope1 = 0.;
                P[i].slope2 = med * ( 1 + err );
            }
            if ( P[i].ivdslope != 0 )
            {
                P[i].ivdslope = 1;
                med = median(nVal, array[iParam]);
                std = limit * stddev(nVal, array[iParam++]);
                std = std > fabs(med) ? fabs(med) : std;
                err = std / err;
                P[i].vdslope1 = med * ( 1 - err );
                if( P[i].vdslope1 < 0 )   P[i].vdslope1 = 0.;
                P[i].vdslope2 = med * ( 1 + err );
            }
            if ( P[i].ivdscat != 0 )
            {
                P[i].ivdscat = 1;
                med = median(nVal, array[iParam]);
                std = limit * stddev(nVal, array[iParam++]);
                std = std > fabs(med) ? fabs(med) : std;
                err = std / err;
                P[i].vdscat1 = med * ( 1 - err );
                if( P[i].vdscat1 < 0 )   P[i].vdscat1 = 0.;
                P[i].vdscat2 = med * ( 1 + err );
            }
            if ( P[i].ircutscat != 0 )
            {
                P[i].ircutscat = 1;
                med = median(nVal, array[iParam]);
                std = limit * stddev(nVal, array[iParam++]);
                std = std > fabs(med) ? fabs(med) : std;
                err = std / err;
                P[i].rcutscat1 = med * ( 1 - err );
                if( P[i].rcutscat1 < 0 )   P[i].rcutscat1 = 0.;
                P[i].rcutscat2 = med * ( 1 + err );
            }
        }
    
    // rescale the noise
    if ( sigposAs.bk != 0 )
    {
        med = median(nVal, array[iParam]);
        sigposAs.min = med;
        sigposAs.max = stddev(nVal, array[iParam++]);
    }
    if ( I.dsigell != -1. )
    {
        med = median(nVal, array[iParam++]);
        I.dsigell = med;
    }
}
