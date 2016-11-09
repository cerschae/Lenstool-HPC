#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<float.h>
#include "fonction.h"
#include "constant.h"
#include "dimension.h"
#include "structure.h"

/****************************************************************/
/*                 Program  : grille                */
/*                 Version  : 1 mai 1992            */
/*                 Location : Obs. Toulouse         */
/*                 Auteur   : jean-paul             */
/*  EJ 30/11/05 Modify indentation and add void flire           */
/*  EJ 2/12/05  Rearrange the coef in rcut/sigma... equations*/
/****************************************************************/

static void convertXY_compat(double *x, double *y, int isAbsWCS);

int set_potfile(int nplens, struct g_pot *pot)
{
    extern struct g_mode   M;
    extern struct g_grille G;
    extern struct pot        lens[];
    struct pot *ilens;
    double  aa, bb;
    char      line[256];
    long int  i;
    FILE      *POT;

    double  ref_ra, ref_dec;
    int     iref;

    iref = 0;   // Default, absolute wcs coordinates
    i = nplens;
    ilens = &lens[i];

    if ( pot->ftype == 0 )
        return(0);

    // Check redshift for ftype 1 or 3
    if( pot->ftype == 1 || pot->ftype == 3 )
        if( pot->zlens <= 0. )
        {
            NPRINTF(stderr, "ERROR: z_lens keyword not defined in the potfile section\n");
            exit(-1);
        }

    //********************************************************************
    //  Open and read  the potfile
    //********************************************************************
    NPRINTF(stderr, "READ: %s nlens:%d ftype:%d\n", pot->potfile, nplens, pot->ftype);

    POT = fopen(pot->potfile, "r");

    if ( POT == NULL || ferror(POT) )
    {
        fprintf( stderr, "ERROR: Error reading %s\n", pot->potfile);
        exit(-1);
    }

    // Read the potfile
    fgets(line, 128, POT);
    while ( !feof(POT) && !ferror(POT) && i < NLMAX )
    {
        // Check for a WCS header
        if ( strstr(line, "#REFERENCE" ) != NULL )
        {
            getRADEC( line, &iref, &ref_ra, &ref_dec );
            goto NEXT;
        }

        // Skip commented lines with #
        if ( line[0] == '#' )
            goto NEXT;

        // Default initialisation of clump ilens
        ilens->type = pot->type;
        ilens->z = pot->zlens;
        ilens->epot = ilens->emass = 0.;
        ilens->alpha = ilens->beta = 0.;
        ilens->rckpc = ilens->rc = 0.;
        ilens->rcutkpc = ilens->rcut = DBL_MAX;
        ilens->effradius = 0.;

        // Read a line of the catalog
        if ( pot->ftype == 1 || pot->ftype == 3 )
        {
            if ( sscanf( line, "%s%lf%lf%lf%lf%lf%lf%lf",
                         ilens->n, &ilens->C.x, &ilens->C.y,
                         &aa, &bb, &ilens->theta, &ilens->mag, &ilens->lum) != 8 )
                goto NEXT;
        }
        else if ( pot->ftype == 2 || pot->ftype == 4 )
        {
            sprintf(ilens->n, "%ld", i);
            if ( sscanf( line, "%d%lf%lf%lf%lf%lf%lf%lf%lf",
                         &ilens->type, &ilens->C.x, &ilens->C.y,
                         &ilens->emass, &ilens->theta, &ilens->rckpc,
                         &ilens->rcutkpc, &ilens->sigma, &ilens->z ) != 9 )
                goto NEXT;
        }
        else if ( pot->ftype == 8 )
        {
            if ( sscanf( line, "%s%lf%lf%lf%lf%lf%lf%lf%lf",
                         ilens->n, &ilens->C.x, &ilens->C.y,
                         &aa, &bb, &ilens->theta, &ilens->mag, &ilens->lum,
                         &ilens->effradius) != 9 )
                goto NEXT;
        }
        else if ( pot->ftype >= 5 )
        {
            if ( sscanf( line, "%s%lf%lf%lf%lf%lf%lf%lf",
                         ilens->n, &ilens->C.x, &ilens->C.y,
                         &aa, &bb, &ilens->theta, &ilens->mag, &ilens->z ) != 8 )
                goto NEXT;
        }

        // Convert input from relative to absolute coordinates if needed
        convertXY( &ilens->C.x, &ilens->C.y, iref, ref_ra, ref_dec );
        convertXY_compat( &ilens->C.x, &ilens->C.y, (iref == 0 && pot->ftype == 1) );

        // convert to output relative coordinates
        if (  M.iref == 1 || M.iref == 3 )
        {
            ilens->C.x -= M.ref_ra;
            ilens->C.x *= -3600 * cos(M.ref_dec * DTR);
            ilens->C.y -= M.ref_dec;
            ilens->C.y *= 3600;
        }
        else if ( M.iref == 2 )
        {
            ilens->C.x -= M.ref_ra;
            ilens->C.y -= M.ref_dec;
        }

        // Eventually modify the ID of the clump
        if ( i < G.no_lens )
        {
            strcpy( line, ilens->n );
            sprintf(ilens->n, "O%s", line);
        }

        // Initialize the lens parameters
        ilens->theta *= DTR;

        if ( pot->ftype != 2 && pot->ftype != 4 )
        {
            ilens->emass = (aa * aa - bb * bb) / (aa * aa + bb * bb);
            if ( ilens->emass < 0 )
            {
                fprintf( stderr, "ERROR: The potfile clump %s has a negative ellipticity.\n", ilens->n );
                exit(-1);
            }
        }

        ilens++;
        i++;

NEXT:
        // get next line
        fgets(line, 128, POT);

    } // endof while loop over the potfile lines.

    fclose(POT);

    if ( i == NLMAX )
    {
        fprintf(stderr, "ERROR: Too many clumps. Maximum allowed %d.\n", NLMAX);
        exit(-1);
    }

    // Return the number of clumps read in the potfile
    return(i - nplens);
}

/* In case of potfile with header, perform an additional conversion
 * for compatibility with previous file format.
 * Perform conversion to absolute coordinates even when the galaxy catalog has no header
 *
 * isAbsWCS is true if galaxy catalog has no #REFERENCE header, but has ftype = 1
 *
 */
static void convertXY_compat(double *x, double *y, int isAbsWCS )
{
    extern  struct g_mode   M;

    // For compatibility with the former syntax
    if ( isAbsWCS )
    {
        *x /= -3600.*cos(M.ref_dec * DTR);
        *x += M.ref_ra;
        *y /= 3600.;
        *y += M.ref_dec;
    }
}


/*
 * Initialise the scaling relations parameters
 */
void setScalingRelations(struct g_pot *pot)
{
    extern struct g_grille  G;
    const extern struct g_cosmo   C;
    extern struct pot       lens[];

    //*********************************************************************
    // Check if the scaling relations are defined
    //*********************************************************************
    if ( pot->ftype == 1 || pot->ftype == 3 )
    {
        if ( pot->sigma1 == -1 )
        {
            fprintf(stderr, "ERROR: potfile: sigma not defined\n");
            exit(-1);
        }

        if ( pot->cutkpc1 == DBL_MAX && pot->cut1 == DBL_MAX )
        {
            fprintf(stderr, "ERROR: potfile: cut length not defined\n");
            exit(-1);
        }

        if ( pot->corekpc == -1 && pot->core == -1 )
        {
            fprintf(stderr, "ERROR: potfile: core length not defined\n");
            exit(-1);
        }
    }
    else if ( pot->ftype == 5 || pot->ftype == 6 )
    {
        if ( pot->a1 == DBL_MAX || pot->b1 == DBL_MAX )
        {
            fprintf(stderr, "ERROR: potfile: a or b scaling parameter not defined\n");
            exit(-1);
        }
    }


    //*********************************************************************
    // Set the Potfile current and limiting values
    //*********************************************************************
    if ( pot->ftype <= 4 )
    {
        // Scale potfile SIGMA
        pot->sigma = pot->sigma1;

        // ... and potfile RCUT
        if ( pot->cut1 == DBL_MAX && pot->cutkpc1 != DBL_MAX )
        {
            pot->cut1 = pot->cutkpc1 / (d0 / C.h * distcosmo1(pot->zlens));
            pot->cut2 = pot->cutkpc2 / (d0 / C.h * distcosmo1(pot->zlens));
        }
        pot->cut = pot->cut1;

        // ... and potfile RCORE
        if ( pot->core == -1.0 && pot->corekpc != -1 )
            pot->core = pot->corekpc / (d0 / C.h * distcosmo1(pot->zlens));

        // ... and potfile RCUT SLOPE
        pot->slope = pot->slope1;

        // ... and potfile VDSLOPE
        pot->vdslope = pot->vdslope1;
    }
    else if ( pot->ftype >= 5 )
    {
        // log(m200) = a * log(mstar) + b (cf Leauthaud2009)
        pot->a = pot->a1;
        pot->b = pot->b1;
        pot->slope = pot->slope1;
        pot->vdslope = pot->vdslope1;
    }

    // set potfile VDSCAT for all potfile scaling relations
    pot->vdscat = pot->vdscat1;

    // ... and potfile RCUTSCAT
    pot->rcutscat = pot->rcutscat1;

    //************************************************************
    // Scale potfile parameters
    //************************************************************
    long int i, ilensmin, ilensmax;
    ilensmin = G.nplens[pot->potid];
    ilensmax = G.nplens[pot->potid+1];
 
    if ( pot->ftype == 1 || pot->ftype == 3 )
        for ( i = ilensmin ; i < ilensmax ; i++ )
            if ( lens[i].mag != 0 )
                lens[i].rc = pot->core * pow(10., 0.4 * (pot->mag0 - lens[i].mag) / 2.);

    scale_pot(pot); // ...and RCUT and SIGMA of the clumps.

}


/*
 * Scale the sigma and rcut of a potfile clump according to the potfile parameters
 */
void scale_pot(struct g_pot *pot)
{
    extern struct g_grille  G;
    extern struct g_source  S;
    extern struct pot       lens[];

    struct pot *ilens;
    double lm200, lc200;
    long int i, ilensmin, ilensmax;
    ilensmin = G.nplens[pot->potid];
    ilensmax = G.nplens[pot->potid+1];

    // loop over the potfile clumps to scale
    if ( pot->ftype <= 4 )
    {
        for ( i = ilensmin; i < ilensmax; i++ )
        {
            ilens = &lens[i];

            if ( ilens->mag != 0 )
            {
                ilens->sigma = pot->sigma *
                               pow(10., 0.4 * (pot->mag0 - ilens->mag) / pot->vdslope);
                /* The factor of 2 so that with slope1 = 4, we have
                 * 2/slope1=1/2, then Brainerd, Blandford, Smail, 1996 */
                ilens->rcut = pot->cut *
                              pow(10., 0.4 * (pot->mag0 - ilens->mag) * 2. / pot->slope);
            }

            if ( pot->ivdscat != 0 )
                ilens->sigma += d_gauss(pot->vdscat, &S.rand);

            // Convert sigma to b0
            set_dynamics(i);

            if ( pot->ircutscat != 0 )
                ilens->rcut += d_gauss(pot->rcutscat, &S.rand);

        }
    }
    else if ( pot->ftype == 5 )
    {
        for ( i = ilensmin; i < ilensmax; i++ )
        {
            ilens = &lens[i];

            // scaling relation from Alexie
            lm200 = ilens->mag * pot->a + pot->b; // here ilens->mag = log10(stellar mass)

            // concentration from Maccio2008 for WMAP2005
            lc200 = -0.098  * (lm200 - 12 ) + 0.830;

            // convert to values for lenstool NFW
            ilens->masse = exp( lm200 * LOG10 );
            ilens->beta = exp( lc200 * LOG10 );
            set_dynamics(i);
        }
    }
    else if ( pot->ftype == 6 )
    {
        for ( i = ilensmin; i < ilensmax; i++ )
        {
            ilens = &lens[i];

            // scaling relation from Alexie
//            lm200 = ilens->mag * pot->a + pot->b; // here ilens->mag = log10(stellar mass)
            lm200 = 0.4 * (pot->mag0 - ilens->mag) * pot->a + pot->b;

            // concentration from Gao2008 (zl=0.5)
            lc200 = -0.125  * lm200 + 2.372;

            // convert to values for lenstool NFW
            ilens->masse = exp( lm200 * LOG10 );
            ilens->beta = exp( lc200 * LOG10 );
            set_dynamics(i);
        }
    }
    else if ( pot->ftype == 61 )
    {
        double lsigma;
        for ( i = ilensmin; i < ilensmax; i++ )
        {
            ilens = &lens[i];

            // scaling relation for SIS
            ilens->sigma =pow(ilens->mag / pot->mag0, pot->a) * pot->b;
            set_dynamics(i);
        }
    }
    else if ( pot->ftype == 62 )
    {
        double lsigma;
        for ( i = ilensmin; i < ilensmax; i++ )
        {
            ilens = &lens[i];

            // scaling relation for NFW and c(M) relation (Mandelbaum et al. 2008)
            ilens->masse = pow(ilens->mag / pot->mag0, pot->slope) * pot->a;
            ilens->beta = pow(ilens->masse / pot->a, -pot->vdslope) * pot->b / pow(1. + ilens->z, 0.45);
            set_dynamics(i);
        }
    }
    else if ( pot->ftype == 63 )
    {
        double lsigma;
        for ( i = ilensmin; i < ilensmax; i++ )
        {
            ilens = &lens[i];

            // scaling relation for NFW and c(M) relation (Mandelbaum et al. 2008)
            // same as 62 but for point mass, we don't need beta
            ilens->masse = 0.01 * pow(ilens->mag / pot->mag0, pot->slope) * pot->a;
            ilens->masse /= 1e12; // definition of point mass in Lenstool
            set_dynamics(i);
        }
    }
    else if ( pot->ftype == 7 )
    {
        const extern struct g_cosmo  C;
        double lmvir, mvir, cvir, rvir;
        double rho_c, x;
        double rs, rhos, delta_c, sigma_s2;

        for ( i = ilensmin; i < ilensmax; i++ )
        {
            ilens = &lens[i];

            // scaling relation from Bundy 2007 (stellar mass to mvir)
            //lmvir = ilens->mag * pot->a + pot->b - 0.2 * (1. + ilens->z) ; 
            lmvir = (ilens->mag - 11. + 2. * log(C.h)) + 12.83 - 0.2 * (1. + ilens->z) ; 
            mvir = exp(lmvir * LOG10);

            // concentration from Bullock 2001
            cvir = 9. * pow(mvir / 8.12e12 * C.h, -0.14) / (1. + ilens->z);

            // Delta_c overdensity for a flat universe and wX = -1, Bryan & Norman 1998
            x = C.omegaM * pow(1. + ilens->z, 3) * chiz(ilens->z) * chiz(ilens->z);
            x = x - 1.;
            delta_c = 18. * PI * PI + 60. * x -  32 * x * x;  // Overdensity 

            // rs in arcsec
            rho_c = rho_cri(ilens->z);   // in Msun/Mpc^3
            rvir = pow(3. * mvir / 4. / PI / rho_c / delta_c, 0.333333);  // in Mpc
            rs = rvir / cvir;  // in Mpc
            ilens->rc = rs * 1e3 / (d0 / C.h * distcosmo1(ilens->z));
            
            // sigma_s from Golse 2002
            delta_c = delta_c / 3. * pow(cvir, 3) / (log(1. + cvir) - cvir / (1. + cvir));
            rhos = delta_c * rho_c;  // in Msun/Mpc^3
            sigma_s2 = 8. / 3. / INVG * rs * rs * rhos;
            ilens->sigma = sqrt(sigma_s2);
            ilens->b0 = 6.*pia_c2 * ilens->sigma * ilens->sigma;
        }
    }
    else if ( pot->ftype == 8 )
    {
        for ( i = ilensmin; i < ilensmax; i++ )
        {
            ilens = &lens[i];//lens[i] contsains parameters for potential i

            if ( ilens->mag != 0 )
            {
                if ( ilens->effradius != 0 ) {
                    ilens->sigma = pot->sigma *
                                   pow(10., log(ilens->effradius * 0.05) / 1.49 -
                                          0.74/1.49 * (ilens->mag + 2.5 * log(2. * PI * (0.05*ilens->effradius)*(0.05*ilens->effradius)))/2.5  + 8.779/1.49); ///Bernardi:2003  R0 = sigma^1.49 I^-0.75 --> sigma = R0^1./1.49 . I^0.75/1.49

                    ilens->rcut = pot->cut *
                                  pow(10., 0.4 * (pot->mag0 - ilens->mag) * 2. / pot->slope);
                    }
                else {
                    ilens->sigma = pot->sigma *
                               pow(10., 0.4 * (pot->mag0 - ilens->mag) / pot->vdslope);
                    /* The factor of 2 so that with slope1 = 4, we have
                     * 2/slope1=1/2, then Brainerd, Blandford, Smail, 1996 */
                    ilens->rcut = pot->cut *
                                  pow(10., 0.4 * (pot->mag0 - ilens->mag) * 2. / pot->slope);
                }
            }

            if ( pot->ivdscat != 0 )
                ilens->sigma += d_gauss(pot->vdscat, &S.rand);

            // Convert sigma to b0
            set_dynamics(i);//converts velocity dispersion into different units

            if ( pot->ircutscat != 0 )
                ilens->rcut += d_gauss(pot->rcutscat, &S.rand);

        }
    }


}
