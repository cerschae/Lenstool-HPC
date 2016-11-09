#include<stdio.h>
#include<signal.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*      nom:        g_mass              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

typedef void (*sighandler_t)(int);

static void computeKmap(int ny, int nx, double dx, double dy, double **map, double zl, double zs);
static void signalReset(int);
static int  interrupt;  // Global variable

void    g_mass(int imass, int np, double zl, double zs, char *file)
{
    const extern struct g_mode    M;
    const extern struct g_frame   F;
    const extern struct g_cosmo   C;
    const extern struct pot       lens[];

    double  dl0s, dos;
    double  **mass;     // mass map
    int     nx, ny;     // image size
    double  dx, dy;
    double  dlsds, dl, dcrit, dcritA, conv;
    register int    i, j;

    if (imass == 1)
    {
        if( zl == 0 )
            NPRINTF(stderr, "COMP: projected kappa map for z_s=%.3lf =>%s\n", zs, file);
        else
            NPRINTF(stderr, "COMP: kappa map for z_l=%.3lf and z_s=%.3lf =>%s\n", zl, zs, file);
    }
    else if (imass == 2)
    {
        if( zl == 0 )
            NPRINTF(stderr, "COMP: projected absolute mass map in g/cm2 (normalized at z_l=%.3lf) =>%s\n", lens[0].z, file);
        else
            NPRINTF(stderr, "COMP: absolute mass map in g/cm2 for z_l=%.3lf =>%s\n", zl, file);
    }
    else if ( imass == 3)
    {
        if( zl == 0 )
            NPRINTF(stderr, "COMP: projected absolute mass map in 10^12 Msol/pixel (normalized at z_l=%.3lf) =>%s\n", lens[0].z, file);
        else
            NPRINTF(stderr, "COMP: absolute mass map in 10^12 Msol/pixel for z_l=%.3lf =>%s\n", zl, file);
    }
    else if (imass == 4)
    {
        if( zl == 0 )
            NPRINTF(stderr, "COMP: projected absolute mass map in 10^12 Msol/kpc2 (normalized at z_l=%.3lf) =>%s\n", lens[0].z, file);
        else
            NPRINTF(stderr, "COMP: absolute mass map in 10^12 Msol/kpc2 for z_l=%.3lf =>%s\n", zl, file);
    }
    else
    {
        NPRINTF(stderr,
                "ERROR: mass keyword modifier %d in runmode section unknown\n", imass);
        return;
    }

    // Parse errors
    if ( np == 0 )
    {
        fprintf(stderr, "ERROR: mass map size set to 0\n");
        exit(1);
    }

    if( zl >= zs )
    {
        fprintf(stderr, "ERROR: lens plane redshift z_l=%.3lf larger than source plane redshift z_s=%.3lf\n", zl, zs);
        exit(1);
    }

    // Initialise the mass map(ny,nx)
    nx = np;
    ny = (int) (nx / (F.xmax - F.xmin) * (F.ymax - F.ymin));
    dx = (F.xmax - F.xmin) / (nx - 1);
    dy = (F.ymax - F.ymin) / (ny - 1);
    NPRINTF(stderr, "dx:%lf dy:%lf\n", dx, dy);

    mass = (double **) alloc_square_double(ny, nx);

    // Compute the convergence map
    computeKmap(ny, nx, dx, dy, mass, zl, zs);

    // Convert from convergence to the final map type
    if( zl == 0 )
        zl = lens[0].z;

    dlsds = dratio(zl, zs);
    dl = distcosmo1(zl);
    dcrit = cH4piG * C.h / dl / dlsds;  // in g/cm^2
    dcritA = cH0_4piG * C.h / dl / dlsds;  // in 10^12 M_sol/kpc^2
    conv = MCRIT12 / C.h * dx * dy * dl / dlsds;  // in  10^12 M_sol/pixel
    for (j = 0; j < ny; j++)
        for (i = 0; i < nx; i++)
        {
            if ( imass == 2 )
                mass[j][i] *= dcrit;
            else if ( imass == 3 )
                mass[j][i] *= conv;
            else if ( imass == 4 )
                mass[j][i] *= dcritA;
        }

    // Write the mass maps
    if ( M.iref > 0 )
        wrf_fits_abs(file, mass, nx, ny, F.xmin, F.xmax, F.ymin, F.ymax, M.ref_ra, M.ref_dec);
    else
        wrf_fits(file, mass, nx, ny, F.xmin, F.xmax, F.ymin, F.ymax);

    // Free the arrays
    free_square_double(mass, ny);
}

/* Return a convergence map of size npxnp in pixels
 */
static void computeKmap(int ny, int nx, double dx, double dy, double **map, double zl, double zs)
{
    const extern struct g_grille    G;
    const extern struct g_frame     F;
    const extern struct g_mode      M;
    const extern struct pot    lens[]; 
    struct point pi;
    struct matrix grad2;
    int i, j;
    long int ilens;
    double oldz, dls;

    // Compute dls for lens[0].z (default for projected kappa)
    dls = distcosmo2(lens[0].z, zs); 
    oldz = lens[0].z;
    
    // For kappa at a given redshift
    if( zl != 0 )
    {
        dls = distcosmo2(zl, zs);
        oldz = zl;
    }

    for ( ilens = 0; ilens < G.nlens; ilens ++)
    {
        if( lens[ilens].z >= zs ) 
        {
            NPRINTF(stderr, "WARN: Skip potential %ld at z_l=%.3lf larger than z_s=%.3lf\n", ilens, lens[ilens].z, zs);
            continue;
        }

        if( zl != 0 && lens[ilens].z != zl )
            continue;

        printf( "INFO: Compute map for lens %ld/%ld\r", ilens, G.nlens); fflush(stdout);

        // This is for projected kappa
        if( lens[ilens].z != oldz )
        {
            dls = distcosmo2(lens[ilens].z, zs);
            oldz = lens[ilens].z;
        }

        for (j = 0; j < ny; j++)
        {
            pi.y = j * dy + F.ymin;
            for (i = 0; i < nx; i++)
            {
                pi.x = i * dx + F.xmin;
                grad2 = e_grad2_pot(&pi, ilens);
                map[j][i] += 0.5 * (grad2.a + grad2.c) * dls; 
            }
        }
    }

    // Renormalize everything by DOS
    double dos = distcosmo1(zs);
    for( j = 0; j < ny; j++ )
        for( i = 0; i < nx; i++ )
            map[j][i] /= dos;
}

static void signalReset(int param)
{
    interrupt = 1;
    signal(SIGINT, SIG_DFL);
}
