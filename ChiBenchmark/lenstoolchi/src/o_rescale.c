#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "constant.h"
#include "dimension.h"
#include "structure.h"
#include "fonction.h"

/********************************************************/
/*      fonction: recale            */
/*      auteur: EJ              */
/********************************************************
 * Rescale the cube values to their corresponding parameter values
 * in Lenstool. Each atom corresponds to a specific clump.
 * Fill the lens global variable.
 *
 * Return 1 if there are the number of clumps in the cubes is equal to the
 * number of lenses, 0 otherwise which means the object has to be rejected.
 *
 * Parameters :
 * - cubes[clump][param] : the cubes containing the parameters for all clumps
 * - nclumps : the number of clumps in the cubes
 * - npar : the number of parameters
 *
 * Global variables used :
 */

#define POTENTIAL 0
#define COSMO 1
#define ZMLIMIT 2
#define ZALIMIT 3
#define SIGPOS 4
#define SOURCES 5
#define VFIELD 6
#define POT 7   // WARNING keep POT at the bottom of the list


static double   prior(int k, double val, double x1, double x2);

/* Rescale the lens parameter ipx for the lens ilens and test the
 * physical validy of the rescaled values.
 *
 * Parameters :
 * -ilens : lens index in lens global variable
 * -ipx : parameter index in the Param enumerated type
 * -val : value from the cube in the [0..1] range
 *
 * Return 1 if it's ok, 0 if it's not physicaly valid
 *
 * Global variables used :
 * -lmin, lmax
 */
int rescaleParam(int id, int type, int ipx, double *pval)
{

    double val = *pval;

    int valid = 1;

    if ( val == -1 )
    { 
        *pval = 0.; 
        return valid; 
    }

    if ( type == POTENTIAL || type == COSMO )
    {
        extern int block[][NPAMAX];
        extern int cblock[NPAMAX];
        extern struct pot         lmin[], lmax[];
        extern struct g_cosmo         clmin, clmax;

        switch (ipx)
        {
            case(CX):
                val = prior( block[id][CX], val,
                             lmin[id].C.x, lmax[id].C.x );
                break;
            case(CY):
                val = prior( block[id][CY], val,
                             lmin[id].C.y, lmax[id].C.y );
                break;
            case(EPOT):
                val = prior( block[id][EPOT], val,
                             lmin[id].epot, lmax[id].epot );
                if ( val < 0 ) valid = 0;
            break;
            case(EMASS):
                val = prior( block[id][EMASS], val,
                             lmin[id].emass, lmax[id].emass );
                if ( val < 0 ) valid = 0;
                break;
            case(THETA):
                val = prior( block[id][THETA], val,
                             lmin[id].theta, lmax[id].theta );
                break;
			case(PHI):
				val = prior( block[id][PHI], val, 
					lmin[id].phi, lmax[id].phi );
				break;
            case(RC):
                val = prior( block[id][RC], val,
                             lmin[id].rc, lmax[id].rc );
                if ( val < 0 ) valid = 0;
                break;
            case(B0):
                // rescaling to sigma instead of B0 to
                // keep the shape of the input prior.
                val = prior( block[id][B0], val,
                             lmin[id].sigma, lmax[id].sigma );
                if ( val < 0 ) valid = 0;
                break;
            case(ALPHA):
                val = prior( block[id][ALPHA], val,
                             lmin[id].alpha, lmax[id].alpha );
                break;
            case(BETA):
                val = prior( block[id][BETA], val,
                             lmin[id].beta, lmax[id].beta );
                break;
            case(RCUT):
                val = prior( block[id][RCUT], val,
                             lmin[id].rcut, lmax[id].rcut );
                if ( val < 0 ) valid = 0;
                break;
            case(PMASS):
                val = prior( block[id][PMASS], val,
                             lmin[id].pmass, lmax[id].pmass );
                if ( val < 0 ) valid = 0;
                break;
            case(MASSE):
                val = prior( block[id][MASSE], val,
                             lmin[id].masse, lmax[id].masse );
                if ( val < 0 ) valid = 0;
                break;
            case(ZLENS):
                val = prior( block[id][ZLENS], val,
                             lmin[id].z, lmax[id].z );
                if ( val < 0 ) valid = 0;
                break;
            case(RCSLOPE):
                val = prior( block[id][RCSLOPE], val,
                             lmin[id].rcslope, lmax[id].rcslope );
                if ( val < 0 ) valid = 0;
                break;
            case(OMEGAM):
                val = prior( cblock[OMEGAM], val,
                             clmin.omegaM, clmax.omegaM );
                // Test physical validity. beause: H(z)^2 can become negative
                // in the corner OmegaX > 1 and OmegaM ~ 0, ie Universes with "rebond"
                if ( val < 0 ) valid = 0;
                break;
            case(OMEGAX):
                val = prior( cblock[OMEGAX], val,
                             clmin.omegaX, clmax.omegaX );
                if ( val > 1 ) valid = 0;
                break;
            case(WX):
                val = prior( cblock[WX], val,
                             clmin.wX, clmax.wX );
                break;
            case(WA):
                val = prior( cblock[WA], val,
                             clmin.wa, clmax.wa );
                break;
        };
        //if ( valid && ipx != B0 )
        //    o_set_lens( id, ipx, val );
        //if ( ipx == THETA ) val *= RTD; // for output purpose through the Cube array

    }

// zm_limit section
    if ( type == ZMLIMIT )
    {
        extern struct z_lim     zlim[];
        val = prior( zlim[ipx].bk, val, zlim[ipx].min, zlim[ipx].max );
        if ( val < 0 ) valid = 0;
    }

// z_a_limit section
    if ( type == ZALIMIT )
    {
        extern struct z_lim     zalim;
        val = prior( zalim.bk, val, zalim.min, zalim.max );
        if ( val < 0 ) valid = 0;
    }

// vfield section
    if ( type == VFIELD )
    {
        extern struct vfield  vfmin,vfmax;
        extern int vfblock[NPAMAX];
        switch(ipx)
        {
            case(VFCX):
               val = prior ( vfblock[VFCX], val, vfmin.C.x, vfmax.C.x );
               break;
            case(VFCY):
               val = prior ( vfblock[VFCY], val, vfmin.C.y, vfmax.C.y );
               break;
            case(VFVT):
               val = prior ( vfblock[VFVT], val, vfmin.vt, vfmax.vt );
               if ( val < 0 ) valid = 0;
               break;
            case(VFRT):
               val = prior ( vfblock[VFRT], val, vfmin.rt , vfmax.rt );
               if ( val < 0 ) valid = 0;
               break;
            case(VFI):
               val = prior ( vfblock[VFI], val, vfmin.i, vfmax.i );
               break;
            case(VFTHETA):
               val = prior ( vfblock[VFTHETA], val, vfmin.theta, vfmax.theta );
               break;
            case(VFLCENT):
               val = prior ( vfblock[VFLCENT], val, vfmin.lcent, vfmax.lcent );
               if ( val < 0 ) valid = 0;
               break;
            case(VFSIGMA):
               val = prior ( vfblock[VFSIGMA], val, vfmin.sigma, vfmax.sigma );
               if ( val < 0 ) valid = 0;
               break;
        };
    }
 
// rescale the sigposAs noise
    if ( type == SIGPOS )
    {
        extern struct g_image   I;
        extern struct sigposStr sigposAs;

        switch (ipx)
        {
            case(0): // sigposArsec 
                val = prior( sigposAs.bk, val, sigposAs.min, sigposAs.max );
                if ( val <= 0 ) valid = 0;
                break;
            case(1): // sigell
                val = prior( 3, val, I.sigell, I.dsigell );
                if ( val <= 0 ) valid = 0;
                break;
        }
    }

// source positions for iclean == 2
    if ( type == SOURCES )
    {
        extern int sblock[NFMAX][NPAMAX];
        extern struct galaxie smin[NFMAX], smax[NFMAX];

        switch( ipx )
        {
            case(SCX):
                val = prior(sblock[id][SCX], val, smin[id].C.x, smax[id].C.x);
                break;
            case(SCY):
                val = prior(sblock[id][SCY], val, smin[id].C.y, smax[id].C.y);
                break;
            case(SA):
                val = prior(sblock[id][SA], val, smin[id].E.a, smax[id].E.a);
                break;
            case(SB):
                val = prior(sblock[id][SB], val, smin[id].E.b, smax[id].E.b);
                break;
            case(SEPS):
                val = prior(sblock[id][SEPS], val, smin[id].eps, smax[id].eps);
                break;
            case(STHETA):
                val = prior(sblock[id][STHETA], val, smin[id].E.theta, smax[id].E.theta);
                break;
            case(SFLUX):
                val = prior(sblock[id][SFLUX], val, smin[id].mag, smax[id].mag);
                break;
        }
    }

// rescale the potfile parameters rcut and sigma
    if ( type >= POT ) 
    {
        extern struct g_pot    P[NPOTFILE];
        struct g_pot *pot = &P[type - POT];
	const double max_dmag = 10.;  // difference between brightest gal and mag0 (see set_potfile.c)

        switch (ipx)
        {
            case(0):
                val = prior( pot->ircut, val, pot->cut1, pot->cut2 );
                if ( val <= 0 ) valid = 0;
                break;
            case(1):
                val = prior( pot->isigma, val, pot->sigma1, pot->sigma2 );
                if ( val <= 0 ) valid = 0;
                break;
            case(2):
                val = prior( pot->islope, val, pot->slope1, pot->slope2 );
        	if ( 0.8 * max_dmag / val > DBL_MAX_10_EXP || val <= 0. ) valid = 0;
                break;
            case(3):
                val = prior( pot->ivdslope, val, pot->vdslope1, pot->vdslope2 );
        	if ( 0.4 * max_dmag / val > DBL_MAX_10_EXP || val <= 0. ) valid = 0;
                break;
            case(4):
                val = prior( pot->ivdscat, val, pot->vdscat1, pot->vdscat2 );
                if ( val <= 0 ) valid = 0;
                break;
            case(5):
                val = prior( pot->ircutscat, val, pot->rcutscat1, pot->rcutscat2 );
                if ( val <= 0 ) valid = 0;
                break;
            case(6):
                val = prior( pot->ia, val, pot->a1, pot->a2 );
                break;
            case(7):
                val = prior( pot->ib, val, pot->b1, pot->b2 );
                break;
        }
    }

    *pval = val;
    return valid;
}

static double uniformPrior(double val, double min, double max);
static double gaussianPrior(double val, double mu, double sigma);

/* Prior distribution functions : val is a deviate from the unit Cube.
 *
 * Return the corresponding deviate drawn from the desired distribution
 * (labelled k).
 */
static double prior(int k, double val, double x1, double x2)
{
    switch (k)
    {
        case(1) :
            return uniformPrior(val, x1, x2);
            break;
        case(3) :
            return gaussianPrior(val, x1, x2);
            break;
        default :
            fprintf( stderr, "ERROR : prior index %d undefined.\n", k );
            exit(-1);
    }

    return 0; //to avoid the compiler warning.
}

/* Uniform [0:1] -> Uniform[x1:x2]
 */
static double uniformPrior(double val, double min, double max)
{
    return  min + val * (max - min);
}

double dierfc(double y);
/* Uniform [0:1] -> Gaussian[mean=mu, variance=sigma**2]
 */
static double gaussianPrior(double val, double mu, double sigma)
{
    double  sqrtTwo = 1.414213562;

    return mu + sigma*sqrtTwo*dierfc(2. * (1. - val) );
}

/* Inverse of errror function in double precision
 */
static const int n_ck = 60;
static const double ck[60] = { 1, 1, 1.16666666666667, 1.41111111111111,
             1.73373015873016, 2.14858024691358, 2.67716623911068,
             3.34814636128128, 4.19849396342683, 5.27542268646125,
             6.63899509461546, 8.36550450191292, 10.5518020210761,
             13.3208081702262, 16.8285215357839, 21.2729253053784,
             26.9053027832427, 34.0446123102593, 43.0957486538823,
             54.5727422435001, 69.1282326361197, 87.5909148253449, 
             111.013117434205, 140.731257124519, 178.442657611487,
             226.303167593896, 287.051214503983, 364.165459938403,
             462.065166575814, 586.364858016593, 744.197995615538, 
             944.628392281556, 1199.17316418125, 1522.46748209153,
             1933.10959970608, 2454.73508332157, 3117.38245243559,
             3959.22933516056, 5028.7997269694, 6387.77026383145,  
             8114.53816822392, 10308.7577173092, 13097.108284167,
             16640.6284810277, 21144.0418418296, 26867.6151038299,
             34142.2372053783, 43388.5941585717, 55141.5528563632,
             70081.1694675158, 89072.1229570381, 113213.863830027,
             143904.390909929, 182921.361053908, 232525.244262235,
             295590.518280687, 375772.5271056, 477719.701659528,
             607343.479014971, 772161.612462455};

double dierfc(double y)
{
#define sqrtpio2 0.886226925452758
    int k;
    double val = 0.;

    y = 1. - y;
    for( k = 0; k < n_ck; k++ )
        val += ck[k] / (2 * k + 1) * pow(sqrtpio2 * y, 2 * k + 1);

    return val;
}

/* Inverse of error function in double precision
 */
//static double dierfc(double y)
//{
//#define  qa  9.16461398268964e-01
//#define  qb  2.31729200323405e-01
//#define  qc  4.88826640273108e-01
//#define  qd  1.24610454613712e-01
//#define  q0  4.99999303439796e-01
//#define  q1  1.16065025341614e-01
//#define  q2  1.50689047360223e-01
//#define  q3  2.69999308670029e-01
//#define  q4  -7.28846765585675e-02
//
//#define  pa  3.97886080735226000e+00
//#define  pb  1.20782237635245222e-01
//#define  p0  2.44044510593190935e-01
//#define  p1  4.34397492331430115e-01
//#define  p2  6.86265948274097816e-01
//#define  p3  9.56464974744799006e-01
//#define  p4  1.16374581931560831e+00
//#define  p5  1.21448730779995237e+00
//#define  p6  1.05375024970847138e+00
//#define  p7  7.13657635868730364e-01
//#define  p8  3.16847638520135944e-01
//#define  p9  1.47297938331485121e-02
//#define  p10  -1.05872177941595488e-01
//#define  p11  -7.43424357241784861e-02
//#define  p12  2.20995927012179067e-03
//#define  p13  3.46494207789099922e-02
//#define  p14  1.42961988697898018e-02
//#define  p15  -1.18598117047771100e-02
//#define  p16  -1.12749169332504870e-02
//#define  p17  3.39721910367775861e-03
//#define  p18  6.85649426074558612e-03
//#define  p19  -7.71708358954120939e-04
//#define  p20  -3.51287146129100025e-03
//#define  p21  1.05739299623423047e-04
//#define  p22  1.12648096188977922e-03
//
//    double infinity = 5.0;
//    double z, w, u, s, t, x;
//
//    if ( y < 1.e-40 )
//        return infinity;
//
//    z = y;
//    if ( y > 1. ) z = 2 - y;
//    w = qa - log10(z);
//    u = sqrt(w);
//    s = (qc + log10(u)) / w;
//    t = 1. / (u + qb);
//    x = u * (1. - s * (0.5 + s * qd)) -
//        ((((q4 * t + q3) * t + q2) * t + q1) * t + q0) * t;
//    t = pa / (pa + x);
//    u = t - 0.5;
//    s = (((((((((p22 * u + p21) * u + p20) * u +
//               p19) * u + p18) * u + p17) * u + p16) * u +
//           p15) * u + p14) * u + p13) * u + p12;
//    s = ((((((((((((s * u + p11) * u + p10) * u +
//                  p9) * u + p8) * u + p7) * u + p6) * u + p5) * u +
//             p4) * u + p3) * u + p2) * u + p1) * u + p0) * t -
//        z * exp(x * x - pb);
//    x = x + s * (1. + x * s);
//    x = x + s * (1. + x * s);
//    if ( y > 1.) x = -x;
//
//    return x;
//}

/********************************************************
* Rescale the cube values to their corresponding parameter values
* in Lenstool. A single atom contains all the clumps
* Fill the lens global variable.
*
* Return 1 if the rescaled parameters have significant meaning,
* 0 otherwise.
*
* Parameters :
* - cube[param] : the cubes containing the parameters for all clumps
* - npar : the number of parameters
*
* Global variables used :
* -
*/
int rescaleCube_1Atom(double *cube, int npar)
{
    extern struct g_mode M;
    extern struct g_grille  G;
    extern struct g_source  S;
    extern struct g_cosmo   C;
    extern struct g_image   I;
    extern struct g_pot P[NPOTFILE];
    extern struct pot       lens[];
    extern struct galaxie   multi[NFMAX][NIMAX];
    extern struct z_lim zlim[];
    extern struct z_lim zalim;
    extern struct pot   lens[], lmin[], lmax[];
    extern struct sigposStr sigposAs;

    extern int block[][NPAMAX];
    extern int cblock[NPAMAX];
    extern int sblock[NFMAX][NPAMAX];
    extern int vfblock[NPAMAX];

    extern double *np_b0;  // non parametric accelerator (o_global.c)

    int ipar;   /*index of the current parameter in the atom*/
    int ipx;    /*index of the current parameter in the ipot structure*/
    long int ilens;  /*index of the clump in the lens[] global variable*/
    int valid;      // Accept or not this sample
    long int i;
    int j, k;
    double tsig, tlmax, tlmin, tmp;

    valid = 1;      // OK

// Rescale de potentials
    ipar = 0;
    for ( ilens = 0; ilens < G.no_lens; ilens++ )
    {
        // switch from rectangular to circular prior on position
        if( block[ilens][CX] == 4 && block[ilens][CX] == 4 )
        {
            double r = 0.5 * cube[ipar];
            double theta = cube[ipar+1] * 2. * M_PI;
            cube[ipar] = 0.5 + r * cos(theta);
            cube[ipar+1] = 0.5 + r * sin(theta);
        }
		
		// normal scaling for all the parameters
        for ( ipx = CX; ipx <= PMASS; ipx++ )
            if ( block[ilens][ipx] != 0 )
                valid *= rescaleParam(ilens, POTENTIAL, ipx, &cube[ipar++]);
    }

// Rescale the clumps of the multiscale grid
    for ( ilens = G.nmsgrid ; ilens < G.nlens ; ilens++ )
       rescaleParam(ilens, POTENTIAL, B0, &cube[ipar++]); // we don't care if sigma < 0

// Rescale the source parameters
    if ( M.iclean == 2 )
        for ( k = 0; k < S.ns; k++ )
            for ( ipx = SCX; ipx <= SFLUX; ipx++ )
                if ( sblock[k][ipx] != 0 )
                    valid *= rescaleParam(k, SOURCES, ipx, &cube[ipar++]);

// rescale the cosmological parameters
    for ( ipx = OMEGAM; ipx <= WA; ipx++ )
        if ( cblock[ipx] != 0 )
            valid *= rescaleParam(0, COSMO, ipx, &cube[ipar++]);

// rescale the z_m_limit redshifts
    for ( k = 0; k < I.nzlim; k++ )
        if ( zlim[k].opt != 0 )
            valid *= rescaleParam(0, ZMLIMIT, k, &cube[ipar++]);

// rescale the redshift of the arclets
        if ( zalim.bk != 0 )
            valid *= rescaleParam(0, ZALIMIT, 0, &cube[ipar++]);

// rescale the velocity field parameters
    for ( ipx = VFCX; ipx <= VFSIGMA; ipx++ )
        if ( vfblock[ipx] != 0 )
            valid *= rescaleParam(0, VFIELD, ipx, &cube[ipar++]);

// rescale the pot parameters
    for( k = 0; k < G.npot; k++ )
        if ( P[k].ftype != 0 )
        {
            if ( P[k].ircut != 0 )
                valid *= rescaleParam(0, POT + k, 0,  &cube[ipar++]);
            if ( P[k].isigma != 0 )
                valid *= rescaleParam(0, POT + k, 1,  &cube[ipar++]);
            if ( P[k].islope != 0 )
                valid *= rescaleParam(0, POT + k, 2,  &cube[ipar++]);
            if ( P[k].ivdslope != 0 )
                valid *= rescaleParam(0, POT + k, 3,  &cube[ipar++]);
            if ( P[k].ivdscat != 0 )
                valid *= rescaleParam(0, POT + k, 4,  &cube[ipar++]);
            if ( P[k].ircutscat != 0 )
                valid *= rescaleParam(0, POT + k, 5,  &cube[ipar++]);
            if ( P[k].ia != 0 )
                valid *= rescaleParam(0, POT + k, 6,  &cube[ipar++]);
            if ( P[k].ib != 0 )
                valid *= rescaleParam(0, POT + k, 7,  &cube[ipar++]);
        }
 
// rescale the noise
    if ( sigposAs.bk != 0. )
        valid *= rescaleParam(0, SIGPOS, 0, &cube[ipar++]);
    
    if ( I.dsigell != -1. )
        valid *= rescaleParam(0, SIGPOS, 1,  &cube[ipar++]);

    // no need to go further??
    if ( valid == 0 ) return valid;

    setBayesModel( -5, 1, &cube);
    
    // final checks
    for ( ilens = 0; ilens < G.no_lens; ilens++ )
        if ( lens[ilens].b0 != 0. && lens[ilens].rc > lens[ilens].rcut )   valid = 0;

    if ( M.iclean == 2 )
    {
        extern struct galaxie source[NFMAX];
        for ( k = 0; k < S.ns; k++ )
            if( source[k].E.a < source[k].E.b )   valid = 0;
    }

    return valid;
}


int rescaleCube_NAtoms(double **cubes, int natoms, int npar)
{
    extern struct g_grille  G;
    extern struct g_mode    M;
    extern int block[][NPAMAX];
    extern struct pot   lens[];

    int ipar;   /*index of the current parameter in the atom*/
    int ipx;    /*index of the current parameter in the ipot structure*/
    int valid;      // Accept or not this sample
    int ilens; 

    valid = 1;      // OK

    G.no_lens = natoms;
    for ( ilens = 0; ilens < G.no_lens; ilens++ )
    {
        ipar = 0;
        for ( ipx = CX; ipx <= PMASS; ipx++ )
            if ( block[ilens][ipx] != 0 )
                valid *= rescaleParam(ilens, POTENTIAL, ipx, &cubes[ilens][ipar++]);
    }
 
    // no need to go further??
    if ( valid == 0 ) return valid;


    return valid;
}

