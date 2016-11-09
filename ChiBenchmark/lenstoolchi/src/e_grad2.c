#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
//#include "omp.h"

#define cube(A) A*A*A

static double par1(double x, double y, const struct pot *ilens);
static double par2(double x, double y, const struct pot *ilens);
static struct matrix e_grad2_np(struct galaxie *imag, double * np_b0_thread);

/****************************************************************/
/*      nom:        e_grad2             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Return the projected lenses potential computed in pi. It is
 * the sum over all lenses of their individual potential second
 * derivatives computed in pi.
 *
 * This function doesnt call functions that use and modify global variables.
 *
 * Parameters :
 * - pi : position of computation in the image plane
 *
 * Global variables used :
 * - G, lens
 * - in par1() : lens
 * - in par2() : lens
 * - in ngwg_kappa() : lens_table
 * - in nfwg_gamma() : lens_table
 * - in nfwg_kappa_eps() : lens_table
 * - in nfwg_gamma_eps() : lens_table
 */
struct  matrix  e_grad2(const struct point *pi, double dl0s, double zs)
{
    const extern  struct  g_grille    G;
    const extern struct pot       lens[];

    struct matrix  MA, grad2;
    double dls, oldz;
    long int i;

    MA.a = MA.b = MA.c = MA.d = 0.;

    /*for each lens*/
    oldz = lens[0].z; dls = dl0s;
    for ( i = 0; i < G.nlens; i++ )
    {
        if( lens[i].z >= zs ) 
            continue;
        
        if( lens[i].z != oldz )
        {
            dls = distcosmo2( lens[i].z, zs);
            oldz = lens[i].z;
        }

        grad2 = e_grad2_pot(pi, i);
        MA.a += grad2.a * dls;
        MA.b += grad2.b * dls;
        MA.c += grad2.c * dls;
    }

    MA.d = MA.b;

    return (MA);
}

/* Special e_grad2() function for an image or arclet  */
struct matrix e_grad2_gal(struct galaxie *image, double *np_b0)
{
    const extern struct g_grille  G;
    const extern struct g_pot     P[NPOTFILE];
    const extern struct pot       lens[];
    const extern struct sigposStr sigposAs;
   
    struct matrix MA, grad2;
    struct matrix  *igrad2;  // image's grad2 matrix of not optimised clumps
    double dx, dy, u;
    double dls  = 0;
    double oldz = -1;  //should be initilised in case of G.nmsgrid <= 0!
    int skip_clump; // if 1, skip the gradient computation for a clump
    long int i;    
    int j;

    MA.a = MA.b = MA.c = MA.d = 0;  // final grad2 matrix
    igrad2 = &image->grad2;

    // get the potential of the optimised clumps
    oldz = lens[0].z; dls = image->dl0s;
    for ( i = 0; i < G.no_lens; i++ )
    {
        if( lens[i].z > image->z )
            continue;

        if( lens[i].z != oldz )
        {
	        dls  = distcosmo2(lens[i].z, image->z);
            oldz = lens[i].z;
        }
       	 
        grad2 = e_grad2_pot(&image->C, i);
        MA.a += grad2.a * dls;
        MA.b += grad2.b * dls;
        MA.c += grad2.c * dls;
      
    }

    // Scaling relation clumps (G.no_lens -> G.nmsgrid ) 
    // if not defined, compute the potential of the not optimised clumps
    if ( igrad2->a == igrad2->c )
    {
        igrad2->a = igrad2->b = 0.;
        igrad2->c = 1e-10; // to block images that skip all pots

        // Potentials no optimized but defined individually
        for ( i = G.no_lens; i < G.nplens[0]; i++ )
        {
            if( lens[i].z > image->z )
                continue;

            if( lens[i].z != oldz )
            {
                dls = distcosmo2(lens[i].z, image->z);
                oldz = lens[i].z;
            }
            
            grad2 = e_grad2_pot(&image->C, i);
            igrad2->a += grad2.a * dls;
            igrad2->b += grad2.b * dls;
            igrad2->c += grad2.c * dls;
        }

        // Potentials in potfiles
        for ( j = 0; j < G.npot; j++ )
            for ( i = G.nplens[j]; i < G.nplens[j+1]; i++ )
            {
                skip_clump = 0;  // do not skip the gradient computation (Default)
    
                if ( P[j].select == 1 )
                {
                    // test if the deflexion produced by this clump is detectable
                    // assuming a SIS potential in first approx
                    dx = image->C.x - lens[i].C.x;
                    dy = image->C.y - lens[i].C.y;
                    u = sqrt(dx * dx + dy * dy);
                    dx = lens[i].b0 * dx / u * image->dr;
                    dy = lens[i].b0 * dy / u * image->dr;
                    if ( dx*dx + dy*dy > sigposAs.max*sigposAs.max )
                        skip_clump = 1;
                }
                else if( P[j].select == 2 )
                {
                    // test on the distance to the image in arcsec
                    // work for NFW and PIEMD grid clumps
                    dx = image->C.x - lens[i].C.x;
                    dy = image->C.y - lens[i].C.y;
                    u = sqrt(dx * dx + dy * dy);
                    if ( u > 10 * lens[i].rc ) skip_clump = 1;
                }
		else if( P[j].select > 2 )
                {
                    // test on the distance to the image in arcsec
                    dx = image->C.x - lens[i].C.x;
                    dy = image->C.y - lens[i].C.y;
                    u = sqrt(dx * dx + dy * dy);
                    if ( u > P[j].select ) skip_clump = 1;
                }
                else if( lens[i].z >= image->z )
                    skip_clump = 1;
    
                if ( ! skip_clump )
                {
                    if( lens[i].z != oldz )
                    {
    		            dls = distcosmo2(lens[i].z, image->z);
                        oldz = lens[i].z;
                    }
    
                    grad2 = e_grad2_pot(&image->C, i);
                    igrad2->a += grad2.a * dls;
                    igrad2->b += grad2.b * dls;
                    igrad2->c += grad2.c * dls;
                }
            }
    }
    
    // and those of the multiscale grid (G.nmsgrid -> G.nlens)
    // profile is precomputed so just multiply by b0
    if ( image->np_grad2a != NULL )
    {
        grad2 = e_grad2_np(image, np_b0);
        MA.a += grad2.a;
        MA.b += grad2.b;
        MA.c += grad2.c;
    }

    // add the potential of the not optimised clumps
    MA.a += igrad2->a;
    MA.b += igrad2->b;
    MA.c += igrad2->c;

    MA.d = MA.b;

    return (MA);

}

/* Return the laplacian of all the potentials, when optimising
 * in non-parametric mode for an image or arclet.
 *
 * The laplacian must have been initialised first with prep_non_param().
 */
static struct matrix e_grad2_np(struct galaxie *image, double * np_b0_thread)
{
//    const extern struct pot lens[];
    const extern struct g_grille  G;
    const extern double *np_b0;
    const double *np_b0_local;

    struct matrix MA;
    double a, b, c, *pgrad2;
    long int k;
//    long int k, n, startBin[100];
//    int nzbins, i_z;
//    double oldz, dls[100];

//    n = G.nlens - G.nmsgrid;
    a = b = c = 0.;
    if( np_b0_thread != NULL )
        np_b0_local = np_b0_thread;
    else
        np_b0_local = np_b0;

    // Cut the catalog of grid potentials in redshift bins
//    nzbins = 1;
//    startBin[0] = 0;
   
//    oldz = lens[G.nmsgrid].z; dls[0] = image->dl0s;
//    for ( k = G.nmsgrid; k < G.nlens; k++ )
//        if( oldz != lens[k].z ) 
//        {
//            startBin[nzbins] = k - G.nmsgrid;
//            oldz = lens[k].z;
//            dls[nzbins] = oldz < image->z ? distcosmo2(oldz, image->z): 0.;
//            nzbins ++;
//        }
//    startBin[nzbins] = n;

    pgrad2 = image->np_grad2a;
//    for ( i_z = 0; i_z < nzbins; i_z++ )
//        for ( k = startBin[i_z]; k < startBin[i_z+1]; k++ ) 
        for ( k = 0; k < G.nlens - G.nmsgrid; k++ )
            a += np_b0_local[k] * pgrad2[k];// * dls[i_z];

    pgrad2 = image->np_grad2b;
//    for ( i_z = 0; i_z < nzbins; i_z++ )
//        for ( k = startBin[i_z]; k < startBin[i_z+1]; k++ ) 
        for ( k = 0; k < G.nlens - G.nmsgrid; k++ )
            b += np_b0_local[k] * pgrad2[k];// * dls[i_z];

    pgrad2 = image->np_grad2c;
//    for ( i_z = 0; i_z < nzbins; i_z++ )
//        for ( k = startBin[i_z]; k < startBin[i_z+1]; k++ ) 
        for ( k = 0; k < G.nlens - G.nmsgrid; k++ )
            c += np_b0_local[k] * pgrad2[k];// * dls[i_z];

    MA.a = a; MA.b = MA.d = b; MA.c = c;

    return MA;
}

/***********************************************************/
struct matrix e_grad2_pot_ptr(const struct point *pi, const struct pot *ilens)
{
    const extern struct g_grille  G;
    struct point   R;  /*position of a lens relative to pi*/
    struct point   Q;  /*position of R rotated so that the lens is // to Xaxis*/
    double ep, em, ee; /*special potential ellipticities*/
    double t05, t10, t15, q, z, p, al, be; //potential temporary variables
    double X, Y, RR;
    struct polar   QP;  // for -1 potential
    struct matrix  grad2, g2, g05c, g05cut, g10c, g10cut, g15c, g15cut; //,g21,g20
    struct matrix  gsiemd;
    struct point   Qe;  // elliptical radius NFW/Sersic

    double  kap, tell;
    struct point gamma;

    g2.a = g2.b = g2.c = g2.d = 0.;

    // shortcut for N atoms
    if (ilens->b0 == 0. && ilens->type != 14)
        return g2;

    R.x = pi->x - ilens->C.x;
    R.y = pi->y - ilens->C.y;
    Q = rotation(R, ilens->theta);
    ep = 1. + ilens->epot;  /* a/(a+b) */
    em = 1. - ilens->epot; /* b/(a+b) */
    ee = 1. - ilens->epot*ilens->epot; /* 4ab/(a+b)^2 */

    switch (ilens->type)
    {
        case(1):
            z = ilens->b0 * pow(par1(Q.x, Q.y, ilens), -1.5);
            g2.a = z * ee * Q.y * Q.y;
            g2.b = g2.d = -z * ee * Q.x * Q.y;
            g2.c = z * ee * Q.x * Q.x;
            break;
        case(-1):
            QP = polxy(Q);
            g2.a = 0.;
            g2.b = g2.d = 0.;
            g2.c = ilens->b0 / QP.r / sqrt(1. - ilens->epot * 3.*cos(2.*QP.theta));
            g2 = rotmatrix(&g2, QP.theta);
            break;
        case(-2):
            mdcsiemd(Q.x, Q.y, ilens->epot, ilens->b0, &gsiemd);
            mdci05(Q.x, Q.y, ilens->epot, ilens->rcut, ilens->b0, &g05cut);

            g2.a = gsiemd.a - g05cut.a;
            g2.b = gsiemd.b - g05cut.b;
            g2.c = gsiemd.c - g05cut.c;
            g2.d = gsiemd.d - g05cut.d;

            break;
        case(3): // PL with core
            z = ilens->rc * ilens->rc;
            p = par2(Q.x, Q.y, ilens);
            q = 2 * ilens->alpha * ilens->b0 / ilens->rc * pow(p, ilens->alpha - 2.);
            g2.a = em * q * (p + 2.*em * (ilens->alpha - 1.) * Q.x * Q.x / z);
            g2.b = g2.d = 2.*ee * q * (ilens->alpha - 1.) * Q.x * Q.y / z;
            g2.c = ep * q * (p + 2.*ep * (ilens->alpha - 1.) * Q.y * Q.y / z);
            break;
        case(4):
            z = ilens->rc * ilens->rc;
            X = Q.x * Q.x / z;
            Y = Q.y * Q.y / z;
            RR = X + Y;
            q = ilens->b0 / ilens->rc / pow(1. + RR, 1.5);
            g2.a = q * (1. + Y - ilens->epot / 2. / (1. + RR) * (2. - X + 5.*Y - 3.*Y * (X - Y)));
            g2.b = g2.d = -q * Q.x * Q.y / z * (1. + ilens->epot * 3. / 2.*(X - Y) / (1. + RR));
            g2.c = q * (1. + X + ilens->epot / 2. / (1. + RR) * (2. + 5.*X - Y + 3.*X * (X - Y)));
            break;
        case(41):
            /* same as 4 but with a cut
            */
            z = ilens->rc * ilens->rc;
            X = Q.x * Q.x / z;
            Y = Q.y * Q.y / z;
            RR = X + Y;
            if (sqrt(RR) < ilens->rcut)
            {
                q = ilens->b0 / ilens->rc / pow(1. + RR, 1.5);
                g2.a = q * (1. + Y - ilens->epot / 2. / (1. + RR) * (2. - X + 5.*Y - 3.*Y * (X - Y)));
                g2.b = g2.d = -q * Q.x * Q.y / z * (1. + ilens->epot * 3. / 2.*(X - Y) / (1. + RR));
                g2.c = q * (1. + X + ilens->epot / 2. / (1. + RR) * (2. + 5.*X - Y + 3.*X * (X - Y)));

                z = 2.*ilens->psicut;
                g2.a -= z;
                g2.c -= z;
                break;
            }
            else
            {
                X = Q.x * Q.x;
                Y = Q.y * Q.y;
                RR = X + Y;
                z = ilens->psimcut / RR / RR;
                g2.a = z * (Y - X);
                g2.b = g2.d = -2.*z * Q.x * Q.y;
                g2.c = z * (X - Y);
            }
            break;
        case(5):    // not corrected for ellipticity
            QP = polxy(Q);
            z = ilens->rc * ilens->rc;
            RR = QP.r * QP.r / z;
            q = ilens->b0 / RR / ilens->rc;
            g2.c = q * (log(1. + RR) / 2. + ilens->epot * cos(2.*QP.theta)
                        * (log(1. + RR) / RR - (2. + RR) / (1. + RR)));
            g2.b = g2.d = q * ilens->epot * sin(2.*QP.theta) / 2.
                          *((3. + RR) / (1. + RR) - 3.*log(1. + RR) / RR);
            g2.a = q * ( (RR / (1. + RR) - log(1. + RR) / 2.)
                         + ilens->epot / 2.*cos(2.*QP.theta)
                         * (3.*log(1. + RR) / RR - (3. + 5.*RR) / (1. + RR) / (1. + RR)));
            g2 = rotmatrix(&g2, QP.theta);
            break;
        case(6):    // not corrected for ellipticity
            QP = polxy(Q);
            z = ilens->rc * ilens->rc;
            RR = QP.r * QP.r / z;
            q = 2.*ilens->b0 / ilens->rc;
            al = ilens->alpha;
            be = ilens->beta;
            g2.a = al * q * (1. + (2.*al - 1.) * RR) / pow(1. + RR, 2. - al) + ilens->epot * q
                   * (1. + (2. - 5.*be) * RR + (2.*be * be - 3.*be + 1.) * RR * RR) / pow(1. + RR, be + 2.)
                   * cos(2.*QP.theta);
            g2.b = g2.d = -q * ilens->epot * sin(2.*QP.theta)
                          * (1. + (1. - 2.*be) * RR) / pow(1. + RR, be + 1.);
            g2.c = al * q / pow(1. + RR, 1. - al) - ilens->epot * q * cos(2.*QP.theta)
                   * (1. + (1. + be) * RR) / pow(1. + RR, be + 1.);
            g2 = rotmatrix(&g2, QP.theta);
            break;
        case(7):    // not corrected for ellipticity
            if (ilens->b0 != 0.)
            {
                X = Q.x * Q.x;
                Y = Q.y * Q.y;
                RR = X + Y;
                q = (G.dx * G.dx + G.dy * G.dy) / 4.;
                /* if(RR>q) */
                if ((fabs(Q.x) > G.dx / 2.) || (fabs(Q.y) > G.dy / 2.))
                {
                    z = ilens->b0 / RR / RR;
                    g2.a = z * (Y - X);
                    g2.b = g2.d = -2.*z * Q.x * Q.y;
                    g2.c = z * (X - Y);
                }
                else
                {
                    z = ilens->b0 / q;
                    g2.a = z;
                    g2.b = g2.d = 0.;
                    g2.c = z;
                }
            }
            else
                g2.a = g2.b = g2.d = g2.c = 0.;
            break;

        case(8): /* PIEMD kovner */
            mdci05(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0, &g2);
            break;

        case(81): // trucated PIEMD kovner (HK model1)
            if ( ilens->epot > 2E-4 )
            {
                t05 = ilens->rcut / (ilens->rcut - ilens->rc);
                mdci05(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0, &g05c);
                mdci05(Q.x, Q.y, ilens->epot, ilens->rcut, ilens->b0, &g05cut);
                g2.a = t05 * (g05c.a - g05cut.a);
                g2.b = t05 * (g05c.b - g05cut.b);
                g2.c = t05 * (g05c.c - g05cut.c);
                g2.d = t05 * (g05c.d - g05cut.d);
            }
            else if ( (RR = Q.x * Q.x + Q.y * Q.y) > 0. )
            {
                // Circular dPIE Elliasdottir 2007 Eq A23 slighly modified for t05
                X = ilens->rc;
                Y = ilens->rcut;
                t05 = ilens->b0 * Y / (Y - X); // 1/u because t05/sqrt(u) and normalised Q/sqrt(u)
                z  = sqrt(RR + X * X) - X - sqrt(RR + Y * Y) + Y;  // R*dphi/dR
                X = RR / X;
                Y = RR / Y;
                p  = (1. - 1. / sqrt(1. + X / ilens->rc)) / X - (1. - 1. / sqrt(1. + Y / ilens->rcut)) / Y;  // d2phi/dR2
                X = Q.x * Q.x / RR;
                Y = Q.y * Q.y / RR;
                g2.a = t05 * (p * X + z * Y / RR);
                g2.c = t05 * (p * Y + z * X / RR);
                X = Q.x * Q.y / RR;
                g2.b = g2.d = t05 * (p * X - z * X / RR);
            }
            else
            {
                g2.a = g2.c = ilens->b0 / ilens->rc / 2.;
                g2.b = g2.d = 0.;
            }
            break;

        case(82): /* PIEMD kovner with a shallower central slope*/
            t05 = (ilens->rcut + ilens->rc) / ilens->rcut;
            mdci05(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0, &g05c);
            mdci05(Q.x, Q.y, ilens->epot, ilens->rcut, ilens->b0, &g05cut);
            g2.a = t05 * (g05c.a + g05cut.a);
            g2.b = t05 * (g05c.b + g05cut.b);
            g2.c = t05 * (g05c.c + g05cut.c);
            g2.d = t05 * (g05c.d + g05cut.d);
            break;

        case(83): /* EMD kovner 3/2 */
            g2 = mdci15(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0);
            break;

        case(84): /* EMD kovner I0.5c-I0.5cut + I1.5c */

            t05 = ilens->rcut / (ilens->rcut - ilens->rc);
            t15 = ilens->rcut / ilens->rc;
            g15c = mdci15(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0);
            mdci05(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0, &g05c);
            mdci05(Q.x, Q.y, ilens->epot, ilens->rcut, ilens->b0, &g05cut);

            q = ilens->alpha;
            g2.a = q * t15 * g15c.a + (1 - q) * t05 * (g05c.a - g05cut.a);
            g2.b = q * t15 * g15c.b + (1 - q) * t05 * (g05c.b - g05cut.b);
            g2.c = q * t15 * g15c.c + (1 - q) * t05 * (g05c.c - g05cut.c);
            g2.d = q * t15 * g15c.d + (1 - q) * t05 * (g05c.d - g05cut.d);
            break;

        case(85): /* EMD kovner 1: I1c*/

            g2 = mdci10(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0);
            break;

        case(86): /* EMD kovner 1: I1c-I1cut*/

            t10 = ilens->rc / (1. - ilens->rc * ilens->rc / ilens->rcut / ilens->rcut);
            g10c = mdci10(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0);
            g10cut = mdci10(Q.x, Q.y, ilens->epot, ilens->rcut, ilens->b0);
            g2.a = t10 * (g10c.a - g10cut.a);
            g2.b = t10 * (g10c.b - g10cut.b);
            g2.c = t10 * (g10c.c - g10cut.c);
            g2.d = t10 * (g10c.d - g10cut.d);
            break;

        case(87): /* EMD kovner 1: I1c-I1cut + I0.5c - I0.5cut */

            t05 = ilens->rcut / (ilens->rcut - ilens->rc);
            t10 = ilens->rc / (1. - ilens->rc * ilens->rc / ilens->rcut / ilens->rcut);
            mdci05(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0, &g05c);
            mdci05(Q.x, Q.y, ilens->epot, ilens->rcut, ilens->b0, &g05cut);
            g10c = mdci10(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0);
            g10cut = mdci10(Q.x, Q.y, ilens->epot, ilens->rcut, ilens->b0);
            q = ilens->alpha;
            g2.a = (1 - q) * t05 * (g05c.a - g05cut.a) + q * t10 * (g10c.a - g10cut.a);
            g2.b = (1 - q) * t05 * (g05c.b - g05cut.b) + q * t10 * (g10c.b - g10cut.b);
            g2.c = (1 - q) * t05 * (g05c.c - g05cut.c) + q * t10 * (g10c.c - g10cut.c);
            g2.d = (1 - q) * t05 * (g05c.d - g05cut.d) + q * t10 * (g10c.d - g10cut.d);
            break;

        case(88): /* EMD kovner 1: I1c-I1cut + I1.5cut */

            t10 = ilens->rc / (1. - ilens->rc * ilens->rc / ilens->rcut / ilens->rcut);
            t15 = ilens->rcut / ilens->rc;
            g10c = mdci10(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0);
            g10cut = mdci10(Q.x, Q.y, ilens->epot, ilens->rcut, ilens->b0);
            g15cut = mdci15(Q.x, Q.y, ilens->epot, ilens->rcut, ilens->b0);
            q = ilens->alpha;
            g2.a = (1 - q) * t15 * g15cut.a + q * t10 * (g10c.a - g10cut.a);
            g2.b = (1 - q) * t15 * g15cut.b + q * t10 * (g10c.b - g10cut.b);
            g2.c = (1 - q) * t15 * g15cut.c + q * t10 * (g10c.c - g10cut.c);
            g2.d = (1 - q) * t15 * g15cut.d + q * t10 * (g10c.d - g10cut.d);
            break;

        case(89): /* EMD kovner 1: I1c-I1cut + I0.5c - I0.5cut */

            t05 = 1. / (1. - 1. / ilens->beta);
            t10 = ilens->rc / (1. - 1. / ilens->beta / ilens->beta);
            mdci05(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0, &g05c);
            mdci05(Q.x, Q.y, ilens->epot, ilens->rc*ilens->beta, ilens->b0, &g05cut);
            g10c = mdci10(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0);
            g10cut = mdci10(Q.x, Q.y, ilens->epot, ilens->rc * ilens->beta, ilens->b0);
            q = ilens->alpha;
            g2.a = q * t05 * (g05c.a - g05cut.a) + (1. - q) * t10 * (g10c.a - g10cut.a);
            g2.b = q * t05 * (g05c.b - g05cut.b) + (1. - q) * t10 * (g10c.b - g10cut.b);
            g2.c = q * t05 * (g05c.c - g05cut.c) + (1. - q) * t10 * (g10c.c - g10cut.c);
            g2.d = q * t05 * (g05c.d - g05cut.d) + (1. - q) * t10 * (g10c.d - g10cut.d);
            break;

        case(9):
            g2.a = ilens->b0;
            g2.b = g2.d = 0.;
            g2.c = ilens->b0;
            break;

        case(10):
            grad2 = sp_grad2(*pi);
            break;

        case(11): /* 1/r4 mass */
            z = ilens->rc * ilens->rc;
            q = 2.*ilens->b0 / ilens->rc;
            X = Q.x * Q.x / z;
            Y = Q.y * Q.y / z;
            RR = 1. + (Q.x * Q.x + Q.y * Q.y) / z;
            g2.a = q * (1. - X + Y) / RR;
            g2.b = -2.*q * Q.x * Q.y / z / RR;
            g2.c = q * (1. + X - Y) / RR;
	    g2.d = g2.b;
            break;

        case(12): /* NFW */
            /*  Qe.x=Q.x-ilens->rc*ilens->emass/(1.-ilens->emass);
              Qe.y=Q.y;
              QP=polxy(Qe);
              kap=nfw_kappa(QP.r,ilens->rc,ilens->b0/2.);
              gam=nfw_gamma(QP.r,ilens->rc,ilens->b0/2.);
              g20.a=kap-gam;
              g20.b=g20.d=0.;
              g20.c=kap+gam;
              g20=rotmatrix(g20,QP.theta);
              Qe.x=Q.x+ilens->rc*ilens->emass/(1.-ilens->emass);
              Qe.y=Q.y;
              QP=polxy(Qe);
              kap=nfw_kappa(QP.r,ilens->rc,ilens->b0/2.);
              gam=nfw_gamma(QP.r,ilens->rc,ilens->b0/2.);
              g21.a=kap-gam;
              g21.b=g21.d=0.;
              g21.c=kap+gam;
              g21=rotmatrix(g21,QP.theta);
              g2.a=g21.a+g20.a;
              g2.b=g21.b+g20.b;
              g2.c=g21.c+g20.c;
              g2.d=g21.d+g20.d;
              */
            /*
            QP=polxy(Q);
            kap=nfw_kappa(QP.r,ilens->rc,ilens->b0);
            gam=nfw_gamma(QP.r,ilens->rc,ilens->b0);
            g2.a=kap-gam;
            g2.b=g2.d=0.;
            g2.c=kap+gam;
            g2=rotmatrix(g2,QP.theta);
            */
            Qe.x = sqrt(1. - ilens->emass) * Q.x;
            Qe.y = sqrt(1. + ilens->emass) * Q.y;
            QP = polxy(Qe);
            if ( QP.r == 0 ) QP.r = 1E-8;

            if ( ilens->alpha == 0. )
            {
                kap = nfw_kappa_eps(QP.r, ilens->rc, QP.theta, ilens->b0, ilens->emass);
                gamma = nfw_gamma_eps(QP.r, ilens->rc, QP.theta, ilens->b0, ilens->emass);
            }
            else
            {
                kap = nfwg_kappa_eps(QP.r, ilens->rc, QP.theta, ilens->b0, ilens->emass,
                                     ilens->alpha);
                gamma = nfwg_gamma_eps(QP.r, ilens->rc, QP.theta, ilens->b0, ilens->emass,
                                     ilens->alpha);
            }
//          }
//          else
//          {
//              kap = DBL_MAX;
//              gam = 0.5;
//          }
            /*
            kap=nfw_kappa(QP.r,ilens->rc,ilens->b0);
            gam=nfw_gamma(QP.r,ilens->rc,ilens->b0);
            */
            //QP = polxy(Qe);
            //double gam1 = ilens->emass * kapc + gamc * cos(2. * QP.theta);
            //double gam2 = -sqrt(1. - ilens->emass * ilens->emass) * gamc * sin(2. * QP.theta);
            g2.a = kap - gamma.x;
            g2.b = g2.d = gamma.y;
            g2.c = kap + gamma.x;
            //g2 = rotmatrix(&g2, QP.theta);
            break;
        case(121): // NFW triaxial
            tell = elli_tri(ilens);
			Qe.x=sqrt(1.-tell)*Q.x;
			Qe.y=sqrt(1.+tell)*Q.y;
			QP=polxy(Qe);
			if (ilens->alpha == 0.) // standard NFW
			{
   				kap=nfw_kappa_eps(QP.r,ilens->rc,QP.theta,ilens->b0,tell);
   				gamma=nfw_gamma_eps(QP.r,ilens->rc,QP.theta,ilens->b0,tell);
			}
            g2.a = kap - gamma.x;
            g2.b = g2.d = gamma.y;
            g2.c = kap + gamma.x;
            break;
        case(13): // Sersic
            Qe.x = sqrt(1. - ilens->emass) * Q.x;
            Qe.y = sqrt(1. + ilens->emass) * Q.y;
            QP = polxy(Qe);
            kap = sersic_kappa_eps(QP.r, ilens->rc, ilens->alpha, QP.theta, ilens->b0, ilens->emass);
            gamma = sersic_gamma_eps(QP.r, ilens->rc, ilens->alpha, QP.theta, ilens->b0, ilens->emass);
            g2.a = kap - gamma.x;
            g2.b = g2.d = gamma.y;
            g2.c = kap + gamma.x;
//            QP = polxy(Q);
//            g2 = rotmatrix(&g2, QP.theta);
            break;
        case(14):  // local shear
            g2.a = ilens->emass;
            g2.c = -ilens->emass;
            g2.b = g2.d = 0;
            break;
        case(15)://EInasto
	    Qe.x=sqrt(1.-ilens->emass)*Q.x;
	    Qe.y=sqrt(1.+ilens->emass)*Q.y;
	    QP=polxy(Qe);
	    kap=einasto_kappa_eps(QP.r,ilens->rc,ilens->alpha,QP.theta,ilens->b0,ilens->pmass,ilens->emass); //kappa elliptique
	    gamma=einasto_gamma_eps(QP.r,ilens->rc,ilens->alpha,QP.theta,ilens->b0,ilens->pmass,ilens->emass);
	    g2.a=kap-gamma.x;
	    g2.b=g2.d=gamma.y;
	    g2.c=kap+gamma.x;
	    //QP=polxy(Q);
	    //g2=rotmatrix(&g2,QP.theta);
	    break;
        case(16): // Hernquist model
            Qe.x = sqrt(1. - ilens->emass) * Q.x;
            Qe.y = sqrt(1. + ilens->emass) * Q.y;
            QP = polxy(Qe);
            kap = hern_kappa_eps(QP.r, ilens->rc, QP.theta, ilens->b0, ilens->emass);
            gamma = hern_gamma_eps(QP.r, ilens->rc, QP.theta, ilens->b0, ilens->emass);

            g2.a = kap - gamma.x;
            g2.b = g2.d = gamma.y;
            g2.c = kap + gamma.x;
            break;
        default:
            fprintf(stderr, "ERROR: profil type of clump %s unknown: %d \n",
                    ilens->n, ilens->type);
            exit(-1);
            break;
    };

    if (ilens->type != 10)
        grad2 = rotmatrix(&g2, ilens->theta);

    return grad2;
}

/***********************************************************/

struct matrix e_grad2_pot(const struct point *pi, long int i)
{
    const extern struct pot       lens[];

    const struct pot     *ilens = &lens[i];
    return e_grad2_pot_ptr(pi, ilens);
}



/***********************************************************/

static double par1(double x, double y, const struct pot *ilens)
{
    return( (1. - ilens->epot)*x*x + (1. + ilens->epot)*y*y );
}

/***********************************************************/

static double par2(double x, double y, const struct pot *ilens)
{
    double z;
    z = ilens->rc * ilens->rc;
    return( 1. + (1 - ilens->epot)*x*x / z + (1 + ilens->epot)*y*y / z );
}
