#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<string.h>


static double hypo(double x, double y, long int i);
static struct point e_grad_np(const struct galaxie *image, double *np_b0);

/****************************************************************/
/*      nom:        e_grad              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Return the gradient of the projected lens potential.
 *
 * !!! You have to multiply by dlsds to obtain the true gradient
 *
 * Global variables used (ALL global variables is constant here) : 
 * - lens, G
 * - in hypo() : lens
 * - in nfwg_dpl() : lens_table
 */
struct point e_grad_pot(const struct point *pi, long int i)
{
//  const extern  struct  g_mode      M;
    const extern struct g_grille G;
    const extern struct pot      lens[];
    const struct pot *ilens = &lens[i];
    double X, Y, R;
    double t, u, q, t05, t10, t15, z;
    double za, zb, zas, zbc, dpl_rad;
    struct point P, Q, Grad, g, pg; //,Qnfw
    struct polar QP;
    complex zis, zis_cut, ztot, ztot_cut;
    double tell;

    g.x = g.y = 0.;

    // shortcut for N atoms
    if (ilens->b0 == 0.)
        return g;
   
    /*positionning at the potential center*/
    P.x = pi->x - ilens->C.x;
    P.y = pi->y - ilens->C.y;

    /*rotation to the potential axes*/
    Q = rotation(P, ilens->theta);

    switch (ilens->type)
    {
        case(1):
            u = sqrt(hypo(Q.x, Q.y, i));
            g.x = ilens->b0 * (1 - ilens->epot) * Q.x / u;
            g.y = ilens->b0 * (1 + ilens->epot) * Q.y / u;
            break;
        case(-1):
            QP = polxy(Q);
            z = sqrt(2.*ilens->epot * 3.);
            za = z / sqrt(1 + ilens->epot * 3.);
            zb = z / sqrt(1 - ilens->epot * 3.);
            zas = -za * cos(QP.theta);
            zbc = zb * sin(QP.theta);
            pg.x = ilens->b0 / z * ( -cos(QP.theta) * asin(zas)
                                     + sin(QP.theta) * asinh(zbc) );
            pg.y = ilens->b0 / z * (sin(QP.theta) * asin(zas)
                                    + cos(QP.theta) * asinh(zbc) );
            g = rotation(pg, -QP.theta);
            break;
        case(-2):
            zis = csiemd(Q.x, Q.y, ilens->epot, ilens->b0);
            zis_cut = ci05(Q.x, Q.y, ilens->epot, ilens->rcut);
            g.x = zis.re - ilens->b0 * zis_cut.re;
            g.y = zis.im - ilens->b0 * zis_cut.im;
            break;

        case(3): // power law with core radius
            z = ilens->rc * ilens->rc;
            t = 2 * ilens->alpha * ilens->b0 / ilens->rc;
            u = pow(1. + (1 - ilens->epot) * Q.x * Q.x / z + (1 + ilens->epot) * Q.y * Q.y / z,
                    ilens->alpha - 1.);
            g.x = t * (1 - ilens->epot) * Q.x * u;
            g.y = t * (1 + ilens->epot) * Q.y * u;
            break;
        case(4):
            /* ici l'ellipticite est celle du potentiel */
            z = ilens->rc * ilens->rc;
            X = Q.x * Q.x / z;
            Y = Q.y * Q.y / z;
            R = X + Y;
            t = ilens->b0 / ilens->rc / sqrt(1. + R);
            g.x = t * Q.x * (1 - ilens->epot / 2. / (1. + R) * (2. + X + 3.*Y));
            g.y = t * Q.y * (1 + ilens->epot / 2. / (1. + R) * (2. + 3.*X + Y));
            break;
        case(41):
            /* same as 4 but with a troncated potentiel
            */
            z = ilens->rc * ilens->rc;
            X = Q.x * Q.x / z;
            Y = Q.y * Q.y / z;
            R = X + Y;
            if (sqrt(R) < ilens->rcut)
            {
                t = ilens->b0 / ilens->rc / sqrt(1. + R);
                g.x = t * Q.x * (1 - ilens->epot / 2. / (1. + R) * (2. + X + 3.*Y));
                g.y = t * Q.y * (1 + ilens->epot / 2. / (1. + R) * (2. + 3.*X + Y));
                /* cut */
                t = 2.*ilens->psicut;
                g.x -= t * Q.x;
                g.y -= t * Q.y;
            }
            else
            {
                u = Q.x * Q.x + Q.y * Q.y;
                g.x = ilens->psimcut * Q.x / u;
                g.y = ilens->psimcut * Q.y / u;
            }
            break;
        case(5):
            QP = polxy(Q);
            z = ilens->rc * ilens->rc;
            R = QP.r * QP.r / z;
            t = ilens->b0 * ilens->rc / 2. / QP.r;
            pg.x = t * (log(1. + R)
                        - ilens->epot * cos(2.*QP.theta) * (1. / (1. + R) - log(1. + R) / R));
            pg.y = ilens->epot * t * sin(2.*QP.theta) * (log(1. + R) / R - 1.);
            g = rotation(pg, -QP.theta);
            break;
        case(6):
            QP = polxy(Q);
            z = ilens->rc * ilens->rc;
            R = QP.r * QP.r / z;
            t = ilens->b0 / ilens->rc * QP.r;
            pg.x = 2.*ilens->alpha * t * pow(1. + R, ilens->alpha - 1.) +
                   2.*ilens->epot * t * (1. + (1. - ilens->beta) * R) /
                   pow(1. + R, 1. + ilens->beta) * cos(2.*QP.theta);
            pg.y = -2.*ilens->epot * t / pow(1. + R, ilens->beta) * sin(2.*QP.theta);
            g = rotation(pg, -QP.theta);
            break;
        case(7): /* point mass */
            if (ilens->b0 != 0.)
            {
                u = Q.x * Q.x + Q.y * Q.y;
                z = (G.dx * G.dx + G.dy * G.dy) / 4.;
                /* if(u>z) */
                if ( (fabs(Q.x) > G.dx / 2.) || (fabs(Q.y) > G.dy / 2.))
                {
                    g.x = ilens->b0 * Q.x / u;
                    g.y = ilens->b0 * Q.y / u;
                }
                else
                {
                    g.x = ilens->b0 * Q.x / z;
                    g.y = ilens->b0 * Q.y / z;
                }
            }
            else
                g.x = g.y = 0.;
            break;
        case(8): /* PIEMD kovner  I0.5 */
            zis = ci05(Q.x, Q.y, ilens->epot, ilens->rc);
            g.x = ilens->b0 * zis.re;
            g.y = ilens->b0 * zis.im;
            break;

        case(81): //PIEMD Kassiola & Kovner,1993 I0.5c-I0.5cut
            if ( ilens->epot > 2E-4 )
            {
                t05 = ilens->b0 * ilens->rcut / (ilens->rcut - ilens->rc);
                zis = ci05(Q.x, Q.y, ilens->epot, ilens->rc);
                zis_cut = ci05(Q.x, Q.y, ilens->epot, ilens->rcut);
                g.x = t05 * (zis.re - zis_cut.re);
                g.y = t05 * (zis.im - zis_cut.im);
            }
            else if ((u = Q.x * Q.x + Q.y * Q.y) > 0. )
            {
                // Circular dPIE Elliasdottir 2007 Eq A23 slighly modified for t05
                X = ilens->rc;
                Y = ilens->rcut;
                t05  = sqrt(u + X * X) - X - sqrt(u + Y * Y) + Y;  // Faster and equiv to Elliasdottir (see Golse PhD)
                t05 *= ilens->b0 * Y / (Y - X) / u; // 1/u because t05/sqrt(u) and normalised Q/sqrt(u)
                g.x = t05 * Q.x;
                g.y = t05 * Q.y;
            }
            else
            {
                g.x = 0.;
                g.y = 0.;
            }
            break;

        case(82): /* PIEMD kovner with shallow central slope I0.5c+I0.5cut*/

            t05 = ilens->b0 * (ilens->rcut + ilens->rc) / ilens->rcut;
            zis = ci05(Q.x, Q.y, ilens->epot, ilens->rc);
            zis_cut = ci05(Q.x, Q.y, ilens->epot, ilens->rcut);
            g.x = t05 * (zis.re + zis_cut.re);
            g.y = t05 * (zis.im + zis_cut.im);
            break;

        case(83): /* EMD kovner 3/2: I1.5c*/

            ztot = ci15(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0);
            g.x = ztot.re;
            g.y = ztot.im;
            break;

        case(84): /* EMD kovner isotherme I0.5c-I0.5cut + I1.5c */

            t05 = ilens->b0 * ilens->rcut / (ilens->rcut - ilens->rc);
            t15 = ilens->rcut / ilens->rc;
            zis = ci05(Q.x, Q.y, ilens->epot, ilens->rc);
            zis_cut = ci05(Q.x, Q.y, ilens->epot, ilens->rcut);
            ztot = ci15(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0);
            q = ilens->alpha;
            g.x = q * t15 * ztot.re + (1 - q) * t05 * (zis.re - zis_cut.re);
            g.y = q * t15 * ztot.im + (1 - q) * t05 * (zis.im - zis_cut.im);
            break;

        case(85): /* EMD kovner 1: I1c*/

            ztot = ci10(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0);
            g.x = ztot.re;
            g.y = ztot.im;
            break;

        case(86): /* EMD kovner 1: I1c-I1cut*/

            t10 = ilens->rc / (1. - ilens->rc * ilens->rc / ilens->rcut / ilens->rcut);
            ztot = ci10(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0);
            ztot_cut = ci10(Q.x, Q.y, ilens->epot, ilens->rcut, ilens->b0);
            g.x = t10 * (ztot.re - ztot_cut.re);
            g.y = t10 * (ztot.im - ztot_cut.im);
            break;

        case(87): /* EMD kovner  I1c-I1cut +I0.5c-I0.5cut */

            t05 = ilens->b0 / (1. - ilens->rc / ilens->rcut);
            t10 = ilens->rc / (1. - ilens->rc * ilens->rc / ilens->rcut / ilens->rcut);
            zis = ci05(Q.x, Q.y, ilens->epot, ilens->rc);
            zis_cut = ci05(Q.x, Q.y, ilens->epot, ilens->rcut);
            ztot = ci10(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0);
            ztot_cut = ci10(Q.x, Q.y, ilens->epot, ilens->rcut, ilens->b0);
            q = ilens->alpha;
            g.x = (1 - q) * t05 * (zis.re - zis_cut.re) + q * t10 * (ztot.re - ztot_cut.re);
            g.y = (1 - q) * t05 * (zis.im - zis_cut.im) + q * t10 * (ztot.im - ztot_cut.im);
            break;

        case(88): /* EMD kovner  I1c-I1cut +I1.5cut */

            t10 = ilens->rc / (1. - ilens->rc * ilens->rc / ilens->rcut / ilens->rcut);
            t15 = ilens->rcut / ilens->rc;
            zis_cut = ci15(Q.x, Q.y, ilens->epot, ilens->rcut, ilens->b0);
            ztot = ci10(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0);
            ztot_cut = ci10(Q.x, Q.y, ilens->epot, ilens->rcut, ilens->b0);
            q = ilens->alpha;
            g.x = (1 - q) * t15 * zis_cut.re + q * t10 * (ztot.re - ztot_cut.re);
            g.y = (1 - q) * t15 * zis_cut.im + q * t10 * (ztot.im - ztot_cut.im);
            break;

        case(89): /* EMD kovner  I1c-I1cut +I0.5c-I0.5cut */

            t05 = ilens->b0 / (1. - 1. / ilens->beta);
            t10 = ilens->rc / (1. - 1. / ilens->beta / ilens->beta);
            zis = ci05(Q.x, Q.y, ilens->epot, ilens->rc);
            zis_cut = ci05(Q.x, Q.y, ilens->epot, ilens->rc * ilens->beta);
            ztot = ci10(Q.x, Q.y, ilens->epot, ilens->rc, ilens->b0);
            ztot_cut = ci10(Q.x, Q.y, ilens->epot, ilens->rc * ilens->beta, ilens->b0);
            q = ilens->alpha;
            g.x = q * t05 * (zis.re - zis_cut.re) + (1 - q) * t10 * (ztot.re - ztot_cut.re);
            g.y = q * t05 * (zis.im - zis_cut.im) + (1 - q) * t10 * (ztot.im - ztot_cut.im);
            break;


        case(9): /* plan masse */
            t = ilens->b0;
            g.x = t * Q.x;
            g.y = t * Q.y;
            break;

        case(10): /* spline */
            Grad = sp_grad(*pi);
            break;

        case(11): /* 1/r4 mass */
            z = ilens->rc * ilens->rc;
            t = 2.*ilens->b0 / ilens->rc;
            R = 1. + (Q.x * Q.x + Q.y * Q.y) / z;
            g.x = t * Q.x / R;
            g.y = t * Q.y / R;
            break;

        case(12): /* NFW */
            /*     Qnfw.x=Q.x-ilens->rc*ilens->emass/(1.-ilens->emass);
              u=sqrt(Qnfw.x*Qnfw.x+Q.y*Q.y);
              g.x=nfw_dpl(u,ilens->rc,ilens->b0/2.)*Qnfw.x/u;
              g.y=nfw_dpl(u,ilens->rc,ilens->b0/2.)*Q.y/u;
                  Qnfw.x=Q.x+ilens->rc*ilens->emass/(1.-ilens->emass);
              u=sqrt(Qnfw.x*Qnfw.x+Q.y*Q.y);
              g.x+=nfw_dpl(u,ilens->rc,ilens->b0/2.)*Qnfw.x/u;
              g.y+=nfw_dpl(u,ilens->rc,ilens->b0/2.)*Q.y/u;
             */
            /*
             u=sqrt(Q.x*Q.x+Q.y*Q.y);
             g.x=nfw_dpl(u,ilens->rc,ilens->b0)*Q.x/u;
             g.y=nfw_dpl(u,ilens->rc,ilens->b0)*Q.y/u;
            */
            u = sqrt( (1. - ilens->emass) * Q.x * Q.x + (1. + ilens->emass) * Q.y * Q.y );
            if ( u == 0. )
            {
                g.x = g.y = 0.;
                return g;
            }
            if (ilens->alpha == 0.)
            {
                dpl_rad = nfw_dpl(u, ilens->rc, ilens->b0); // circular potential
                g.x = dpl_rad * (1. - ilens->emass) * Q.x / u;        // corrections to elliptical potential
                g.y = dpl_rad * (1. + ilens->emass) * Q.y / u;
            }
            else
            {
                dpl_rad = nfwg_dpl(u, ilens->rc, ilens->b0, ilens->alpha);
                g.x = dpl_rad * (1. - ilens->emass) * Q.x / u;
                g.y = dpl_rad * (1. + ilens->emass) * Q.y / u;
            }
            break;
        case(121): // NFW triaxial
            tell = elli_tri(ilens);
    		u = sqrt( (1.-tell)*Q.x*Q.x+(1.+tell)*Q.y*Q.y );
    		if( u == 0. ) // central point
    		{
    			g.x = g.y = 0.;
    			return g;
    		}	
    		if (ilens->alpha == 0.) // standard NFW
    		{
        		dpl_rad=nfw_dpl(u,ilens->rc,ilens->b0);	// circular potential
        		g.x=dpl_rad*(1.-tell)*Q.x/u;		// corrections to elliptical potential
        		g.y=dpl_rad*(1.+tell)*Q.y/u;
    		}
    		else  // generalized NFW
    		{
                fprintf(stderr, "ERROR: Triaxial generalized NFW not yet implemented\n");
                exit(1);
    		}
            break;
        case(13): // Sersic
            u = sqrt( (1. - ilens->emass) * Q.x * Q.x + (1. + ilens->emass) * Q.y * Q.y );
            dpl_rad = sersic_dpl(u, ilens->rc, ilens->alpha, ilens->b0); // circular potential
            g.x = dpl_rad * (1. - ilens->emass) * Q.x / u;        // corrections to elliptical potential
            g.y = dpl_rad * (1. + ilens->emass) * Q.y / u;
            break;
        case(14): // local shear
            g.x = ilens->emass * Q.x;
            g.y = -ilens->emass * Q.y;
            break;
        case(15): // Einasto model
	    u=sqrt( (1.-ilens->emass)*Q.x*Q.x+(1.+ilens->emass)*Q.y*Q.y );
	    dpl_rad=einasto_alpha(u,ilens->rc,ilens->alpha, ilens->pmass,ilens->b0);
	    g.x=dpl_rad*(1.-ilens->emass)*Q.x/u;
	    g.y=dpl_rad*(1.+ilens->emass)*Q.y/u;
	    break;
        case(16): // Hernquist model
            u = sqrt( (1. - ilens->emass) * Q.x * Q.x + (1. + ilens->emass) * Q.y * Q.y );
            dpl_rad = hern_dpl(u, ilens->rc, ilens->b0); // circular potential
            g.x = dpl_rad * (1. - ilens->emass) * Q.x / u;        // corrections to elliptical potential
            g.y = dpl_rad * (1. + ilens->emass) * Q.y / u;
            break;
        default:
            fprintf(stderr, "ERROR: profil type of clump %s unknown : %d\n",
                    ilens->n, ilens->type);
            exit(-1);
            break;
    };


    if (ilens->type != 10)
        Grad = rotation(g, -ilens->theta);

    return Grad;
}

struct point e_grad(const struct point *pi)
{
    const extern struct g_grille G;
    struct point Grad, grad;
    int i;

    /* NPRINTF(stderr,"e_grad\n"); */

    /* here elipticity is the potential one */

    grad.x = 0.;
    grad.y = 0.;

    /* for each lens*/
    for (i = 0; i < G.nlens; i++)
    {
        Grad = e_grad_pot(pi, i);
        grad.x += Grad.x;
        grad.y += Grad.y;
    };  /*end of for each lens*/

    return (grad);
}

/* Return the gradient at the position of an image contained in the
 * multi[ifamille][iimage] global variable.
 */
struct point e_grad_gal(struct galaxie *image, double *np_b0)
{
    const extern struct g_grille  G;
    const extern struct g_pot     P[NPOTFILE];
    const extern struct pot       lens[];
    const extern struct sigposStr sigposAs;

    struct point grad, Grad;
    struct point   *igrad = &image->grad;   // image's grad of not optimised clumps
    double dx, dy, u;
    int skip_clump; // if 1, skip the gradient computation for a clump
    long int i;
    int j;


    grad.x = grad.y = 0;

    // get the potential of the optimised clumps
    for ( i = 0; i < G.no_lens; i++ )
    {
        Grad = e_grad_pot(&image->C, i);
        grad.x += Grad.x;
        grad.y += Grad.y;
    }

    // and eventually those of the multiscale grid
    if ( image->np_grad != NULL )
    {
        Grad = e_grad_np(image, np_b0);
        grad.x += Grad.x;
        grad.y += Grad.y;
    }

    // if not defined, compute the gradient of the not optimised clumps
    if ( igrad->x == igrad->y )
    {
        igrad->x = 0.;
        igrad->y = 1e-10;  // to block images that skip to pots

        // Potentials no optimized but defined individually
        for ( i = G.no_lens; i < G.nplens[0]; i++ )
        {
            Grad = e_grad_pot(&image->C, i);
            igrad->x += Grad.x;
            igrad->y += Grad.y;
        }

        // Potentials in potfiles
        for ( j = 0; j < G.npot; j++ )
            for( i = G.nplens[j]; i < G.nplens[j+1]; i++ )
            {
                skip_clump = 0;  // do not skip the gradient computation (Default)
    
                if ( P[j].select == 1 )
                {
                    // test if the deflexion produced by this clump is detectable
                    // assuming a simple SIS potential in first approx
                    dx = image->C.x - lens[i].C.x;
                    dy = image->C.y - lens[i].C.y;
                    u = sqrt(dx * dx + dy * dy);
                    dx = lens[i].b0 * dx / u * image->dr;
                    dy = lens[i].b0 * dy / u * image->dr;
                    if ( dx*dx + dy*dy > sigposAs.max*sigposAs.max ) skip_clump = 1;
                }
                else if( P[j].select == 2 )
                {
                    // test on the distance to the image in arcsec
                    dx = image->C.x - lens[i].C.x;
                    dy = image->C.y - lens[i].C.y;
                    u = sqrt(dx * dx + dy * dy);
                    if ( u > 500. * lens[i].rc ) skip_clump = 1;
                }
                else if( lens[i].z >= image->z )
                    skip_clump = 1;
    
                if ( ! skip_clump )
                {
                    Grad = e_grad_pot(&image->C, i);
                    igrad->x += Grad.x;
                    igrad->y += Grad.y;
                }
            }
    }

    // add the potential of the not optimised clumps
    grad.x += igrad->x;
    grad.y += igrad->y;

    return grad;
}

/* Return the gradient of all the potentials, when optimising
 * in non-parametric mode.
 */
static struct point e_grad_np(const struct galaxie *image, double *np_b0_thread)
{
    const extern struct g_grille  G;
    const extern double *np_b0;
    const double *np_b0_local;

    struct point grad, *pgrad;
    long int k, n;
    double b0, gx, gy;

    gx = gy = 0.;
    if( np_b0_thread != NULL )
        np_b0_local = np_b0_thread;
    else
        np_b0_local = np_b0;

    n = G.nlens - G.nmsgrid;
    for ( k = 0; k < n; k++ )
    {
        pgrad = &image->np_grad[k];
        b0 = np_b0_local[k];

        gx += b0 * pgrad->x;
        gy += b0 * pgrad->y;
    }

    grad.x = gx;
    grad.y = gy;
    return grad;
}

/***********************************************************
 * Global variables used :
 * - lens
 */
static double hypo(double x, double y, long int i)
{
    const extern struct pot lens[];
    return( (1 - lens[i].epot)*x*x + (1 + lens[i].epot)*y*y );
}
