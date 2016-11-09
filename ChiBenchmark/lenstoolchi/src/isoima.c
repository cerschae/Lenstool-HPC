#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

static void rotelip(struct ellipse *el, double theta, struct ellipse *el_r);
static void formatrix(struct ellipse *eli, struct matrix *r);
static void mag(struct matrix *A, struct ellipse *ampli, struct matrix *B);


/****************************************************************/
/*      nom:        isoima              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Convert an ellipse from source to image plane or from
 * image to source plane according to the amplification matrix.
 *
 * When passing from image to source plane, ampli is the eigenvalues of the A^-1 matrix
 * ( 1-k+gamma, 1-k-gamma ) and theta is the magnification axis angle theta_pot.
 *
 * When passing from the source to the image plane, ampli is the eigenvalues of the A matrix
 * ( 1 / (1-k+gamma), 1 / (1-k-gamma) ) and theta is ALSO the magnification axis angle theta_pot.
 *
 * ampli(IN) : not modified
 * Global variables used :
 * - none
 */
void isoima(struct ellipse *es, struct ellipse *ampli, struct ellipse *ei)
{
    struct  ellipse esm, eim;
    struct  matrix  S, I;

    /* on se place dans le repere de magnification*/

    rotelip(es, ampli->theta, &esm);

    /* la matrice des fij, c'est le carre des longueur !!! */

    esm.a = esm.a * esm.a;
    esm.b = esm.b * esm.b;

    /* on determine la matrice correspondante f_i ou f_s (eq 2.92 JPK PhD thesis) */

    formatrix(&esm, &S);

    /* on determine la matrice image */

    mag(&S, ampli, &I);

    /* on determine l'ellipsoide correspondant par diagonalisation */

    eim = formeli(I.a, I.b, I.c);

    /* l'ellipse c'est la racine carre de l'ellipsoide */

    eim.a = sqrt(eim.a);
    eim.b = sqrt(eim.b);

    /* on retourne  dans le repere propre*/

    rotelip(&eim, -ampli->theta, ei);
}
/* Change the ellipse angle from its original reference frame to a new reference frame
 * with the same origin but rotated by an angle theta.
 */
static void rotelip(struct ellipse *el, double theta, struct ellipse *el_r)
{
    el_r->a = el->a;
    el_r->b = el->b;
    el_r->theta = el->theta - theta;
}

/* Build a matrix form of the ellipse. The procedure works in both the source and the image plane.
 * theta_i = 0 means that the ellipse major axis is aligned with the magnification direction
 * theta_pot.
 */
static void formatrix(struct ellipse *eli, struct matrix *r)
{
    struct matrix m;
    m.a = eli->a;    // start with fij = | a^2  0  |
    m.c = eli->b;    //                  |  0  b^2 |
    m.b = m.d = 0.;
    *r = rotmatrix(&m, eli->theta);   // and apply the rotation f = R(theta_i) fij R(-theta_i)
    // if theta_i = 0, f = fij
}

static void mag(struct matrix *A, struct ellipse *ampli, struct matrix *B)
{
    double alpha, beta;
    alpha = ampli->a;   // streching or shrinking factor along major magnification axis
    beta = ampli->b;    // streching or shrinking factor along minor magnification axis

    B->a = alpha * alpha * A->a;      // see eq 2.109 of JPK PhD thesis
    B->b = B->d = alpha * beta * A->b;  // where alpha is 1-K+gamma and beta is 1-K-gamma when passing
    B->c = beta * beta * A->c;  // from image to source plane.
}
