#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


static void petitA(const struct bitriplet *E, const struct bitriplet *I, struct bitriplet *petit);
static void petitB(const struct bitriplet *E, const struct bitriplet *I, struct bitriplet *petit);
static void petitC(const struct bitriplet *E, const struct bitriplet *I, struct bitriplet *petit);
static void zoom(const struct triplet *T, struct triplet *Z);
static void zoom1(const struct triplet *T, struct triplet *Z);
static void P1(const struct triplet *T, struct triplet *P);
static void P2(const struct triplet *T, struct triplet *P);


/****************************************************************/
/*      nom:        inthere             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Return the small couple of triangles of size 1/4 of the E.i triangle,
 * (ie. the size of the I.i triangle) in which P is located.
 * This small triangle can be inside E or a small triangle outside.
 *
 *          Ea                                          Rb1-- Ea -- Rc1
 *         /  \                                           \  /  \  /
 *       Ic -- Ib   the returned triangles can be    Ra1-- Ic -- Ib -- Ra2
 *      /  \  /  \                                     \  /  \  /  \  /
 *    Eb -- Ia -- Ec                                    Eb -- Ia -- Ec
 *                                                        \  /  \  /
 *                                                         Rc2   Rb2
 * Parameters :
 * - E : 2 triangles (simulated images and corresponding sources)
 * - I : 2 triangles inside the corresponding (image, source) 2 E triangles
 * - P : barycenter of the original sources
 * - dlsds : Lens efficiency at the source redshift
 *
 * Global variables used :
 * - in e_transform() : G, lens, lens_table
 * - in e_dpl() : G, lens, lens_table
 */
void e_inthere( const struct bitriplet *E,
                const struct bitriplet *I,
                const struct point *P,
                double dlsds, struct bitriplet *res)
{
    struct triplet tpl; // temporary triplet
    double sa, sb, sc, sAb, sAc, sBa, sBc, sCa, sCb;

    // in the source plane triangle I->s
    sc = signe(&I->s, P);  // >0 if P and Ic corner are on the same side?
    P1(&I->s, &tpl);  // rotate I->s anticlockwise
    sb = signe(&tpl, P);   // >0 if P and Ib corner are on the same side?
    P2(&I->s, &tpl);    // rotate I->s clockwise
    sa = signe(&tpl, P);   // >0 if P and Ib corner are on the same side?

    /*sa, sb and sc are >0 if P is in the small triangle I->s*/
    if ( sa >= 0. && sb >= 0. && sc >= 0. )
    {
        //return(I)
        res->i = I->i;
        res->s = I->s;
        return;
    }
    // P is not in the small triangle I->s
    else if ( sa < 0 )
        // P is in front of Ia on the other side of the [Ib,Ic] segment of I->s triangle*/
    {
        petitA(E, I, res);  // res is the small triangle (Ea,Ic,Ib) (both image&source plane)
        P1(&res->s, &tpl); // rotate res->s anti-clockwise
        sAb = signe(&tpl, P); // P on the same side as Ic?
        sAc = signe(&res->s, P); // P on the same side as Ib?

        if ( sAc >= 0. && sAb >= 0. )
            //P is inside res small triangle (Ea,Ic,Ib)
        {
            zoom1(&res->i, &tpl); // tlp is 2x larger than res in the image plane*/
            res->i = tpl;
            e_transform(&res->i, dlsds, &res->s);
            return;
        }
        // P is not in the small triangle res=(Ea,Ic,Ib)... so ouside of the E->s triangle
        else if ( sAc < 0 )
            // P is in front of Ib
        {
            // res becomes the (Ea,Ib,Rc1) triangle in image plane
            res->i.c.x = res->i.a.x + res->i.b.x - res->i.c.x;
            res->i.c.y = res->i.a.y + res->i.b.y - res->i.c.y;
            zoom(&res->i, &tpl); //... 1.5x larger
            res->i = tpl;
            //e_dpl(&res->s.c,dlsds,&res->s.c);
            e_transform(&res->i, dlsds, &res->s); //... and in source plane
            return;
        }
        else
            // P is in front of Ic
        {
            // res becomes the (Ea,Rb1,Ic) triangle in image plane
            res->i.b.x = res->i.a.x + res->i.c.x - res->i.b.x;
            res->i.b.y = res->i.a.y + res->i.c.y - res->i.b.y;
            zoom(&res->i, &tpl); // 1.5x larger
            res->i = tpl;
            //e_dpl(&res->s.b,dlsds,&res->s.b);
            e_transform(&res->i, dlsds, &res->s); //... and in source plane
            return;
        };
    }
    // P is still not in the small triangle I->s
    else if (sb < 0)
        // P is in front of Ib on the other side of the [Ic,Ia] segment of I->s triangle*/
    {
        petitB(E, I, res);    // res is the small triangle (Ia,Eb,Ic) (both image&source plane)
        sBc = signe(&res->s, P); // >0 if P on the same side as Ic?
        P2(&res->s, &tpl);    // rotate res->s clockwise
        sBa = signe(&tpl, P); // >0 if P on the same side as Ia?
        if ( sBc >= 0. && sBa >= 0. )
            //P is inside res small triangle (Ia,Eb,Ic)
        {
            zoom1(&res->i, &tpl);// tlp is 2x larger than res in the image plane*/
            res->i = tpl;
            e_transform(&res->i, dlsds, &res->s);
            return;
        }
        else if (sBc < 0)
        {
            // res becomes the (Rc2,Eb,Ia) triangle in image plane
            res->i.c.x = res->i.a.x + res->i.b.x - res->i.c.x;
            res->i.c.y = res->i.a.y + res->i.b.y - res->i.c.y;
            zoom(&res->i, &tpl);
            res->i = tpl;
            e_transform(&res->i, dlsds, &res->s);
            return;
        }
        else if (sBa < 0)
        {
            // res becomes the (Ra1,Eb,Ic) triangle in image plane
            res->i.a.x = res->i.c.x + res->i.b.x - res->i.a.x;
            res->i.a.y = res->i.c.y + res->i.b.y - res->i.a.y;
            zoom(&res->i, &tpl);
            res->i = tpl;
            e_transform(&res->i, dlsds, &res->s);
            return;
        };
    }
    // P is still not in the small triangle I->s
    else
        // P is in front of Ic on the other side of the [Ib,Ia] segment of I->s triangle*/
    {
        petitC(E, I, res);  // res is the small triangle (Ia,Ib,Ec) (both image&source plane)
        P2(&res->s, &tpl);
        sCa = signe(&tpl, P);  // P on the same side as Ia?
        P2(&res->s, &tpl);
        sCb = signe(&tpl, P);  // P on the same side as Ib?

        if ( sCa >= 0. && sCb >= 0. )
            //P is inside res small triangle (Ia,Ib,Ec)
        {
            zoom1(&res->i, &tpl);
            res->i = tpl;
            e_transform(&res->i, dlsds, &res->s);
            return;
        }
        else if ( sCa < 0 )
        {
            // res becomes the (Ra2,Ib,Ec) triangle in image plane
            res->i.a.x = res->i.c.x + res->i.b.x - res->i.a.x;
            res->i.a.y = res->i.c.y + res->i.b.y - res->i.a.y;
            zoom(&res->i, &tpl);
            res->i = tpl;
            e_transform(&res->i, dlsds, &res->s);
            return;
        }
        else
        {
            // res becomes the (Ia,Rb2,Ec) triangle in image plane
            res->i.b.x = res->i.c.x + res->i.a.x - res->i.b.x;
            res->i.b.y = res->i.c.y + res->i.a.y - res->i.b.y;
            zoom(&res->i, &tpl);
            res->i = tpl;
            e_transform(&res->i, dlsds, &res->s);
            return;
        }
    }
}

/**************************************************************/
/* Return 2 triangles (source and image planes) corresponding
 * for each plane to a small triangle which corners are [aE,Ib,Ic].
 *
 * Parameters :
 * - E : big triangle
 * - I : small triangle inscribed in the big E triangle
 *
 * Global variables used :
 * - none
 */
static void petitA(const struct bitriplet *E, const struct bitriplet *I, struct bitriplet *petit)
{
    petit->i.a = E->i.a;
    petit->i.b = I->i.b;
    petit->i.c = I->i.c;
    petit->s.a = E->s.a;
    petit->s.b = I->s.b;
    petit->s.c = I->s.c;
}

/**************************************************************/
/* Return 2 triangles (source and image planes) corresponding
 * for each plane to a small triangle which corners are [Ia,Eb,Ic].
 *
 * Parameters :
 * - E : big triangle
 * - I : small triangle inscribed in the big E triangle
 *
 * Global variables used :
 * - none
 */
static void petitB(const struct bitriplet *E, const struct bitriplet *I, struct bitriplet *petit)
{
    petit->i.a = I->i.a;
    petit->i.b = E->i.b;
    petit->i.c = I->i.c;
    petit->s.a = I->s.a;
    petit->s.b = E->s.b;
    petit->s.c = I->s.c;
}

/**************************************************************/
/* Return 2 triangles (source and image planes) corresponding
 * for each plane to a small triangle which corners are [Ia,Ib,Ec].
 *
 * Parameters :
 * - E : big triangle
 * - I : small triangle inscribed in the big E triangle
 *
 * Global variables used :
 * - none
 */
static void petitC(const struct bitriplet *E, const struct bitriplet *I, struct bitriplet *petit)
{
    petit->i.a = I->i.a;
    petit->i.b = I->i.b;
    petit->i.c = E->i.c;
    petit->s.a = I->s.a;
    petit->s.b = I->s.b;
    petit->s.c = E->s.c;
}

/**************************************************************/
/* Return a triangle with almost the same center as T but with
 * its sides multiplied by 2 and its surface by 4.
 *
 * Global variables used :
 * - none
 */
static void zoom(const struct triplet *T, struct triplet *Z)
{
    Z->a.x = ( 5.*T->a.x - T->b.x - T->c.x ) / 3.;
    Z->b.x = ( 5.*T->b.x - T->a.x - T->c.x ) / 3.;
    Z->c.x = ( 5.*T->c.x - T->b.x - T->a.x ) / 3.;
    Z->a.y = ( 5.*T->a.y - T->b.y - T->c.y ) / 3.;
    Z->b.y = ( 5.*T->b.y - T->a.y - T->c.y ) / 3.;
    Z->c.y = ( 5.*T->c.y - T->b.y - T->a.y ) / 3.;
}

/**************************************************************/
/* Return a triangle with almost the same center as T but with
 * its sides multiplied by 3/2 and its surface by 9/4.
 *
 * Global variables used :
 * - none
 */
static void zoom1(const struct triplet *T, struct triplet *Z)
{
    Z->a.x = ( 8.*T->a.x - T->b.x - T->c.x ) / 6.;
    Z->b.x = ( 8.*T->b.x - T->a.x - T->c.x ) / 6.;
    Z->c.x = ( 8.*T->c.x - T->b.x - T->a.x ) / 6.;
    Z->a.y = ( 8.*T->a.y - T->b.y - T->c.y ) / 6.;
    Z->b.y = ( 8.*T->b.y - T->a.y - T->c.y ) / 6.;
    Z->c.y = ( 8.*T->c.y - T->b.y - T->a.y ) / 6.;
}

/* Return the T triangle with anticlockwise rotated corners.
 *
 * Global variables used :
 * - none
 */
static void P1(const struct triplet *T, struct triplet *P)
{
    P->b = T->a;
    P->c = T->b;
    P->a = T->c;
}


/*********************************************************************
 * Return the T triangle with clockwise rotated corners.
 *
 * Global variables used :
 * - none
 */
static void P2(const struct triplet *T, struct triplet *P)
{
    P->c = T->a;
    P->a = T->b;
    P->b = T->c;
}
/*********************************************************************/

