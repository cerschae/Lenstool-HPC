#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<errors.h>

/****************************************************************/
/*      nom:        critic              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            *
 *
 * Compute the critic and caustic lines.
 *
 * Output : tangeant[ntline] and radial[nrline]
 *
 * Global variables used :
 * - G, CL, F, lens, radial, tangent, nrline, ntline, flagr, flagt
 * - in dratio() : C
 * - in chgsigne() : G, lens
 * - in zeroamp() : G, lens, lens_table
 * - in e_dpl() : G, lens, lens_table
 * - in next() : G, lens
 * - in follow() : CL, radial, nrline, flagr, G, lens, lens_table
 * - in followi() : CL, tangent, ntline, flagt, G, lens, lens_table
 */



static void snake(int verbose);
static void marchingSquares(int verbose);

void    criticnew(int verbose)
{
    extern struct   g_cline CL;

    if ( !strcmp(CL.algorithm, "SNAKE") )
        snake(verbose);
    else
        marchingSquares(verbose);

}

/* Variables definitions for the marchingSquare algorithm */
static double limitLow;     // smaller area of a square
static double limitHigh;    // minimal area of a square
static double z_cl;     // redshift of the current CL.nplan[i]
static double dos;              // distance to current CL.nplan[i]
static double dl0s;             // distance between lens[0] and CL.nplan[i]
static double dlsds_cl;     // dlsds ratio of dl0s/dos


/* Compute the amplification at the 4 cornes of the square and
 * return the corresponding marching square type.
 *
 * The amplification is always computed at the center of the square.
 *
 * The points are ccw oriented and the 1st point coordinates are xmin, ymin.
 *
 * If squareNb = 4, force to compute the 4 corners.
 *
 * Parameters :
 * - x, y, width, height : position and size of the square
 * - parentType : marching square type for the parent square (see: cutSquare() )
 * - squareNb : index of the current child square in the parent square
 * - ampli (IN/OUT) : array of double[5] that contains the 5 computed
 *                    amplification or 0. In input, it contains the parent values
 */
static unsigned char getSquareType(
    double x, double y,
    double width, double height,
    unsigned char parentType, char squareNb, double *ampli )
{
    struct point a;
    unsigned char type;

    // initialize type with parentType and squareNb
    if ( squareNb == 0 )
        type = (parentType & 1) | (parentType & 128 ? 8 : 0);
    else if ( squareNb == 1 )
        type = (parentType & 2) | (parentType & 128 ? 4 : 0);
    else if ( squareNb == 2 )
        type = (parentType & 8) | (parentType & 128 ? 1 : 0);
    else if ( squareNb == 3 )
        type = (parentType & 4) | (parentType & 128 ? 2 : 0);
    else
        type = 0;

    // compute only 3 amplifications per square
    if ( squareNb == 1 || squareNb == 3 || squareNb == 4)
    {
        a.x = x;
        a.y = y;
        ampli[0] = e_amp(&a, dl0s, dos, z_cl);
        type |=  ampli[0] > 0 ? 1 : 0;
    }
    if ( squareNb == 0 || squareNb == 2 || squareNb == 4 )
    {
        a.x = x + width;
        a.y = y;
        ampli[1] = e_amp(&a, dl0s, dos, z_cl);
        type |=  ampli[1] > 0 ? 2 : 0;
    }
    if ( squareNb == 0 || squareNb == 2 || squareNb == 4 )
    {
        a.x = x;
        a.y = y + height ;
        ampli[3] = e_amp(&a, dl0s, dos, z_cl);
        type |=  ampli[3] > 0 ? 4 : 0;
    }
    if ( squareNb == 1 || squareNb == 3 || squareNb == 4 )
    {
        a.x = x + width;
        a.y = y + height ;
        ampli[2] = e_amp(&a, dl0s, dos, z_cl);
        type |=  ampli[2] > 0 ? 8 : 0;
    }

    // center of the square... always computed
    a.x = x + width / 2.;
    a.y = y + height / 2.;
    ampli[4] = e_amp(&a, dl0s, dos, z_cl);
    type |=   ampli[4] > 0 ? 128 : 0;

    return type;
}

/* Append a segment to the tangent/radial arrays depending on the sign of
 * the parities on each side of the segment.
 * Parameters :
 * - x1, y1 : head of the segment
 * - x2, y2 : end of the segment
 * - signe_flag : 1 for radial line, 2 for tangeant line (see : chgsigne() )
 */
static void appendSegment( double x1, double y1, double x2, double y2, int signe_flag )
{
    extern struct biline radial[NMAX], tangent[NMAX];
    extern int nrline, ntline, flagr, flagt;

    if ( signe_flag == 1 )
    {
        radial[nrline].i = flagr;
        radial[nrline].I.x = x1;
        radial[nrline].I.y = y1;
        e_dpl(&radial[nrline].I, dlsds_cl, &radial[nrline].S);
        nrline++;

        radial[nrline].i = flagr;
        radial[nrline].I.x = x2;
        radial[nrline].I.y = y2;
        e_dpl(&radial[nrline].I, dlsds_cl, &radial[nrline].S);
        nrline++;
        flagr++;
    }
    else if (signe_flag == 2)
    {
        tangent[ntline].i = flagt;
        tangent[ntline].I.x = x1;
        tangent[ntline].I.y = y1;
        e_dpl(&tangent[ntline].I, dlsds_cl, &tangent[ntline].S);
        ntline++;

        tangent[ntline].i = flagt;
        tangent[ntline].I.x = x2;
        tangent[ntline].I.y = y2;
        e_dpl(&tangent[ntline].I, dlsds_cl, &tangent[ntline].S);
        ntline++;
        flagt++;
    }
#ifdef DEBUG
    else
        fprintf( stderr, "ERROR with chsigne!!!\n" );
#endif
    if( ntline > NMAX || nrline > NMAX )
    {
        fprintf(stderr, "Error: Too many critical line segments. Maximum set to NMAX = %d\n", NMAX);
        exit(E_NMAX_REACHED);
    }
}

/* Given a square, compute its type with getSquareType() function.
 *
 * If the square has at least one corner with a sign different from
 * the other corners, then the square is divided in 4 squares and
 * cutSquare is called for each of the 4 children.
 *
 * Stop when the square size is lower than the limitLow global variable.
 *
 * Parameters :
 * - x, y, width, height : position and size of the square
 * - parentType : marching square type for the parent square (see: cutSquare() )
 * - squareNb : index of the current child square in the parent square
 * - ampli (IN) : array of double[5] that contains the 5 computed
 *                    amplification or 0. In input, it contains the parent values
 */
static  void cutSquare(
    double x, double y, double width, double height,
    unsigned char parentType, int squareNb, double *parentAmpli )
{
#ifdef DEBUG
    extern struct g_mode M;
    FILE * dbg;
#endif

    struct point A, B;  //check points for radial/tangeantial line detection
    double sizeSquare;
    unsigned char type;
    double tmp1, tmp2;  //segment coordinates for the smallest square
    double ampli[5]; // amplifications for this square

#ifdef DEBUG
    dbg = fopen( "cidbg.reg", "a" );
    fprintf( dbg, "fk5;box(%lf,%lf,%lf\",%lf\") # color=green\n",
             M.ref_ra + (x + width / 2.) / -3600 / cos(M.ref_dec*DTR),
             M.ref_dec + (y + height / 2.) / 3600,
             width, height );
    fclose(dbg);
#endif

    sizeSquare = width > height ? width : height;

    // initialise the amplifications of this square
    ampli[squareNb] = parentAmpli[squareNb];
    ampli[(squareNb + 2) % 4] = parentAmpli[4];


    type = getSquareType( x, y, width, height, parentType, squareNb, ampli );

    if ( sizeSquare > limitLow )
    {
        // divide the square in 4 children
        if ( sizeSquare > limitHigh || ( type != 0 && type != 15 && type != 143) )
        {
            width /= 2.;
            height /= 2.;
            // square 1
            cutSquare( x, y, width, height, type, 0, ampli);
            // square 2
            cutSquare( x + width, y, width, height, type, 1, ampli);
            // square 3
            cutSquare( x + width, y + height, width, height, type, 2, ampli);
            //square 4
            cutSquare( x, y + width, width, height, type, 3, ampli);
        }
    }
    else
    {
        // Remove the detection of the central point (bit 7 = 128)
        if ( type != 137 && type != 134 ) type &= 127;

        // append a square marching square edge to the tangent/radial list
        switch ( type )
        {
            case 1:
                tmp1 = -ampli[0] / (ampli[3] - ampli[0]) * height;
                tmp2 = -ampli[0] / (ampli[1] - ampli[0]) * width;
                A.x = x;
                A.y = y;
                B.x = x + width;
                B.y = y;
                appendSegment(x, y + tmp1, x + tmp2, y, chsigne(A, B, dl0s, dos, z_cl) );
                break;
            case 2:
                tmp1 = -ampli[0] / (ampli[1] - ampli[0]) * width;
                tmp2 = -ampli[1] / (ampli[2] - ampli[1]) * height;
                A.x = x + width;
                A.y = y;
                B.x = x;
                B.y = y;
                appendSegment(x + tmp1, y, x + width, y + tmp2, chsigne(A, B, dl0s, dos, z_cl) );
                break;
            case 3:
                tmp1 = -ampli[0] / (ampli[3] - ampli[0]) * height;
                tmp2 = -ampli[1] / (ampli[2] - ampli[1]) * height;
                A.x = x;
                A.y = y;
                B.x = x;
                B.y = y + height;
                appendSegment(x, y + tmp1, x + width, y + tmp2, chsigne(A, B, dl0s, dos, z_cl) );
                break;
            case 4:
                tmp1 = -ampli[0] / (ampli[3] - ampli[0]) * height;
                tmp2 = -ampli[3] / (ampli[2] - ampli[3]) * width;
                A.x = x;
                A. y = y + height;
                B.x = x;
                B.y = y;
                appendSegment(x, y + tmp1, x + tmp2, y + height, chsigne(A, B, dl0s, dos, z_cl) );
                break;
            case 5:
                tmp1 = -ampli[0] / (ampli[1] - ampli[0]) * width;
                tmp2 = -ampli[3] / (ampli[2] - ampli[3]) * width;
                A.x = x;
                A.y = y;
                B.x = x + width;
                B.y = y;
                appendSegment(x + tmp1, y, x + tmp2, y + height, chsigne(A, B, dl0s, dos, z_cl) );
                break;
            case 6:
                tmp1 = -ampli[0] / (ampli[3] - ampli[0]) * height;
                tmp2 = -ampli[3] / (ampli[2] - ampli[3]) * width;
                A.x = x;
                A.y = y + height;
                B.x = x + width / 2.;
                B.y = y + height / 2.;
                appendSegment( x, y + tmp1, x + tmp2, y + height, chsigne(A, B, dl0s, dos, z_cl) );
                tmp1 = -ampli[0] / (ampli[1] - ampli[0]) * width;
                tmp2 = -ampli[1] / (ampli[2] - ampli[1]) * height;
                A.x = x + width;
                A.y = y;
                appendSegment( x + tmp1, y, x + width, y + tmp2, chsigne(A, B, dl0s, dos, z_cl) );
                break;
            case 134:
                tmp1 = -ampli[0] / (ampli[3] - ampli[0]) * height;
                tmp2 = -ampli[0] / (ampli[1] - ampli[0]) * width;
                A.x = x;
                A.y = y;
                B.x = x + width / 2.;
                B.y = y + height / 2.;
                appendSegment( x, y + tmp1, x + tmp2, y, chsigne(A, B, dl0s, dos, z_cl) );
                tmp1 = -ampli[1] / (ampli[2] - ampli[1]) * height;
                tmp2 = -ampli[3] / (ampli[2] - ampli[3]) * width;
                A.x = x + width;
                A.y = y + height;
                appendSegment( x + width, y + tmp1, x + tmp2, y + height, chsigne(A, B, dl0s, dos, z_cl) );
                break;
            case 7:
                tmp1 = -ampli[1] / (ampli[2] - ampli[1]) * height;
                tmp2 = -ampli[3] / (ampli[2] - ampli[3]) * width;
                A.x = x;
                A.y = y;
                B.x = x + width;
                B.y = y + height;
                appendSegment( x + width, y + tmp1, x + tmp2, y + height, chsigne(A, B, dl0s, dos, z_cl) );
                break;
            case 8:
                tmp1 = -ampli[3] / (ampli[2] - ampli[3]) * width;
                tmp2 = -ampli[1] / (ampli[2] - ampli[1]) * height;
                A.x = x;
                A.y = y;
                B.x = x + width;
                B.y = y + height;
                appendSegment( x + tmp1, y + height, x + width, y + tmp2, chsigne(A, B, dl0s, dos, z_cl) );
                break;
            case 9:
                tmp1 = -ampli[0] / (ampli[3] - ampli[0]) * height;
                tmp2 = -ampli[0] / (ampli[1] - ampli[0]) * width;
                A.x = x;
                A.y = y;
                B.x = x + width / 2.;
                B.y = y + height / 2.;
                appendSegment( x, y + tmp1, x + tmp2, y, chsigne(A, B, dl0s, dos, z_cl) );
                tmp1 = -ampli[1] / (ampli[2] - ampli[1]) * height;
                tmp2 = -ampli[3] / (ampli[2] - ampli[3]) * width;
                A.x = x + width;
                A.y = y + height;
                appendSegment( x + width, y + tmp1, x + tmp2, y + height, chsigne(A, B, dl0s, dos, z_cl) );
                break;
            case 137:
                tmp1 = -ampli[0] / (ampli[1] - ampli[0]) * width;
                tmp2 = -ampli[1] / (ampli[2] - ampli[1]) * height;
                A.x = x + width;
                A.y = y;
                B.x = x + width / 2.;
                B.y = y + height / 2.;
                appendSegment( x + tmp1, y, x + width, y + tmp2, chsigne(A, B, dl0s, dos, z_cl) );
                tmp1 = -ampli[0] / (ampli[3] - ampli[0]) * height;
                tmp2 = -ampli[3] / (ampli[2] - ampli[3]) * width;
                A.x = x;
                A.y = y + height;
                appendSegment( x, y + tmp1, x + tmp2, y + height, chsigne(A, B, dl0s, dos, z_cl) );
                break;
            case 10:
                tmp1 = -ampli[0] / (ampli[1] - ampli[0]) * width;
                tmp2 = -ampli[3] / (ampli[2] - ampli[3]) * width;
                A.x = x;
                A.y = y;
                B.x = x + width;
                B.y = y;
                appendSegment( x + tmp1, y, x + tmp2, y + height, chsigne(A, B, dl0s, dos, z_cl) );
                break;
            case 11:
                tmp1 = -ampli[0] / (ampli[3] - ampli[0]) * height;
                tmp2 = -ampli[3] / (ampli[2] - ampli[3]) * width;
                A.x = x;
                A.y = y;
                B.x = x;
                B.y = y + height;
                appendSegment( x, y + tmp1, x + tmp2, y + height, chsigne(A, B, dl0s, dos, z_cl) );
                break;
            case 12:
                tmp1 = -ampli[0] / (ampli[3] - ampli[0]) * height;
                tmp2 = -ampli[1] / (ampli[2] - ampli[1]) * height;
                A.x = x;
                A.y = y;
                B.x = x;
                B.y = y + height;
                appendSegment( x, y + tmp1, x + width, y + tmp2, chsigne(A, B, dl0s, dos, z_cl) );
                break;
            case 13:
                tmp1 = -ampli[0] / (ampli[1] - ampli[0]) * width;
                tmp2 = -ampli[1] / (ampli[2] - ampli[1]) * height;
                A.x = x;
                A.y = y;
                B.x = x + width;
                B.y = y;
                appendSegment( x + tmp1, y, x + width, y + tmp2, chsigne(A, B, dl0s, dos, z_cl) );
                break;
            case 14:
                tmp1 = -ampli[0] / (ampli[3] - ampli[0]) * height;
                tmp2 = -ampli[0] / (ampli[1] - ampli[0]) * width;
                A.x = x;
                A.y = y;
                B.x = x + width;
                B.y = y;
                appendSegment( x, y + tmp1, x + tmp2, y, chsigne(A, B, dl0s, dos, z_cl) );
                break;
        }
    }
}

/* Browse the marching square segments and update the flag to 
 * associate them into independent critical lines
 */
static void reflag()
{
    extern struct biline radial[NMAX], tangent[NMAX];
    extern int nrline, ntline, flagr, flagt;

    struct point pi0, pi1;
    double dx, dy;
    int i, j;

    // for tangent critical lines
    flagt = 0;
    for( i = 0; i < ntline; i+=2 )
    {
        tangent[i].i = flagt;
        tangent[i+1].i = flagt;
        pi0 = tangent[i].I;
        pi1 = tangent[i+1].I;

        j = 0;
        dx = tangent[j].I.x - pi1.x;
        dy = tangent[j].I.y - pi1.y;
        while( dx * dx + dy * dy > limitHigh * limitHigh / 4 && j < ntline )
        {
            dx = tangent[j].I.x - pi1.x;
            dy = tangent[j].I.y - pi1.y;
            j += 2;
        }

        tangent[j].i = flagt;
        tangent[j+1].i = flagt;
    }

    // for radial critical lines
    flagr = 0;
    for( i = 0; i < nrline; i+=2 )
    {
        radial[i].i = flagr;
        radial[i+1].i = flagr;
        pi0 = radial[i].I;
        pi1 = radial[i+1].I;

        j = 0;
        dx = radial[j].I.x - pi1.x;
        dy = radial[j].I.y - pi1.y;
        while( dx * dx + dy * dy > limitHigh * limitHigh / 4 && j < nrline )
        {
            dx = radial[j].I.x - pi1.x;
            dy = radial[j].I.y - pi1.y;
            j += 2;
        }

        radial[j].i = flagr;
        radial[j+1].i = flagr;
    }
}

/* Find the caustic/critical lines with a multiresolution
 * marching squares algorithm.
 *
 * The initial square is defined by the CL or F global variables
 *
 * The critical line is used to compute the caustic one in the
 * source plane.
 *
 */
static void marchingSquares(int verbose)
{
    extern struct g_cline  CL;
    extern struct g_frame  F;
    extern struct pot      lens[];
    extern int nrline, ntline, flagr, flagt;

    double xmin, ymin, xmax, ymax;
    double ampli[5];        // amplifications of the root square

    int k;

    /* Definition de la fenetre de calcul */
    if (CL.dmax != 0.)
    {
        xmin = -CL.dmax;
        xmax = CL.dmax;
        ymin = -CL.dmax;
        ymax = CL.dmax;
    }
    else
    {
        xmin = (F.xmax - F.xmin) / 2.;

        ymax = (F.ymin + F.ymax) / 2. + xmin;
        ymin = (F.ymin + F.ymax) / 2. - xmin;

        xmax = (F.xmin + F.xmax) / 2. + xmin;
        xmin = (F.xmin + F.xmax) / 2. - xmin;
    }

    // Limit definitions

    // we don't want square side lengths smaller than CL.cpas
    limitLow = CL.cpas;

    // we want to divide the field at least in 64x64 squares
    if ( CL.limitHigh != 0 )
        limitHigh = CL.limitHigh;
    else
        limitHigh = (xmax - xmin) / 64.;

    /* initialisation de quelques constantes */
    nrline = 0;
    ntline = 0;

#ifdef DEBUG
    system( "rm -f cidbg.reg" );
#endif

    for (k = 0; k < CL.nplan; k++)
    {

        if (verbose)
        {
            fprintf(stderr, "COMP%d: critic and caustic lines for source plane at z=%.3lf\n",
                    k + 1, CL.cz[k]);
            fprintf(stderr, "limitHigh(in arcsec)=%.3lf limitLow(in arcsec)=%.3lf\n",
                    limitHigh, limitLow);

        }
        dl0s = distcosmo2(lens[0].z, CL.cz[k]);
        dos = distcosmo1(CL.cz[k]);
        dlsds_cl = dl0s / dos;
	z_cl = CL.cz[k]; 
        flagr = 1;
        flagt = 1;
        /* prepare the tree*/
        cutSquare( xmin, ymin, xmax - xmin, ymax - ymin, 15, 4, ampli);
    }
}


/* Find the caustic/critical lines with the snake algorithm.
 *
 * It's an iterative method that start from a point and looks
 * always for the next point in a cone oriented along the
 * previous direction.
 *
 * */
static void snake(int verbose)
{
    extern struct g_cline       CL;
    extern struct biline radial[NMAX], tangent[NMAX];
    extern struct pot lens[];
    extern struct g_frame       F;
    extern struct g_grille  G;
    extern int nrline, ntline, flagr, flagt;

    struct point A, B, DPL, N, O, OI, OS, P, Q;
    struct point LP[NLMAX+1], LDPL[NLMAX+1];
    int lnpas[NLMAX+1];
    int nu[NLMAX+1];
    int fsel[NLMAX+1];
    double  ppas, Dmin, dc;
    double  J;
    double  xmin, ymin, xmax, ymax;
    register int    j, k, ii;
    register int i;
    int signe_flag;

    /* Definition de la fenetre de calcul */

    if (CL.dmax != 0.)
    {
        xmin = -CL.dmax;
        xmax = CL.dmax;
        ymin = -CL.dmax;
        ymax = CL.dmax;
    }
    else
    {
        xmin = F.xmin;
        xmax = F.xmax;
        ymin = F.ymin;
        ymax = F.ymax;
    }

    /* initialisation de quelques constantes */

    nrline = 0;
    flagr = 1;
    ntline = 0;
    flagt = 1;
    J = NPOINT;
    ppas = sqrt((xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin)) / NPOINT;

    /* liste des points a suivre */

    if (G.nlens_crit == 1)
    {
        LP[0] = lens[0].C;
        LDPL[0].x = Max(xmax - lens[0].C.x, xmin - lens[0].C.x) / J;
        LDPL[0].y = Max(ymax - lens[0].C.y, ymin - lens[0].C.y) / J;
        if ((lens[0].type == 5) || (lens[0].type == 7) ||
            (lens[0].type == 0) || (lens[0].type == -1) || (lens[0].type == 1))
        {
            LP[0].x += LDPL[0].x;
            LP[0].y += LDPL[0].y;
        };

        lnpas[0] = NPOINT;
    }
    else
    {
        LP[0].x = xmin;
        LP[0].y = ymin;
        Dmin = 9999.;
        for (i = 0; i < G.nlens_crit; i++)
        {
            nu[i] = 99;
            fsel[i] = 0;
        };

        /* on cherche le centre le plus pres du bord en bas a gauche*/
        for (i = 0; i < G.nlens_crit; i++)
        {
            dc = dist(LP[0], lens[i].C);
            if (dc < Dmin)
            {
                Dmin = dc;
                nu[0] = i;
            };
        };
        LP[1] = lens[nu[0]].C;
        fsel[nu[0]] = 1;

        for (i = 1; i < G.nlens_crit; i++)
        {
            Dmin = 9999.;
            for (j = 0; j < G.nlens_crit; j++)
            {
                if (fsel[j] != 1)
                {
                    dc = dist(LP[i], lens[j].C);
                    if (dc < Dmin)
                    {
                        Dmin = dc;
                        nu[i] = j;
                    };
                };
            };
            LP[i+1] = lens[nu[i]].C;
            fsel[nu[i]] = 1;
        };


        for (i = 0; i < G.nlens_crit; i++)
        {
            if ((lens[nu[i]].type == 5) || (lens[nu[i]].type == 7) ||
                (lens[nu[i]].type == 0) || (lens[nu[i]].type == 1))
            {
                LP[i+1].x += 0.01;
                LP[i+1].y += 0.015;
            };

            lnpas[i] = (int) (dist(LP[i], LP[i+1]) / 2. / ppas);
            LDPL[i].x = (LP[i+1].x - LP[i].x) / lnpas[i];
            LDPL[i].y = (LP[i+1].y - LP[i].y) / lnpas[i];
        };
    };

    for (k = 0; k < CL.nplan; k++)
    {
        z_cl = CL.cz[k];
        dl0s = distcosmo2(lens[0].z, z_cl);
        dos = distcosmo1(z_cl);
        dlsds_cl = dl0s / dos;

        if (verbose)
            fprintf(stderr, "COMP%d: critic and caustic lines for source plane at z=%.3lf\n",
                    k + 1, z_cl);

        for (ii = 0; ii < G.nlens_crit; ii++)
        {
            A = LP[ii];
            DPL = LDPL[ii];

            if (verbose)
                fprintf(stderr, "\tnlens=%d npas=%d x0=%.3lf y0=%.3lf\n",
                        ii + 1, lnpas[ii], A.x, A.y);

            for (i = 1; i < lnpas[ii]; A = B, i++)
            {
                B.x = A.x + DPL.x;
                B.y = A.y + DPL.y;
                signe_flag = chsigne(A, B, dl0s, dos, z_cl);

                if (signe_flag == 1)
                {
                    radial[nrline].i = flagr;
                    OI = O = radial[nrline].I = e_zeroamp(A, B, dl0s, dos, z_cl);
                    e_dpl(&O, dlsds_cl, &OS);
                    radial[nrline++].S = OS;

                    N.x = O.x - DPL.y / 2.;
                    N.y = O.y + DPL.x / 2.;
                    radial[nrline].i = flagr;
                    P = radial[nrline].I = next(N, O, CL.cpas / 2., dl0s, dos, z_cl);
                    e_dpl(&P, dlsds_cl, &radial[nrline++].S);

                    radial[nrline].i = flagr;
                    Q = radial[nrline].I = next(O, P, CL.cpas / 2., dl0s, dos, z_cl);
                    e_dpl(&Q, dlsds_cl, &radial[nrline++].S);

                    follow(P, Q, O, dl0s, dos, z_cl);
                    radial[nrline].i = flagr;
                    radial[nrline].I = OI;
                    radial[nrline++].S = OS;
                    flagr++;
                }
                else if (signe_flag == 2)
                {
                    tangent[ntline].i = flagt;
                    OI = O = tangent[ntline].I = e_zeroamp(A, B, dl0s, dos, z_cl);
                    e_dpl(&O, dlsds_cl, &OS);
                    tangent[ntline++].S = OS;
                    N.x = O.x - DPL.y / 2.;
                    N.y = O.y + DPL.x / 2.;
                    tangent[ntline].i = flagt;
                    P = tangent[ntline].I = next(N, O, CL.cpas, dl0s, dos, z_cl);
                    e_dpl(&P, dlsds_cl, &tangent[ntline++].S);
                    tangent[ntline].i = flagt;
                    Q = tangent[ntline].I = next(O, P, CL.cpas, dl0s, dos, z_cl);
                    e_dpl(&Q, dlsds_cl, &tangent[ntline++].S);
                    followi(P, Q, O, dl0s, dos, z_cl);
                    tangent[ntline].i = flagt;
                    tangent[ntline].I = OI;
                    tangent[ntline++].S = OS;
                    flagt++;
                };
            };
        };
    };

}
