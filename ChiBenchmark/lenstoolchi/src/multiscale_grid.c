#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<float.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<fonction.h>

/****************************************************************/
/*      nom:        multiscale_grid             */
/*      auteur:     Eric Jullo          */
/*      date:       22/02/08            */
/*      place:      Marseille           */
/****************************************************************
 * Build a multiscale grid according to the density of constraints
 * read from the multfile and the field defined by champ section.
 *
 * Increment the number of clumps in the lens[] variable.
 *
 * Exit if multfile is not defined.
 * The Champ section default is defined in set_default()
 *
 * If there are no constraints in the field, split the field only once, ie. 6 triangles.
 */

#define COPY(p1,p2) {p1.x=p2.x; p1.y=p2.y;}

static struct galaxie mult[NFMAX*NIMAX]; // LOCAL catalog of multiple images
static long int   nmult;    // number of elements in mult[]

static long int   ilens;    // index of the last clump of the grid in lens[]
static struct pot **gridID;  // grid index containing pointers to registered clumps
static double     dx, dy;   // size of the smallest triangle in the grid at the finest resolution
static double     width;    // width of the grid
static unsigned int maxLens;// maximum number of lens in the grid for a given number of splitting

static void divide_tri(struct triplet *tri, double rci, int leveli);
static double dens_tri(struct triplet *tri);
static void readMultfile();
static void buildgrid(long int *pilens);
static void circle();
static void regLens(struct point *ppoint, double rc);
static unsigned int getID(struct point *ppoint);
static void initClump(struct point *ppoint, double rc, double lmin, double lmax );
static void readGridfile(long int *pilens);

void multiscale_grid(long int *pilens)
{
    // Define variables
    extern struct g_msgrid H;
    extern struct g_image  I;

    if ( strcmp(H.gridfile, ""))
    {
        readGridfile(pilens);
        return;
    }

    // Check that the multfile is defined
    if ( I.n_mult == 0 )
    {
        fprintf(stderr, "ERROR[multiscale grid] multfile undefined in .par file.");
        exit(-1);
    }

    // Read the catalog of multiple images
    readMultfile();

    // Place the PIEMD clumps and define their rcore/rcut
    buildgrid(pilens);

#ifdef DEBUG
    // Create a DS9 file with the clump positions
    circle();
#endif
}

/* Read a grid file and initialise it
 */
static void readGridfile(long int *pilens)
{
    extern struct g_mode   M;
    extern struct g_grille G;
    extern struct g_msgrid H;
    extern struct pot      lens[], lmin[], lmax[];
    extern int    block[][NPAMAX];

    char   line[128];
    double ra, dec;
    int    iref;
    FILE   *in;

    in = fopen(H.gridfile, "r");
    if ( in == NULL )
    {
        fprintf( stderr, "ERROR: reading %s\n", H.gridfile);
        exit(-1);
    }

    ilens = G.nmsgrid;
    while ( !feof(in) && !ferror(in) && ilens < NLMAX )
    {
        fgets(line, 128, in);
        if ( strstr(line, "#REFERENCE" ) != NULL )
        {
            getRADEC(line, &iref, &ra, &dec );
        }
        else if ( sscanf( line, "%d%lf%lf%lf%lf%lf%lf%lf",
                          &lens[ilens].type, &lens[ilens].C.x, &lens[ilens].C.y,
                          &lens[ilens].rc, &lens[ilens].rcut, &lens[ilens].sigma,
                          &lmin[ilens].sigma, &lmax[ilens].sigma) == 8 )
        {
            // Convert input from relative to absolute coordinates if needed
            convertXY( &lens[ilens].C.x, &lens[ilens].C.y, iref, ra, dec);

            // Give an ID to each clump
            sprintf(lens[ilens].n, "%ld", ilens);

            // Initialize the lens parameters
            lens[ilens].epot = 0.;
            lens[ilens].theta = 0;
            lens[ilens].alpha = 0;

            // convert to output relative coordinates
            if ( M.iref == 1 || M.iref == 3 )
            {
                lens[ilens].C.x -= M.ref_ra;
                lens[ilens].C.x *= -3600 * cos(M.ref_dec * DTR);
                lens[ilens].C.y -= M.ref_dec;
                lens[ilens].C.y *= 3600;
            }
            else if ( M.iref == 2 )
            {
                lens[ilens].C.x -= M.ref_ra;
                lens[ilens].C.y -= M.ref_dec;
            }

            // set uniform prior
            block[ilens][B0] = 1;

            ilens++;
        }
    }

    *pilens = ilens;
}

/* Create a catalog of the clump positions to be plotted with pelli
 * DEBUG USE ONLY
 */
static void circle()
{
    extern struct g_mode   M;
    extern struct g_grille G;
    extern struct pot      lens[];
    struct pot  *plens;
    FILE   *ds9;
    int    i;

    ds9 = fopen("msgrid.dat", "w");
    fprintf( ds9, "#REFERENCE  3  %lf  %lf\n ", M.ref_ra, M.ref_dec);
    for ( i = G.nmsgrid; i < ilens; i++ )
    {
        plens = &lens[i];
        fprintf( ds9, "%s  %lf %lf %lf %lf 0. 0. 0.\n",
                 plens->n, plens->C.x, plens->C.y, plens->rc, plens->rc);
    }
    fclose(ds9);
}

/* Place the PIEMD clumps and define their rcore/rcut
 * Parameters:
 * - pilens : global number of lenses in lens[] list
 */
static void buildgrid(long int *pilens)
{
    // Define variables
    extern struct g_frame  F;
    extern struct g_msgrid H;

    int    id[] = {0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 5, 6, 0, 6, 1};
    struct triplet tri;
    struct point   node[7];// summits of the hexagon and center
    double rc;             // parameter of the PIEMD clumps
    double s3 = sqrt(3.) / 2.; // const
    unsigned int    i;

    // find the limits of the field in arcsec
    width = F.xmax - F.xmin;

    // Initialise the grid parameter
    dx = width / 2 / (1 << (H.levels - 1));
    dy = s3 * dx;

    // Max number of nodes in the grid
    maxLens = 1 << (H.levels - 1);  // 2^(H.levels-1)
    maxLens = 1 + 6 * maxLens * ( maxLens + 1 ) / 2;  // arithmetic series
    gridID = (struct pot **) malloc( maxLens * sizeof( struct pot *) );
    for ( i = 0; i < maxLens ; i++ )
        gridID[i] = 0;

    rc = width / 2.;  // because we are sure that we go at least to level1

    // 7 initial points in relative arcsec
    ilens = *pilens;

    node[0].x = 0.;
    node[0].y = 0.;
    node[1].x = width / 2.;
    node[1].y = 0.;
    node[2].x = width / 4.;
    node[2].y = s3 * width / 2.;
    node[3].x = -width / 4.;
    node[3].y = s3 * width / 2.;
    node[4].x = -width / 2.;
    node[4].y = 0.;
    node[5].x = -width / 4.;
    node[5].y = -s3 * width / 2.;
    node[6].x = width / 4.;
    node[6].y = -s3 * width / 2.;

    // Divide subtriangles 0 -> 5
    if ( H.levels > 1 )
        for ( i = 0; i < 6; i++ )
        {
            COPY(tri.a, node[id[3*i]]);
            COPY(tri.b, node[id[3*i+1]]);
            COPY(tri.c, node[id[3*i+2]]);
            if ( dens_tri(&tri) >= H.threshold )
                divide_tri(&tri, rc, 1);
        }

    // free gridID
    free(gridID);

    *pilens = ilens;
}

/* Map a grid position to a location in the gridID[] array.
 * Return the location in the gridID[] array
 */
static unsigned int getID(struct point *ppoint)
{
    int i, j, n;

    j = ppoint->y / dy;
    i = ( ppoint->x + j * dx / 2. ) / dx;
    n = (int) (abs(j) * width / dx + abs(j) - abs(j) * (abs(j) + 1) / 2);
    n = j > 0 ? n : -n;
    return( (unsigned int) (n + i + (maxLens - 1) / 2) );
}

/* Initialise a lens[] clumps
 */
static void initClump(struct point *ppoint, double rc, double min, double max )
{
    extern struct g_msgrid H;
    extern struct pot      lens[];
    extern struct pot      lmin[], lmax[];
    extern int    block[][NPAMAX];

    struct pot   *plens;

    // Initialise the new clump
    plens = &lens[ilens];

    sprintf(plens->n, "G%ld", ilens);
    plens->type = 81;    // PIEMD
    plens->C = *ppoint;
    plens->emass = 0.;
    plens->theta = 0.;
    plens->rc = rc;
    plens->rcut = rc * H.param;
    plens->rckpc = 0;       // for compatibility with set_lens()
    plens->rcutkpc = DBL_MAX;
    plens->mag = 0;
    plens->sigma = 100.;
    plens->z = lens[0].z;


    // Priors
    block[ilens][B0] = 1;
    lmin[ilens].sigma = min;
    lmax[ilens].sigma = max;

}

/* Register a clumps in the gridID[] array
 */
static void regLens(struct point *ppoint, double rc)
{
    extern struct pot  lens[];

    // Initialise the new clump
    initClump(ppoint, rc, 0., 1200.);

    // Register lens in gridID[] array
    gridID[getID(ppoint)] = &lens[ilens];

    ilens++;
}

/* Divide a triangle into 4 subtriangles
 * Append the new subtriangles positions to the x and y list
 * Increment the nlens global variable
 *
 *            2
 *          /   \
 *         / [2] \
 *        /       \
 *     5 ---------- 4
 *    /   \       /   \
 *   / [0] \ [3] / [1] \
 *  /       \   /       \
 * 0 -------  3 -------- 1
 *
 *  - tri : corners position (0,1,2) in arcsec of this triangle
 *  - rci : scale radius of the clumps in this triangle
 *  - leveli : level of this triangle
 */
static void divide_tri( struct triplet *tri, double rci, int leveli)
{
    extern struct g_msgrid H;

    struct point   node[6];   // triangle and subtriangle nodes
    struct pot     *plens;  // pointer to a lens[] element
    struct triplet trio;
    double rcs;             // subtriangle PIEMD params
    double dens;            // density of constraints in subtriangle
    int    nsplit = 0;      // boolean to say if this triangle is splitted or not
    int    i;

    // this triangle corners
    COPY(node[0], tri->a);
    COPY(node[1], tri->b);
    COPY(node[2], tri->c);

    if ( leveli + 1 <= H.levels )
    {
        // subtriangle corners and size
        node[3].x = (node[0].x + node[1].x) / 2.;
        node[4].x = (node[1].x + node[2].x) / 2.;
        node[5].x = (node[0].x + node[2].x) / 2.;
        node[3].y = (node[0].y + node[1].y) / 2.;
        node[4].y = (node[1].y + node[2].y) / 2.;
        node[5].y = (node[0].y + node[2].y) / 2.;
        rcs = rci / 2.;

        // subtriangle 0
        COPY(trio.a, node[0]);
        COPY(trio.b, node[3]);
        COPY(trio.c, node[5]);
        dens = dens_tri(&trio);
        //printf("%lf  %lf\n",dens,H.threshold);
        if ( dens >= H.threshold )
        {
            nsplit = 1;
            divide_tri(&trio, rcs, leveli + 1);
        }
        // subtriangle 1
        COPY(trio.a, node[3]);
        COPY(trio.b, node[1]);
        COPY(trio.c, node[4]);
        dens = dens_tri(&trio);
        //printf("%lf  %lf\n",dens,H.threshold);
        if ( dens >= H.threshold )
        {
            nsplit = 1;
            divide_tri(&trio, rcs, leveli + 1);
        }
        // subtriangle 2
        COPY(trio.a, node[5]);
        COPY(trio.b, node[4]);
        COPY(trio.c, node[2]);
        dens = dens_tri(&trio);
        //printf("%lf  %lf\n",dens,H.threshold);
        if ( dens >= H.threshold )
        {
            nsplit = 1;
            divide_tri(&trio, rcs, leveli + 1);
        }

        // subtriangle 3 (do not add its point 5)
        COPY(trio.a, node[3]);
        COPY(trio.b, node[4]);
        COPY(trio.c, node[5]);
        dens = dens_tri(&trio);
        //printf("%lf  %lf\n",dens,H.threshold);
        if ( dens >= H.threshold )
        {
            nsplit = 1;
            divide_tri(&trio, rcs, leveli + 1);
        }
    }

    if ( nsplit == 1 ) rci = rcs;
    // if not splitting triangle or leave triangle --> register
    //if( leveli+1 > H.levels || nsplit == 0 )
    {
        // register subtriangle corners or rescale them
        for ( i = 0; i <= 2 ;  i++ )
        {
            if ( (plens = gridID[getID(&node[i])]) == NULL )
                regLens(&node[i], rci);
            else if ( plens->rc > rci )
                plens->rc = rci;
        }
    }
}

/* Return the number of constraints inside a triangle
 */
static double dens_tri(struct triplet *tri)
{
    long int i;
    int sum = 0;

    for ( i = 0; i < nmult; i++ )
        sum += inside(&mult[i].C, tri);

    return sum;
}

/* Read the multfile catalog and save the images in mult variable
 */
static void readMultfile()
{
    extern struct g_image I;
    f_shape(&nmult, mult, I.multfile, 1);
}
