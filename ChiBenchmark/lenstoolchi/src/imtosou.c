#include <stdlib.h>
#include <math.h>
#include <constant.h>
#include <fonction.h>
#include <lt.h>
/*
*       nom:        imtosou
*       auteur:     Jean-Paul Kneib
*       creation date:  10/02/92
*       place:      Toulouse
* major modification
*   oct 94
*   IoA Cambridge
*
* aim: transform a collection of points of an image into the source plane
*
*/

void    imtosou(double zimage, char *sname )
{
//  register    int i,j,k;
//  int is,js,size[4];
    const extern  struct g_mode   M;
    extern  struct g_pixel  imFrame, ps;
    const extern  struct  pot lens[];
    double  dlsds;
//  struct  point   A,B;
    double  **im, **source;
    double  **erreur;
    int **imult;

//  char    comment[1024];
    int npl;
    struct  pixlist *pl;    // list of points to send to source plane
    double  ra, dec, width, height;

    /*
    * compute  dlsds
    */

    dlsds = dratio(lens[0].z, zimage);

    /*
    * read image
    */

    im = (double **) readimage(&imFrame);

    /* Allocate the maximum number of points to send to source plane */
    pl = (struct pixlist*)malloc((unsigned int)
                                 imFrame.nx * imFrame.ny * sizeof(struct pixlist));

    /*  keep only the relevant points */
    if (imFrame.ncont == 0)
        s_pixlist(im, pl, &npl);
    else
        keep_cl(im, pl, &npl);

    /*
    * resize the pixel-size in the source plan and copy the wcs information
    */

    ps.pixelx = imFrame.pixelx / ps.ech;
    ps.pixely = imFrame.pixely / ps.ech;


    NPRINTF(stderr, "INFO: Resolution I(%.5lf %.5lf) S(%.5lf %.5lf) (\"/pix) \n",
            imFrame.pixelx, imFrame.pixely, ps.pixelx, ps.pixely);

    /*
    *  allocate square map
    */

    source = alloc_square_double(ps.ny, ps.nx);
    erreur = alloc_square_double(ps.ny, ps.nx);
    imult = alloc_square_int(ps.ny, ps.nx);

    /*
    * define the source box as the barycenter of the source image-center points
    */

    NPRINTF(stderr, "DO: Define source box around c_image barycenter\n");
    s_sourcebox(&ps, M.centerfile, dlsds);

    /*
    *  do the transform image -> source + compute the error
    */

    NPRINTF(stderr, "COMP: image(%d %d) -> source frame (%d %d)\n",
            imFrame.nx, imFrame.ny, ps.nx, ps.ny);

    do_itos(im, pl, npl, dlsds, source, erreur, imult);


    /*
    *  save the results
    */

    // Compute the source image center and relative limits from ps.wcsinfo
    wcsfull(ps.wcsinfo, &ra, &dec, &width, &height);
    ps.xmin -= -3600.*cos(M.ref_dec * DTR) * (ra - M.ref_ra);
    ps.xmax -= -3600.*cos(M.ref_dec * DTR) * (ra - M.ref_ra);
    ps.ymin -= 3600.*(dec - M.ref_dec);
    ps.ymax -= 3600.*(dec - M.ref_dec);


    wrf_fits_abs(sname, source, ps.nx, ps.ny, ps.xmin, ps.xmax, ps.ymin, ps.ymax, ra, dec);
    wri_fits_abs("imult.fits", imult, ps.nx, ps.ny, ps.xmin, ps.xmax, ps.ymin, ps.ymax, ra, dec);
    wrf_fits_abs("erreur.fits", erreur, ps.nx, ps.ny, ps.xmin, ps.xmax, ps.ymin, ps.ymax, ra, dec);

    /*
    *  free the square maps
    */

    free_square_double(source, ps.ny);
    free_square_double(erreur, ps.ny);
    free_square_int(imult, ps.ny);

    free_square_double(im, imFrame.ny);

}
