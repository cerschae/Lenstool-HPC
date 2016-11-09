#include <fonction.h>


/*
*       nom:        o_invim
*       auteur:     Jean-Paul Kneib
*       date:       10/02/92
*       place:      Toulouse
*
* major changes:
*       oct 94
*       IoA cambridge
*
* preparation of the image-inversion
* Parameters :
* - zim(in) is the redshift of the image to inverse
* - drim(out) is the distance ratio between the image and the main lens
* - plo(out) is the polygon around the image to inverse
* - nplot(out) is the number of points of the polygon around the image
*/

double  **o_invim(  double zim,
                    double *drim,
                    struct pixlist *plo,
                    int *nplo )
{
    extern  struct g_pixel  imFrame, ps;

    const extern  struct  g_mode          M;
    const extern  struct  pot lens[];
    double  **im;
//  struct  point   A,B;

    /* calcul du dlsds */
    *drim = dratio(lens[0].z, zim);

    /*
    * read image
    */

    NPRINTF(stderr, "READ: image frame\n");

    im = (double **) readimage(&imFrame);

    /*
    *  keep only the relevant points
    */

    if (imFrame.ncont == 0)
        s_pixlist(im, plo, nplo);
    else
        keep_cl(im, plo, nplo);

    /*
    * resize the pixel-size in the source plan
    */

    ps.pixelx = imFrame.pixelx / ps.ech;
    ps.pixely = imFrame.pixely / ps.ech;

    NPRINTF(stderr, "SIZE: pixel I[%.5lf,%.5lf] S[%.5lf,%.5lf]\n",
            imFrame.pixelx, imFrame.pixely, ps.pixelx, ps.pixely);

    return((double **) im);
}
