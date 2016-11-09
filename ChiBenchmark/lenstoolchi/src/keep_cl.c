#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

static int  readlist(struct point I[NPOINT], char file[NIMAX]);

/*
*                 Program  : keep_cl
*                 Version  : oct 94
*                 Location : IoA Cambridge
*                 Auteur   : jean-paul
*
* Extract the points of image that are inside the polygonal contour pl
* and save the result in the imFrame.outfile FITS file in double format.
* This FITS file is 0 ouside of the imFrame.ncont polygonal contours.
*
* The polygonal contours are read in the files imFrame.contfile[k].
* k varying from 0 to imFrame.ncont. Each contour can contain one of
* the image of the multiple images system.
*
* imFrame is a global variable.
*
* Parameters :
*  - image is the original full image of the field
*
* Return :
*  - pl.{i=row,j=col} is the list of the pixels inside the contour
*  - npl is the number of pixels inside the polygonal contour
*  - imFrame.outfile in FITS file format is written on disk
*
*/

void    keep_cl(double **image, struct pixlist *pl, int *npl)
{
    extern  struct  g_mode      M;
    extern struct   g_pixel     imFrame;

    register    int k, kk, ii, jj;
    double  **clean;

//  int *size;
    struct  point   P, I[NPOINT];
    int np;

    double xpixel, ypixel, xw, yw;
    int cxmin, cxmax, cymin, cymax, offscl;

    /*clean is in [nb_row][nb_col] format*/
    clean = (double **) alloc_square_double(imFrame.ny, imFrame.nx);

    NPRINTF(stderr, "COMP: clean image (%d,%d)[%.5lf,%.5lf]\n",
            imFrame.nx, imFrame.ny, imFrame.pixelx, imFrame.pixely);

    kk = 0;
    for (k = 0; k < imFrame.ncont; k++)
    {
        np = readlist(I, imFrame.contfile[k]);
        NPRINTF(stderr, "KEEP: contour %d #%s# %d points \n", k + 1, imFrame.contfile[k], np);
        /*Limit the search region for the contour in the image*/
        if (imFrame.wcsinfo != NULL)
        {
            cxmin = imFrame.nx;
            cxmax = 0;
            cymin = imFrame.ny;
            cymax = 0;

            for (ii = 0; ii < np; ii++)
            {
                xw = (double)I[ii].x;
                yw = (double)I[ii].y;
		// points in contfile[k] are in relative arcsec
                xw = xw / (-3600.) / cos(M.ref_dec * DTR) + M.ref_ra;
                yw = yw / 3600. + M.ref_dec;
                wcs2pix(imFrame.wcsinfo, xw, yw, &xpixel, &ypixel, &offscl);
                if (offscl == 0)
                {
                    offscl = (int)xpixel;
                    if (cxmin > offscl)     cxmin = offscl;
                    if (cxmax < offscl)     cxmax = offscl;
                    offscl = (int)ypixel;
                    if (cymin > offscl)     cymin = offscl;
                    if (cymax < offscl)     cymax = offscl;
                }
            } /*end of for each point of the contour*/
        }
        else /*imFrame.wcsinfo is NULL*/
        {
            /*TODO*/
            cxmin = 0;
            cymin = 0;
            cxmax = imFrame.nx;
            cymax = imFrame.ny;
        }

        /*Look for the pixels inside the contour*/
        for (ii = cymin; ii < cymax; ii++)
        {
            P.y = imFrame.ymin + imFrame.pixely * ii;
            for (jj = cxmin; jj < cxmax; jj++)
            {
                if (imFrame.wcsinfo != NULL)
                {
                    xpixel = (double)(jj + 1);
                    ypixel = (double)(ii + 1);
                    pix2wcs(imFrame.wcsinfo, xpixel, ypixel, &xw, &yw);
                    P.x = xw - M.ref_ra;
                    P.x *= -3600.*cos(M.ref_dec * DTR);
                    P.y = yw - M.ref_dec;
                    P.y *= 3600.;
                }
                else
                    /*increment along x image axis*/
                    P.x = imFrame.xmin + imFrame.pixelx * jj;

                if (inconvexe(P, np, I) > 0)
                {
                    pl[kk].i = P.y;
                    pl[kk].j = P.x;
                    pl[kk++].flux = image[ii][jj];
                    clean[ii][jj] = image[ii][jj];
                }

            }
        }
    }
    NPRINTF(stderr, "KEEP: %d pixel\n", kk);
    *npl = kk;

//    if( imFrame.wcsinfo == NULL )
        wrf_fits(imFrame.outfile, clean, imFrame.nx, imFrame.ny, imFrame.xmin, imFrame.xmax, imFrame.ymin, imFrame.ymax);
//    else
//    {
//        double ra, dec, width, height;
//        wcsfull(imFrame.wcsinfo, &ra, &dec, &width, &height);
//        wrf_fits_abs(imFrame.outfile, clean, imFrame.nx, imFrame.ny, 
//              imFrame.xmin, imFrame.xmax, imFrame.ymin, imFrame.ymax, ra, dec);
//    }
}




/*
*                 Program  : s_pixlist
*                 Version  : oct 94
*                 Location : IoA Cambridge
*                 Auteur   : jean-paul
*/

void    s_pixlist(double **image, struct pixlist *pl, int *npl)
{
    extern  struct  g_mode          M;
    extern struct   g_pixel imFrame;
    register    int kk, ii, jj;

    kk = 0;
    for (ii = 0; ii < imFrame.ny; ii++)
        for (jj = 0; jj < imFrame.nx; jj++)
            if (image[ii][jj] != 0.)
            {
                pl[kk].i = ii;
                pl[kk++].j = jj;
            }

    NPRINTF(stderr, "KEEP: %d pixel\n", kk);
    *npl = kk;
}

/****************************************************************
 * In the cleanlens mode, read the contour file
 * Return
 *  - I the array of point filled with the contour
 *  - k the number of elements in I
 */

static int  readlist(struct point I[NPOINT], char file[NIMAX])
{
    extern  struct  g_mode          M;
    int n, k = 0;
    int stop = 0;
    FILE    *IN;

    IN = fopen(file, "r");
    if ( IN != NULL )
        while (stop != -1)
        {
            stop = fscanf(IN, "%d%lf%lf\n", &n, &I[k].x, &I[k].y);

            if (stop == 3)
                k++;
            else if (stop > 0)
            {
                fprintf(stderr, "ERROR: %s is not in the right format\n", file);
                exit(-1);
            }
        }
    else
    {
        fprintf( stderr, "ERROR: Opening the contour file %s\n", file);
        exit(-1);
    }
    fclose(IN);
    return(k);
}

