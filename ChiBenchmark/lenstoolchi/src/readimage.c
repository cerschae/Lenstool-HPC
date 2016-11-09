#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "wcs.h"
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"


/****************************************************************/
/*                 Program  : grille                */
/*                 Version  : 1 mai 1992            */
/*                 Location : Obs. Toulouse         */
/*                 Auteur   : jean-paul             */
/*  EJ 26/12/05 : Modified in format 3 to include WCS
 ****************************************************************
 *
 * Read an image (eventually FITS) and extract the pixels that are inside
 * the contour P.
 *
 * Parameter
 *  - P contains the contour of the image to be inversed
 *
 *  Return an array of long int containing the image[nb_lin=ny][nb_col=nx]
 * of the P->pixfile FITS file.
 */

double **readimage(struct g_pixel *P)
{
    const extern  struct  g_mode          M;
    double  x1, x2, y1, y2;
    double  xmin, xmax, ymin, ymax;

    int i, j;
    register    int ii; //,jj;
    FILE    *IN;

    int     nx, ny;
    double  **zf;
    char    line[50];
    char    type[8];
    char    mode[4];
    char    nature[8];
    char    comments[1024];
    double  zz[10]; //zmin,zmax,
    char    *header;

    IN = fopen(P->pixfile, "r");

    if ( IN != NULL )
    {
        NPRINTF(stderr, "READ: image: %s - format: %d\n", P->pixfile, P->format);

        if ((P->format == 0) || (P->format == -1))
        {
            rewind(IN);
            if (P->header != 0)
                for (ii = 0; ii < P->header; flire(IN, line), ii++);
            if (P->format == 0)
                fscanf(IN, "%d%d", &P->ny, &P->nx);
            else
            {
                fscanf(IN, "%d%d%lf%lf%lf%lf", &P->ny, &P->nx, &x1, &x2, &y1, &y2);
                P->pixelx = (x2 - x1) / ((double)(P->nx - 1));
                P->pixely = (y2 - y1) / ((double)(P->ny - 1));
                P->xmin = x1;
                P->ymin = y1;
                P->xmax = x2;
                P->ymax = y2;
            };
            flire(IN, line);

            for (ii = 0; ii < P->column; fscanf(IN, "%lf", &zz[ii]), ii++);

            rewind(IN);

            if (P->header != 0)
                for (ii = 0; ii < P->header; flire(IN, line), ii++);

            fscanf(IN, "%d%d", &P->ny, &P->nx);
            flire(IN, line);

            zf = (double **) alloc_square_double(P->ny, P->nx);
            i = 0;
            do
            {
                j = 0;
                do
                {
                    for (ii = 0; ii < P->column; fscanf(IN, "%lf", &zz[ii]), ii++);

                    flire(IN, line);
                    zf[i][j] = zz[P->column-1];
                    j++;
                }
                while (j < P->nx);
                i++;
            }
            while (i < P->ny);
        } /*end of if((P->format == 0)||(P->format == -1))*/

        else if ((P->format == 1) || (P->format == 2))
        {
            if (P->format == 1)
            {
                zf = (double **)rdf_ipxs(P->pixfile, &nx, &ny,
                                         type, mode, nature, comments);
                P->ny = nx;
                P->nx = ny;
                NPRINTF(stderr, "\tdim %d %d, %s %s %s\n\t%s\n", nx, ny,
                        type, mode, nature, comments);
            }
            else
            {
                zf = (double **)rdf_ipx(P->pixfile, &ny, &nx,
                                        type, mode, nature, comments, &xmin, &xmax, &ymin, &ymax);
                P->ny = ny;
                P->nx = nx;
                P->pixelx = (xmax - xmin) / ((double)(P->nx - 1));
                P->pixely = (ymax - ymin) / ((double)(P->ny - 1));
                P->xmin = xmin;
                P->ymin = ymin;
                P->xmax = xmax;
                P->ymax = ymax;

                NPRINTF(stderr, "\tdim %d %d, %lf %lf %lf %lf\n\t%s %s %s\n\t%s\n",
                        nx, ny, xmin, xmax, ymin, ymax, type, mode, nature, comments);
            };
        } /*end of if ((P->format ==1)||(P->format ==2))*/
        else if (P->format == 3)
        {
            zf = (double **)rdf_fits(P->pixfile, &nx, &ny, &header);
            P->wcsinfo = wcsinit(header);

            if (P->wcsinfo != NULL && M.iref != 0)
            {
                /*Permute axis if CRTYPE1 is DEC*/

                pix2wcs(P->wcsinfo, 1., 1., &xmin, &ymin);
                pix2wcs(P->wcsinfo, (double)nx, (double)ny, &xmax, &ymax);
            }
            else
            {
                NPRINTF(stderr, "WARN: No astrometrical data in %s\n", P->pixfile);
                P->wcsinfo = NULL;
                xmin = 1;
                xmax = nx;
                ymin = 1;
                ymax = ny;
            }

            NPRINTF(stderr, "INFO: %s absolute bounds in degrees (RA %lf:%lf, DEC %lf:%lf)\n",
                    P->pixfile, xmin, xmax, ymin, ymax);

            P->ny = ny;
            P->nx = nx;
            if (P->pixelx == 0. && P->pixely == 0. && M.iref != 0)
            {
                /*Convert to relative position in arcsec*/
                P->xmin = xmin - M.ref_ra;
                P->xmin *= -3600. * cos(M.ref_dec * DTR);
                P->xmax = xmax - M.ref_ra;
                P->xmax *= -3600. * cos(M.ref_dec * DTR);
                P->ymin = ymin - M.ref_dec;
                P->ymin *= 3600.;
                P->ymax = ymax - M.ref_dec;
                P->ymax *= 3600.;


                P->pixelx = (P->xmax - P->xmin) / ((double)(P->nx - 1));
                P->pixely = (P->ymax - P->ymin) / ((double)(P->ny - 1));
            }
            else
            {
                P->xmax = P->xmin + P->pixelx * (P->nx - 1);
                P->ymax = P->ymin + P->pixely * (P->ny - 1);
            }

            NPRINTF(stderr, "INFO: %s relative bounds in arcsec (RA %lf:%lf, DEC %lf:%lf)\n",
                    P->pixfile, P->xmin, P->xmax, P->ymin, P->ymax);
            NPRINTF(stderr, "INFO: Image size (%d %d) Resolution (%lf %lf) (\"/pix) \n",
                    nx, ny, P->pixelx, P->pixely);
        } /*end of if (P->format==3)*/

        else
        {
            fprintf(stderr, "FATAL ERROR: format %d not known\n", P->format);
            exit(-1);
        };

    } /*end of if( IN != NULL )*/

    else
    {
        fprintf(stderr, "ERROR: image %s not found\n", P->pixfile);
        exit(-1);
    };


    fclose(IN);
    if ((P->pixelx == 0.) || (P->pixely == 0.) || (P->xmin == P->xmax) || (P->ymin == P->ymax))
    {
        fprintf(stderr, "ERROR: image is not scaled !\n");
        exit(-1);
    };

    // Compute sum of pixels
    double fluxtot = 0;
    for( i = 0; i < P->nx; i++ )
        for( j = 0; j < P->ny; j++ )
            fluxtot += zf[j][i];

    P->meanFlux = fluxtot / P->nx / P->ny;
    NPRINTF(stderr, "READ: done\n");
    return(zf);
}
