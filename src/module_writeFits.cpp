/**
* @file   module_writeFits.c
* @Author Markus Rexroth, EPFL
* @date   July 2015
* @version 0,1
* @brief  Write fits files
*
* Markus Rexroth, EPFL (changes to function names, variables names, a few function calls etc. to comply with
		modular construction of new LENSTOOL). Main function coding etc. was done by:

		Author: JPK
		Date: April 2001
		Write a double array:
		wrf_fits(filename,ima,nx,ny,xmin,xmax,ymin,ymax)
		wrf_fits_abs(filename,ima,nx,ny,xmin,xmax,ymin,ymax,ra,dec)
		Write a int array:
		wri_fits(filename,ima,nx,ny,xmin,xmax,ymin,ymax)
		wri_fits_abs(filename,ima,nx,ny,xmin,xmax,ymin,ymax,ra,dec)

		Update: Johan Richard
		Date: 15 May 2013
		Write a double 3D array in a datacube:
		wrf_cube_fits(filename,cube,nx,ny,nz,xmin,xmax,ymin,ymax,lmin,lmax)
		wrf_cube_fits_abs(filename,cube,nx,ny,nz,xmin,xmax,ymin,ymax,lmin,lmax,ra,dec)
*/







/// Include
///==========================================================================================================
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "fitsio.h"
#include "module_writeFits.hpp"







/// Functions
///==========================================================================================================

void module_writeFits(std::string path, std::string filename, int ii, type_t *map, const struct runmode_param* runmode, int size, const struct grid_param* frame, type_t ra, type_t dec ){

	std::string file;
	file = path;
	file.append("/");
	file.append(filename);
	file.append("_");
    std::ostringstream ss;
    ss << ii;
	file.append(ss.str());
	file.append(".fits");
	ss.clear();
	char file_char[file.length()+1];
	strcpy(file_char,file.c_str());

	if(ra == 0. and dec == 0.)
		module_writeFits_Image(file_char,map,size,size,frame->xmin,frame->xmax,frame->ymin,frame->ymax);
	else
		module_writeFits_ImageAbsoluteCoordinates(file_char,map,size,size,frame->xmin,frame->xmax,frame->ymin,frame->ymax,ra,dec);
}


/** @brief CPU Print error message
 *
* @param status	status
*/
// Print error message
void printerror( int status)
{

    if (status)
    {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
}



/** @brief CPU Create a FITS primary array containing a 2-D image
 * Write in double array
 *
* @param filename	name of fits file
* @param ima		images 2D array
* @param nx			x dimension
* @param ny			y dimension
* @param xmin		x min value
* @param xmax		x max value
* @param ymin		y min value
* @param ymax		y max value
*/
int module_writeFits_Image(char *filename, double *ima, int nx,int ny, double xmin,double xmax,double ymin,double ymax)
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, ii, jj, k;
    long  fpixel, nelements;
    float crval1,crpix1,cdelt1;
    float crval2,crpix2,cdelt2;
    float bzero=0.,bscale=1.;
    char  cunit[10];
    char  ctype[10];
    char  origin[10];
    float *p_ima;
    /* initialize FITS image parameters */
    int bitpix=FLOAT_IMG;
    long naxis=2;                 /* 2-dimensional image                */
    long naxes[2] = { nx, ny };   /* image is nx pixels wide by ny rows */

	printf("WRITE: FITS file: %s[%d,%d] (%.2lf:%.2lf,%.2lf:%.2lf)\n",
		filename,nx,ny,xmin,xmax,ymin,ymax);

	/* Transform the 2D array to a 1D array */
    p_ima=(float *)malloc((unsigned) (nx*ny)*sizeof(float));
    k=0;
    for (jj=0;jj<ny;jj++)
    	for (ii=0;ii<nx;ii++)
        {
        	p_ima[k]=(float)ima[jj*nx+ii];
        	k++;
        }


    remove(filename);               /* Delete old file if it already exists */

    status = CREATE_DISK_FILE;         /* initialize status before calling fitsio routines */
    if (fits_create_file(&fptr, filename, &status)) /* create new FITS file */
         printerror( status );           /* call printerror if error occurs */

    if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
         printerror( status );


    fpixel = 1;                               /* first pixel to write      */
    nelements = naxes[0] * naxes[1];          /* number of pixels to write */


    crval1=xmin;
    crpix1=1;
    cdelt1=(xmax-xmin)/(nx-1);
    crval2=ymin;
    crpix2=1;
    cdelt2=(ymax-ymin)/(ny-1);

	strcpy(cunit,"arcsec");
	strcpy(ctype,"LINEAR");
	strcpy(origin,"JPK-SOFT");

    /* write the array to the FITS file */
    if ( fits_write_img(fptr, TFLOAT, fpixel, nelements, p_ima, &status) )
        printerror( status );

    /* write another optional keyword to the header */
    if ( fits_update_key(fptr, TFLOAT, "BZERO", &bzero, "BZERO", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "BSCALE", &bscale, "BSCALE", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "CRVAL1", &crval1, "CRVAL1", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "CRPIX1", &crpix1, "CRPIX1", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "CDELT1", &cdelt1, "CDELT1", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TSTRING, "CUNIT1", &cunit, "CUNIT1", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TSTRING, "CTYPE1", &ctype, "CTYPE1", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "CRVAL2", &crval2, "CRVAL2", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "CRPIX2", &crpix2, "CRPIX2", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "CDELT2", &cdelt2, "CDELT2", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TSTRING, "CUNIT2", &cunit, "CUNIT2", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TSTRING, "CTYPE2", &ctype, "CTYPE2", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TSTRING, "ORIGIN", &origin, "ORIGIN", &status) )
         printerror( status );
    if ( fits_write_date(fptr, &status) )
         printerror( status );

    if ( fits_close_file(fptr, &status) )                /* close the file */
         printerror( status );

    return 0;
}





/** @brief CCreate a FITS primary array containing a 2-D image
 * Write in double array and in absolute coordinates
 * EJ 27/12/05 Modified Keyword dectype="DEC---TAN"
 * NOT IMPLEMENTED YET
 *
* @param filename	name of fits file
* @param ima		images 2D array
* @param nx			x dimension
* @param ny			y dimension
* @param xmin		x min value
* @param xmax		x max value
* @param ymin		y min value
* @param ymax		y max value
* @param ra			ra
* @param dec		dec
*/

int module_writeFits_ImageAbsoluteCoordinates(char *filename, double *ima, int nx,int ny, double xmin,double xmax,double ymin,double ymax, double ra,double dec)
{
	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	int status, ii, jj, k;
	long  fpixel, nelements;
	double crval1,crpix1,cdelt1;
	double crval2,crpix2,cdelt2;
	double bzero=0.,bscale=1.,equinox;
	char  cunit[10];
	char  system[10];
	char  ratype[10];
	char  dectype[10];
	char  origin[10];
	float *p_ima;

	/* initialize FITS image parameters */
	int bitpix=FLOAT_IMG;
	long naxis=2;             /* 2-dimensional image                */
	long naxes[2] = { nx, ny };   /* image is nx pixels wide by ny rows */

	printf("WRITE: FITS file: %s[%d,%d] (%.2lf:%.2lf %.2lf:%.2lf)\n",
		filename,nx,ny,xmin,xmax,ymin,ymax);

	/* Transform the 2D array to a 1D array */

	p_ima=(float *)malloc((unsigned) (nx*ny)*sizeof(float));
	k=0;
	for (jj=0;jj<ny;jj++)
		for (ii=0;ii<nx;ii++)
	{
		p_ima[k] =(float)ima[jj*nx+ii];
		k++;
	}


	remove(filename);   /* Delete old file if it already exists */

	status = 0; /* initialize status before calling fitsio routines */

	if (fits_create_file(&fptr, filename, &status)) /* create new FITS file */
		printerror( status );   /* call printerror if error occurs */

	if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
		printerror( status );


	fpixel = 1;   /* first pixel to write  */
	nelements = naxes[0] * naxes[1];  /* number of pixels to write */


	crval1=ra;
	crpix1=-(nx - 1)*xmin/(xmax-xmin) + 1;  // nx-1 because if xmin=0, xmax=3, nx=4
	cdelt1=-(xmax-xmin)/(nx-1)/3600;        // then delta=(xmax - xmin) / (nx-1)
	crval2=dec;
	crpix2=-(ny - 1)*ymin/(ymax-ymin) + 1;
	cdelt2=(ymax-ymin)/(ny-1)/3600;

	strcpy(cunit,"deg");
	strcpy(origin,"JPK-SOFT");
	equinox=2000.0;
	strcpy(cunit,"deg");
	strcpy(origin,"JPK-SOFT");
	strcpy(system,"FK5");
	strcpy(ratype,"RA---TAN");
	strcpy(dectype,"DEC--TAN");


	/* write the array to the FITS file */
	if ( fits_write_img(fptr, TFLOAT, fpixel, nelements, p_ima, &status) )
		printerror( status );

	/* write another optional keyword to the header */
	if ( fits_update_key(fptr, TDOUBLE, "BZERO", &bzero, "BZERO", &status) )
		printerror( status );
	if ( fits_update_key(fptr, TDOUBLE, "BSCALE", &bscale, "BSCALE", &status) )
		printerror( status );
	if ( fits_update_key(fptr, TDOUBLE, "EQUINOX", &equinox, "EQUINOX", &status) )
		printerror( status );
	if ( fits_update_key(fptr, TSTRING, "RADECSYS", &system, "RADECSYS", &status) )
		printerror( status );
	if ( fits_update_key(fptr, TSTRING, "CTYPE1", &ratype, "CTYPE1", &status) )
		printerror( status );
	if ( fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1, "CRVAL1", &status) )
		printerror( status );
	if ( fits_update_key(fptr, TDOUBLE, "CRPIX1", &crpix1, "CRPIX1", &status) )
		printerror( status );
	if ( fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1, "CDELT1", &status) )
		printerror( status );
	if ( fits_update_key(fptr, TSTRING, "CUNIT1", &cunit, "CUNIT1", &status) )
		printerror( status );
	if ( fits_update_key(fptr, TSTRING, "CTYPE2", &dectype, "CTYPE2", &status) )
		printerror( status );
	if ( fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2, "CRVAL2", &status) )
		printerror( status );
	if ( fits_update_key(fptr, TDOUBLE, "CRPIX2", &crpix2, "CRPIX2", &status) )
		printerror( status );
	if ( fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2, "CDELT2", &status) )
		printerror( status );
	if ( fits_update_key(fptr, TSTRING, "CUNIT2", &cunit, "CUNIT2", &status) )
		printerror( status );
	if ( fits_update_key(fptr, TSTRING, "ORIGIN", &origin, "ORIGIN", &status) )
		printerror( status );
	if ( fits_write_date(fptr, &status) )
		printerror( status );

	if ( fits_close_file(fptr, &status) )/* close the file */
		printerror( status );

	return 0;
}





/** @brief Write fits cube
 *
* @param filename	name of fits file
* @param cube		3D array
* @param nx			x dimension
* @param ny			y dimension
* @param nz			z dimension
* @param xmin		x min value
* @param xmax		x max value
* @param ymin		y min value
* @param ymax		y max value
* @param lmin		l min value
* @param lmax		l max value
*/

int module_writeFits_cube(char *filename, double ***cube, int nx,int ny, int nz, double xmin,double xmax,double ymin,double ymax, double lmin, double lmax)
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, ii, jj, kk, k;
    long  fpixel, nelements;
    float crval1,crpix1,cdelt1;
    float crval2,crpix2,cdelt2;
    float crval3,crpix3,cdelt3;
    float bzero=0.,bscale=1.;
    char  cunit[10];
    char  cunitw[10];
    char  ctype[10];
    char  origin[10];
    float *p_ima;
    /* initialize FITS image parameters */
    int bitpix=FLOAT_IMG;
    long naxis=3;                 /* 3-dimensional image                */
    long naxes[3] = { nx, ny, nz };   /* image is nx pixels wide by ny rows by nz slices*/

	printf("WRITE: FITS file: %s[%d,%d,%d] (%.2lf:%.2lf,%.2lf:%.2lf,%.2lf:%.2lf)\n",
		filename,nx,ny,nz,xmin,xmax,ymin,ymax,lmin,lmax);

	/* Transform the 3D array to a 1D array */
    p_ima=(float *)malloc((unsigned) (nx*ny*nz)*sizeof(float));
    k=0;
    for (kk=0;kk<nz;kk++)
        for (jj=0;jj<ny;jj++)
    	    for (ii=0;ii<nx;ii++)
            {
        	p_ima[k]=(float)cube[jj][ii][kk];
        	k++;
            }


    remove(filename);               /* Delete old file if it already exists */

    status = CREATE_DISK_FILE;         /* initialize status before calling fitsio routines */
    if (fits_create_file(&fptr, filename, &status)) /* create new FITS file */
         printerror( status );           /* call printerror if error occurs */

    if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
         printerror( status );


    fpixel = 1;                               /* first pixel to write      */
    nelements = naxes[0] * naxes[1] * naxes[2];          /* number of pixels to write */


    crval1=xmin;
    crpix1=1;
    cdelt1=(xmax-xmin)/(nx-1);
    crval2=ymin;
    crpix2=1;
    cdelt2=(ymax-ymin)/(ny-1);
    crval3=lmin;
    crpix3=1;
    cdelt3=(lmax-lmin)/(nz-1);

	strcpy(cunit,"arcsec");
	strcpy(cunitw,"Angstroms");
	strcpy(ctype,"LINEAR");
	strcpy(origin,"JPK-SOFT");

    /* write the array to the FITS file */
    if ( fits_write_img(fptr, TFLOAT, fpixel, nelements, p_ima, &status) )
        printerror( status );

    /* write another optional keyword to the header */
    if ( fits_update_key(fptr, TFLOAT, "BZERO", &bzero, "BZERO", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "BSCALE", &bscale, "BSCALE", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "CRVAL1", &crval1, "CRVAL1", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "CRPIX1", &crpix1, "CRPIX1", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "CDELT1", &cdelt1, "CDELT1", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TSTRING, "CUNIT1", &cunit, "CUNIT1", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TSTRING, "CTYPE1", &ctype, "CTYPE1", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "CRVAL2", &crval2, "CRVAL2", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "CRPIX2", &crpix2, "CRPIX2", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "CDELT2", &cdelt2, "CDELT2", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TSTRING, "CUNIT2", &cunit, "CUNIT2", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TSTRING, "CTYPE2", &ctype, "CTYPE2", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "CRVAL3", &crval3, "CRVAL3", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "CRPIX3", &crpix3, "CRPIX3", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TFLOAT, "CDELT3", &cdelt3, "CDELT3", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TSTRING, "CUNIT3", &cunitw, "CUNIT3", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TSTRING, "CTYPE3", &ctype, "CTYPE3", &status) )
         printerror( status );
    if ( fits_update_key(fptr, TSTRING, "ORIGIN", &origin, "ORIGIN", &status) )
         printerror( status );
    if ( fits_write_date(fptr, &status) )
         printerror( status );

    if ( fits_close_file(fptr, &status) )                /* close the file */
         printerror( status );

    return 0;
}
