#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "fitsio.h"
#include "lt.h"


/*
Author: JPK
Date: April 2001
Read into a double array:
rdf_fits(filename,ima,nx,ny,xmin,xmax,ymin,ymax)


double ** rdf_fits(filename,nx,ny,xmin,xmax,ymin,ymax)
int     *nx,*ny;
double   *xmin,*xmax,*ymin,*ymax;
char *filename;

{
	struct matrix {
		double	a;
		double	b;
		double	c;
		double	d;
	} cd;
	double	cdelt[2];
	double	crpix[2];
	double	crval[2];
	char	ctype[2][4];
	
	
	return(rdf_fits_abs(filename,nx,ny,xmin,xmax,ymin,ymax,cd,cdelt,crpix,crval,ctype)
}
*/
/********************************************************************/
/* Simple read double fits subroutine based on the CFITSIO package
 * Return an array[ny][nx] containing the image in FITS file filename
 * xmin, xmax, ymin and ymax are computed from the header keywords
 */
/*double ** rdf_fits(filename,nx,ny,cd,cdelt,crpix,crval,ctype)
int     *nx,*ny;
double   *cdelt,*crpix,*crval;
char	**ctype; 
struct	matrix {
	double	a;
	double	b;
	double	c;
	double	d;
}	*cd;
*/

static int rdf_get_astro(fitsfile *fptr, char **pheader);


double ** rdf_fits(char *filename,int *nx,int *ny,char **pheader)
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, ii, jj, kk;
    int nfound, anynull;
    long naxes[2];   
    long fpixel, nbuffer, npixels;
	int	nkeys;
//    double crval1,crpix1,cdelt1;
//    double crval2,crpix2,cdelt2;
    double **image;
#define buffsize 1000
    double datamin, datamax, nullval;
    float buffer[buffsize];

    status = 0;         /* initialize status before calling fitsio routines */

    if (fits_open_file(&fptr, filename, READONLY, &status)) 
         printerror( status );           /* call printerror if error occurs */

    /* read the NAXIS1 and NAXIS2 keyword to get image size */
    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
         printerror( status );          


	fpixel = 1;                             /* first pixel in the image    */
	npixels = naxes[0] * naxes[1];          /* number of pixels in the image */
	*nx = naxes[0];
	*ny = naxes[1];
	
	/*image[nb_lin][nb_col]*/
	image =(double **) alloc_square_double(*ny,*nx);
	ii=jj=0;
	
	nullval  = 0;                /* don't check for null values in the image */
	datamin  = 1.0E30;
	datamax  = -1.0E30;
	
		
	while (npixels > 0)
	{
		nbuffer = npixels;
		if (npixels > buffsize)
			nbuffer = buffsize;     /* read as many pixels as will fit in buffer */
	
	  /* Note that even though the FITS images contains unsigned integer */
	  /* pixel values (or more accurately, signed integer pixels with    */
	  /* a bias of 32768),  this routine is reading the values into a    */
	  /* double array.   Cfitsio automatically performs the datatype      */
	  /* conversion in cases like this.                                  */
	
		if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval,
				buffer, &anynull, &status) )
			printerror( status );

		for (kk = 0; kk < nbuffer; kk++)  
		{
			image[jj][ii++]=buffer[kk];
			if (ii>= *nx)
			{	jj++; ii=0; }
	
			if ( buffer[kk] < datamin )
				datamin = buffer[kk];
	
			if ( buffer[kk] > datamax )
				datamax = buffer[kk];
		}
		
		npixels -= nbuffer;    /* increment remaining number of pixels */
		fpixel  += nbuffer;    /* next pixel to be read in image */
	}
	
	
	printf("Min and max image pixels =  %lf, %lf\n", datamin, datamax);
	
	if( fits_hdr2str( fptr, 0, NULL, 0, pheader, &nkeys, &status ) )
          printerror( status );
	
	if ( fits_close_file(fptr, &status) )                /* close the file */
	         printerror( status );           
	
//	printf("%.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n",crval[0],cdelt[0],crpix[0],
//												crval[1],cdelt[1],crpix[1]);
//	printf("%.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n",crval1,cdelt1,crpix1,
//												crval2,cdelt2,crpix2);												
/*	*xmin=crval1 - (crpix1 -1)*(cdelt1);
	*xmax=crval1 + (*nx-crpix1)*(cdelt1);
	*ymin=crval2 - (crpix2 -1)*(cdelt2);
	*ymax=crval2 + (*ny-crpix2)*(cdelt2);
	
	printf("%.3lf %.3lf %.3lf %.3lf \n", *xmin, *xmax, *ymin, *ymax);
	*/
	return ((double **) image);
}

double *** rdf_cube_fits(char *filename,int *nx,int *ny,int *nz, char **pheader)
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, ii, jj, kk, zz;
    int nfound, anynull;
    long naxes[3];   
    long fpixel, nbuffer, npixels;
	int	nkeys;
//    double crval1,crpix1,cdelt1;
//    double crval2,crpix2,cdelt2;
    double ***cubeimage;
#define buffsize 1000
    double datamin, datamax, nullval;
    float buffer[buffsize];

    status = 0;         /* initialize status before calling fitsio routines */

    if (fits_open_file(&fptr, filename, READONLY, &status)) 
         printerror( status );           /* call printerror if error occurs */

    /* read the NAXIS1 NAXIS2 and NAXIS3 keywords to get image size */
    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 3, naxes, &nfound, &status) )
         printerror( status );          


	fpixel = 1;                             /* first pixel in the cube image */
	npixels = naxes[0] * naxes[1] * naxes[2];          /* number of pixels in the cube image */
	*nx = naxes[0];
	*ny = naxes[1];
	*nz = naxes[2];
	
	/*image[nb_lin][nb_col]*/
	cubeimage =(double ***) alloc_cubic_double(*ny,*nx,*nz);
	ii=jj=0;
	
	nullval  = 0;                /* don't check for null values in the image */
	datamin  = 1.0E30;
	datamax  = -1.0E30;
	
		
	while (npixels > 0)
	{
		nbuffer = npixels;
		if (npixels > buffsize)
			nbuffer = buffsize;     /* read as many pixels as will fit in buffer */
	
	  /* Note that even though the FITS images contains unsigned integer */
	  /* pixel values (or more accurately, signed integer pixels with    */
	  /* a bias of 32768),  this routine is reading the values into a    */
	  /* double array.   Cfitsio automatically performs the datatype      */
	  /* conversion in cases like this.                                  */
	
		if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval,
				buffer, &anynull, &status) )
			printerror( status );

		for (kk = 0; kk < nbuffer; kk++)  
		{
			cubeimage[jj][ii][zz++]=buffer[kk];
                        if (zz>= *nz)
                        {       ii++; zz=0; }
			if (ii>= *nx)
			{	jj++; ii=0; }
	
			if ( buffer[kk] < datamin )
				datamin = buffer[kk];
	
			if ( buffer[kk] > datamax )
				datamax = buffer[kk];
		}
		
		npixels -= nbuffer;    /* increment remaining number of pixels */
		fpixel  += nbuffer;    /* next pixel to be read in image */
	}
	
	
	printf("Min and max image pixels =  %lf, %lf\n", datamin, datamax);
	
	if( fits_hdr2str( fptr, 0, NULL, 0, pheader, &nkeys, &status ) )
          printerror( status );
	
	if ( fits_close_file(fptr, &status) )                /* close the file */
	         printerror( status );           
	
//	printf("%.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n",crval[0],cdelt[0],crpix[0],
//												crval[1],cdelt[1],crpix[1]);
//	printf("%.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n",crval1,cdelt1,crpix1,
//												crval2,cdelt2,crpix2);												
/*	*xmin=crval1 - (crpix1 -1)*(cdelt1);
	*xmax=crval1 + (*nx-crpix1)*(cdelt1);
	*ymin=crval2 - (crpix2 -1)*(cdelt2);
	*ymax=crval2 + (*ny-crpix2)*(cdelt2);
	
	printf("%.3lf %.3lf %.3lf %.3lf \n", *xmin, *xmax, *ymin, *ymax);
	*/
	return ((double ***) cubeimage);
}







/* Retrieve the astrometric information contained in a header of a FITS file
 * Return -1 if no astrometric information is present, 0 otherwise
 * */
static int rdf_get_astro(fitsfile *fptr, char **pheader)
{
	int	status,kk;
    int nkeys, keypos;
    char card[FLEN_CARD];
    char header[1500];
	double fval;
	printf("todo bien\n");
	/* read keywords */
	if (fits_get_hdrpos(fptr, &nkeys, &keypos, &status) )
		printerror( status );
	printf("todo bien\n");	
	for (kk = 1; kk <= nkeys; kk++)  
	{
		if ( fits_read_record(fptr, kk, card, &status) )
			printerror( status );

		printf("%s\n", card); /* print the keyword card */
	
		if (strncmp(card,"CRVAL1",6)==0)
		{
			sscanf(card, "CRVAL1  = %lf", &fval);
			strcat(header, card);
//	    	crval[0]=fval;
//			crval1=fval;
	    }
		else if (strncmp(card,"CRPIX1",6)==0)
		{
			sscanf(card, "CRPIX1  = %lf", &fval);
			strcat(header, card);
//	    	crpix[0]=fval;
//			crpix1=fval;
	    }
		else if (strncmp(card,"CDELT1",6)==0)
		{
			sscanf(card, "CDELT1  = %lf", &fval);
			strcat(header, card);
//	    	cdelt[0]=fval;
//			cdelt1=fval;
	    }
		else if (strncmp(card,"CRVAL2",6)==0)
		{
			sscanf(card, "CRVAL2  = %lf", &fval);
			strcat(header, card);
//	    	crval[1]=fval;
//			crval2=fval;
	    }
		else if (strncmp(card,"CRPIX2",6)==0)
		{
			sscanf(card, "CRPIX2  = %lf", &fval);
			strcat(header, card);
//	    	crpix[1]=fval;
//			crpix2=fval;
	    }
		else if (strncmp(card,"CDELT2",6)==0)
		{
			sscanf(card, "CDELT2  = %lf", &fval);
			strcat(header, card);
//	    	cdelt[1]=fval;
//			cdelt2=fval;
	    }
		else if (strncmp(card,"CD1_1",6)==0)
		{
			sscanf(card, "CD1_1  = %lf", &fval);
			strcat(header, card);
//	    	cd[0][0]=fval;
	   }
		else if (strncmp(card,"CD1_2",6)==0)
		{
			sscanf(card, "CD1_2  = %lf", &fval);
			strcat(header, card);
//	    	cd[0][1]=fval;
	    }
		else if (strncmp(card,"CD2_1",6)==0)
		{
			sscanf(card, "CD2_1  = %lf", &fval);
			strcat(header, card);
//	    	cd[1][0]=fval;
	    }
		else if (strncmp(card,"CD2_2",6)==0)
		{
			sscanf(card, "CD2_2  = %lf", &fval);
			strcat(header, card);
//	    	cd[1][1]=fval;
	    }
		else if (strncmp(card,"PC1_1",6)==0)
		{
			sscanf(card, "PC1_1  = %lf", &fval);
			strcat(header, card);
//	    	cd[0][0]=fval;
	    }
		else if (strncmp(card,"PC1_2",6)==0)
		{
			sscanf(card, "PC1_2  = %lf", &fval);
			strcat(header, card);
//	    	cd[0][1]=fval;
	    }
		else if (strncmp(card,"PC2_1",6)==0)
		{
			sscanf(card, "PC2_1  = %lf", &fval);
			strcat(header, card);
//	    	cd[1][0]=fval;
	    }
		else if (strncmp(card,"PC2_2",6)==0)
		{
			sscanf(card, "PC2_2  = %lf", &fval);
			strcat(header, card);
//	    	cd[1][1]=fval;
	    }	
		else if (strncmp(card,"CRTYPE1",6)==0)
		{
			strcat(header, card);
//			sscanf(card, "CRTYPE1  = %s", &sval);
//	    	strcpy(crtype[0],sval);
	    }
		else if (strncmp(card,"CRTYPE2",6)==0)
		{
			strcat(header, card);
//			sscanf(card, "CRTYPE2  = %s", &sval);
//	    	strcpy(crtype[1],sval);
	    }
	}

/*	pcd=(long int)&cd;
	pcdelt=(long int)&cdelt;
	pcrpix=(long int)&crpix;
	pcrval=(long int)&crval;
	pcrtype=(long int)&crtype;*/
	
	*pheader = header;
	return(0); 
}
