#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "fitsio.h"
#include "lt.h"

/*
Simple write binary table fits subroutines based on the CFITSIO package

Author: EJ
Date: March 2010


 ********************************************************************
 * Create a FITS primary data unit containing a binary table
 * Write in double array
 * Parameters :
 * - filename 
 * - tbl[nb_col][nb_row] : image to be written
 */
int wrf_tfits(char *filename,
			 double **tbl, char* ttype[], int nx, long int ny)
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, ii, jj;
    LONGLONG  frow, felem;
    float bzero=0.,bscale=1.;
    char  origin[10];

    // initialize FITS table parameters 
    int tbltype=BINARY_TBL;
    LONGLONG naxis2=0;   // initial number of rows  
    int  tfields=nx; 
    char **tform;
    char chr;

    // Create the TFORM list
    tform = (char **) malloc(nx * sizeof(char *));
    for( ii = 0; ii < nx; ii++ )
    {
        tform[ii] = (char *) malloc(8 * sizeof(char));
        strcpy(tform[ii], "E");
    }

    // Replace non FITS standard characters in TTYPE list [0-9] [A-Z] [a-z]
    for( ii = 0; ii < nx; ii++ )
    {
        chop(ttype[ii]);  // remove trailing blanks 
        for( jj = 0; jj < strlen(ttype[ii]); jj++ )
        {
            chr = ttype[ii][jj];
            if( ! ((chr>=48 && chr<=57) || (chr>=65 && chr<=90) || (chr>=97 && chr<=122)) )
                ttype[ii][jj] = '_';
        }
    }

	printf("WRITE: TFITS file: %s[%d,%ld]\n", filename,nx,ny);

    remove(filename);               // Delete old file if it already exists 

    status = CREATE_DISK_FILE;         // initialize status before calling fitsio routines
    if (fits_create_file(&fptr, filename, &status)) // create new FITS file 
         printerror( status );           // call printerror if error occurs 

    if ( fits_create_tbl(fptr,  tbltype, naxis2, tfields, ttype, tform, NULL, NULL, &status) )
         printerror( status );          

    // Write columns
    for( ii = 0; ii < nx; ii++ )
    {
        frow = 1;  // first row to write
        felem = 1;  // first element to write
        if( fits_write_col(fptr, TDOUBLE, ii+1, frow, felem, ny, tbl[ii], &status) )
            printerror( status );
    }

	strcpy(origin,"JPK-SOFT");

    /* write another optional keyword to the header */
    if ( fits_update_key(fptr, TFLOAT, "BZERO", &bzero, "BZERO", &status) )
         printerror( status );           
    if ( fits_update_key(fptr, TFLOAT, "BSCALE", &bscale, "BSCALE", &status) )
         printerror( status );           
    if ( fits_update_key(fptr, TSTRING, "ORIGIN", &origin, "ORIGIN", &status) )
         printerror( status );           
    if ( fits_write_date(fptr, &status) )
         printerror( status );           

    if ( fits_close_file(fptr, &status) )                /* close the file */
         printerror( status );           

    return 0;
}
