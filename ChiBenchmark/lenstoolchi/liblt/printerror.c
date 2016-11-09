/********************************************************************/
/* Print out cfitsio error messages and exit program */
/*****************************************************/
#include<stdlib.h>
#include "fitsio.h"
#include "lt.h"

void printerror( int status)
{

    if (status)
    {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
}
