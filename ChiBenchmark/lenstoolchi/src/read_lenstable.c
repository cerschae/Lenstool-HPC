/*this function will read a binary file containing gnfw slope (alpha),
  x (r/rsc), kappa and dpl.  These values will have to be rescaled
  before being used by lenstool.  The function will return the number
  of complete entries made*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "wcs.h"
#include <dimension.h>
#include <structure.h>

extern lensdata *lens_table;

void read_lenstable()
{
    long int i, tot_count;
    char infile[200];
    char *dir;
    long int alphasteps, xsteps;
    lensdata array_tmp[2];  // Read the first 2 rows of lenstable[]

    FILE *inf;

    // BUG: The lenstool.tab datafile only needs to exist in ONE place -
    // preferably the table_src directory. lenstool needs to use environment
    // variables for this!! Suggest:

    dir = getenv("LENSTOOL_DIR");
    if ( dir == NULL )
    {
        fprintf( stderr, "ERROR: LENSTOOL_DIR environment variable is not defined. Unable to find lenstool.tab file\n");
        exit(-1);

    }
    strcpy(infile, dir);
    strcat(infile, "/table_src/lenstool.tab");
    inf = fopen(infile, "rb+");
    if ( inf == NULL )
    {
        fprintf( stderr, "ERROR: Unable to find %s\n", infile);
        exit(-1);
    }

    fprintf(stderr, "READ: lens_table binary file into array...");

    //read in the first two rows (0:alpha,1:x) to find out how many more lines to read
    fread(&array_tmp[0], sizeof(lensdata), 1, inf);
    alphasteps = (long int)array_tmp[0].x_now;

    fread(&array_tmp[1], sizeof(lensdata), 1, inf);
    xsteps = (long int)array_tmp[1].x_now;

    tot_count = xsteps * alphasteps + 2;

    if (tot_count > ARRAY_SIZE)
    {
        fprintf(stderr, "\nREAD ERROR: lens_table, Need to increase output array size!!!\n");
        exit(-1);
    }

    // Allocate lens_table
    lens_table = (lensdata*) malloc( (unsigned int) tot_count * sizeof(lensdata) );

    // Fill the first 2 rows of lens_table with array_tmp
    memcpy( &lens_table[0], &array_tmp[0], sizeof(lensdata) );
    memcpy( &lens_table[1], &array_tmp[1], sizeof(lensdata) );

    for (i = 2; i < tot_count; i++)
        fread(&lens_table[i], sizeof(lensdata), 1, inf);

    fprintf(stderr, "done\n");
    fprintf(stderr, "INFO: NFWg limits (alpha[%lf,%lf], r/rs[%lf,%lf])\n",
            lens_table[2].alpha_now, lens_table[tot_count-1].alpha_now,
            lens_table[2].x_now, lens_table[tot_count-1].x_now );

    fclose(inf);

}
