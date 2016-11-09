#include<stdio.h>
#include<signal.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "dimension.h"
#include "structure.h"
#include "constant.h"
#include "fonction.h"
#include "bayesChires.h"

/******************************************************************************
 * Name: bayes_chires.c
 * Authors: EJ
 * Date: 31/08/07
 *
 * Analyse a bayes.dat file and compute a chires.dat file for each line.
 * Finally, it produces a mean chires.dat file from all individual chires.dat
 * file.
 *
 * syntax : bayesChires <.par>
 *
 * The bayes.dat file in the current directory is used.
 * Output : chires*.dat
 *
 ******************************************************************************/
typedef void (*sighandler_t)(int);
static void signalReset();
int optInterrupt;

int main( int argc, char** argv )
{
    double **array; // contains the bayes.dat data
    int nParam;
    long int  nVal;  // size of array
    long int i, iVal;
    char fname[50]; // chires<ival>.dat
    long int *index;  // list of bayes.dat lines
    int    seed;   // random seed
    long int tmp;
    FILE   *pFile;

    // Check the arguments
    if ( argc < 2 ||
         strstr( argv[1], ".par" ) == NULL )
    {
        fprintf(stderr, "Syntax : bayes_chires <.par>  [bayes.dat]\n");
        return -1;
    }

    // Read the .par file
    init_grille( argv[1], 1);

    // Read the image file
    readConstraints();

    if ( G.nmsgrid != G.nlens )
    {
        prep_non_param();
    }

    // Read the bayes.dat file
    array = readBayesModels(&nParam, &nVal);
    if ( array == NULL )
    {
        fprintf(stderr, "ERROR: bayes.dat file not found\n");
        return -1;
    }

    // Create the ./tmp directory
    i = system("mkdir -p chires");

    // Prepare the index list
    index = (long int *) malloc((unsigned) nVal * sizeof(double));
    for ( i = 0 ; i < nVal ; i++ ) index[i] = i;
    seed = -2;

    // Handle CTRL-C to interrupt the optimisation but not lenstool
    signal(SIGINT, signalReset);
    optInterrupt = 0;

    // Loop over each line
    for ( i = 0; i < nVal && !optInterrupt; i++ )
    {
        // Randomly draw a line from index array
        tmp = (long int) floor(d_random(&seed) * (nVal - i));
        iVal = index[i + tmp];
        // and swap the indexes in the index list
        index[i+tmp] = i;

        // skip file is already processed
        sprintf( fname, "chires/chires%ld.dat", iVal );

        pFile = fopen(fname, "r");
        if( pFile == NULL )
        {
            // Set the lens parameters from <array>
            setBayesModel( iVal, nVal, array );
            printf( "INFO: Compute the chires.dat file %ld/%ld  [CTRL-C to interrupt]\r", i + 1, nVal);
            fflush(stdout);
    
            // Execute the chires()  --> chires<iVal>.dat
            o_chires(fname);
        }
        else
            fclose(pFile);

    }
    printf( "\n" );
    free( array );
    return 0;
}


static void signalReset()
{
    signal(SIGINT, SIG_DFL);
    optInterrupt = 1;
    printf( "\nINFO: Optimisation interrupted by CTRL-C\r");
}
