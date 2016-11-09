#include<stdio.h>
#include<signal.h>
#include<stdlib.h>
#include<string.h>
#include "dimension.h"
#include "structure.h"
#include "constant.h" 
#include "fonction.h"
#include "bayesChires.h"

/******************************************************************************
 * Name: bayesDlsDs.c
 * Authors: EJ
 * Date: 30/10/07
 * 
 * Analyse a bayes.dat file and print the Dls/Dos ratio for each optimised 
 * system.
 *
 * syntax : bayesDlsDs <.par> 
 *
 * The bayes.dat file in the current directory is used.
 * Output : none
 *
 ******************************************************************************/
typedef void (*sighandler_t)(int);
static void signalReset();
int optInterrupt;

int main( int argc, char** argv )
{
	double **array; // contains the bayes.dat data
	long int nVal;  // size of array
    int nParam;
	int iVal, i;

	// Check the arguments
	if( argc < 2 || 
	    strstr( argv[1], ".par" ) == NULL )
	{
		fprintf(stderr, "Syntax : bayesDlsDs <.par>\n");
		return -1;
	}

	// Read the .par file
	init_grille( argv[1], 1);
	
	// Read the image file
	readConstraints();

	// Print header with redshifts
	printf( "#Nsample\n#OmegaM\n#wX\n" );
	for( i = 0; i < I.n_mult; i++ )
	{
		printf( "#Redshift of %s : %f\n", 
			multi[i][0].n, multi[i][0].z );
	}

	// Read the bayes.dat file
	array = readBayesModels(&nParam, &nVal);
	if( array == NULL )
	{
		fprintf(stderr, "ERROR: bayes.dat file not found\n");
		return -1; 
	}

	// Handle CTRL-C to interrupt the optimisation but not lenstool
	signal(SIGINT, signalReset);
	optInterrupt = 0;

	// Loop over each line
	for( iVal = 0; iVal < nVal && !optInterrupt; iVal++ )
	{
		// Set the lens parameters from <array>
		setBayesModel( iVal, nVal, array );

		printf( "%d %f %f ", iVal, C.omegaM, C.wX );
		// Print the Dls/Ds ratio for each system
		for( i = 0; i < I.n_mult; i++ )
		{
			printf( "%f ", dratio( lens[0].z, multi[i][0].z ) );
		}
		printf("\n");
	}
	free( array );
	return 0;
}


static void signalReset()
{
	signal(SIGINT, SIG_DFL);
	optInterrupt = 1;
	printf( "\nINFO: Optimisation interrupted by CTRL-C\r");
}
