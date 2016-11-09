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
 * Name: bayesAmp.c
 * Authors: EJ
 * Date: 30/10/07
 * 
 * Analyse a bayes.dat file and print the amplification for each optimised 
 * system.
 *
 * syntax : bayesAmp <.par> 
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
	double amp;
	int nParam;  // size of array
	long int iVal, nVal; 
    int i, j;

	// Check the arguments
	if( argc < 2 || 
	    strstr( argv[1], ".par" ) == NULL )
	{
		fprintf(stderr, "Syntax : bayesAmp <.par>\n");
		return -1;
	}

	// Read the .par file
	init_grille( argv[1], 1);
	
	// Read the image file
	readConstraints();

	// Print header with redshifts
	printf( "#Nsample\n" );
	for( i = 0; i < I.n_mult; i++ )
	{
		printf( "#Amplification of %s : %f\n", 
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

		printf( "%ld ", iVal );
		// Print the averaged amplification for each system
		for( i = 0; i < I.n_mult; i++ )
		{
			amp = 0.;
			for( j = 0; j < I.mult[i]; j++ )
				amp += 1./e_amp(&multi[i][j].C,multi[i][j].dl0s,multi[i][j].dos,multi[i][j].z);

			amp /= I.mult[i];
			printf( "%f ", amp );
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
