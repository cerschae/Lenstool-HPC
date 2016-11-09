#include<stdio.h>
#include<signal.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "dimension.h"
#include "structure.h"
#include "constant.h" 
#include "fonction.h"
#include "bayesChires.h"
#include "lt.h"

/******************************************************************************
 * Name: bayesSource.c
 * Authors: EJ
 * Date: 12/09/07
 * 
 * Analyse a bayes.dat file and compute a source.dat file for each line. 
 *
 * syntax : bayesSource <.par> 
 *
 * The bayes.dat file in the current directory is used.
 * Output : source??.dat
 *
 ******************************************************************************/
typedef void (*sighandler_t)(int);
static void signalReset();
static void barysource(struct galaxie *source,long int ns, struct point *bsource);

int optInterrupt;

int main( int argc, char** argv )
{
    extern long int narclet;
	struct point bsource[100];
	double **array; // contains the bayes.dat data
	int nParam;     // size of array 
    long int iVal, nVal;  
	char fname[40]; // source<ival>.dat
	int nbs; // number of systems (barycenter source)
	int i,j;
	double dx, dy;
	FILE *bayes;

	// Check the arguments
	if( argc < 2 || 
	    strstr( argv[1], ".par" ) == NULL )
	{
		fprintf(stderr, "Syntax : bayesSource <.par>\n");
		return 1;
	}

	// Read the .par file
	init_grille( argv[1], 1);
	
	// Read catalog of multiple images 
	if( M.imafile[0] == '\0' ) 
	{
		fprintf(stderr, "ERROR: keyword image not defined\n");
		return 1;
	}
	narclet = 0;
	f_shape( &narclet, arclet, M.imafile, 0 );
	pro_arclet(narclet, arclet);
    for( i = 0; i < narclet; i++ )
        dratio_gal(&arclet[i], lens[0].z);

	// Initialise the grid
	if( G.pol != 0 )
		gridp();
	else
		grid();
	
	// Switch to silent mode
	M.verbose = 0;
	M.image = 1; // to cheat the e_lensing print out

	if( G.nmsgrid != G.nlens ) 
	{
		prep_non_param();
	}
	
	// Read the bayes.dat file
	array = readBayesModels(&nParam, &nVal);
	if( array == NULL )
	{
		fprintf(stderr, "ERROR: bayes.dat file not found\n");
		return 1; 
	}

	i = system("mkdir -p source");
	
	// Write the header of the bayesSource.dat file
	bayes = fopen( "bayesSource.dat", "w" );

	fprintf( bayes, "#Nsample\n");
	fprintf( bayes, "#Chi2\n");
	for( i = 0; i < narclet; i++ )
		fprintf( bayes, "#Distance(I_%s - <S>)\n", arclet[i].n);

	// Handle CTRL-C to interrupt the optimisation but not lenstool
	signal(SIGINT, signalReset);
	optInterrupt = 0;

	// Loop over each line
	for( iVal = 0; iVal < nVal && !optInterrupt; iVal++ )
	{
		// Set the lens parameters from <array>
		setBayesModel( iVal, nVal, array );
		printf( "INFO: Compute the source.dat file %ld/%ld  [CTRL-C to interrupt]\r", iVal+1, nVal);
		fflush(stdout);

    	// arclet --> source list
		S.ns = 0; e_unlens(narclet, arclet, &S.ns, source);

		// source list --> barycenter source list
		barysource(source, S.ns, bsource);

		// print a line in bayesSource.dat
		fprintf( bayes, "%ld %lf ", iVal, array[0][iVal] );
		for( i = 0;  i < narclet; i++ )
		{
			dx = arclet[i].C.x-bsource[i].x;
			dy = arclet[i].C.y-bsource[i].y;
			fprintf(bayes, "%lf ", sqrt( dx*dx + dy*dy)/arclet[i].dr );
		}
		fprintf( bayes, "\n" );
		fflush(bayes);
	
		// write source array
		sprintf( fname, "source/source%ld.dat", iVal );
		ecrire_r(0,S.ns,source,fname,2);
	}
	printf( "\n" );
	fclose(bayes);
	free( array );
	return 0;
}


static void signalReset()
{
	signal(SIGINT, SIG_DFL);
	optInterrupt = 1;
	printf( "\nINFO: Optimisation interrupted by CTRL-C\r");
}

// For each system return in bsource the barycenter source
static void barysource(struct galaxie *source, long int ns,
	       	struct point *bsource)
{
	int i,j;
	double xs, ys, n;

	xs = ys = 0.;
	n = 0.;
	for( i = 0; i < ns-1; i++ )
	{
		xs += source[i].C.x;
		ys += source[i].C.y;
		n++;

		if( indexCmp(source[i].n,source[i+1].n) )
		{
			for( j = i-n+1; j <= i; j ++) 
			{
				bsource[j].x = xs / n;
				bsource[j].y = ys / n;
			}
			n = 0.;
			xs = ys = 0.;
		}
	}
	
	// Fill the last system
	xs += source[i].C.x;
	ys += source[i].C.y;
	n++;
	for( j = ns-n; j < ns; j ++) 
	{
		bsource[j].x = xs / n;
		bsource[j].y = ys / n;
	}
	
}
