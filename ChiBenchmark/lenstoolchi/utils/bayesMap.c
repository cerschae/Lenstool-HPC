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
 * Name: bayesMap.c
 * Authors: EJ
 * Date: 12/09/07
 * 
 * Analyse a bayes.dat file and for each line, compute the maps specified
 * in the runmode section.
 *
 * syntax : bayesMap <.par> 
 *
 * Input : The bayes.dat file in the current directory is used.
 * Output in tmp/ directory
 *
 * Existing files are skipped.
 *
 ******************************************************************************/
typedef void (*sighandler_t)(int);
static void signalReset();
int optInterrupt;
struct point gsource_global[NGGMAX][NGGMAX]; //(for g_ampli!)

int main( int argc, char** argv )
{
	double **array; // contains the bayes.dat data
	int nParam;
    long int iVal, nVal;  // size of array
	char fname[50]; // <map><ival>.fits
	char fname2[50]; // <map><ival>.fits
	FILE *pFile;
	int i;
	double *index;  // list of bayes.dat lines
	int    seed;   // random seed
	int tmp;
	
	// Check the arguments
	if( argc < 2 || 
	    strstr( argv[1], ".par" ) == NULL )
	{
		fprintf(stderr, "Syntax : bayesMap <.par>\n");
		return -1;
	}

	// Read the .par file
	init_grille( argv[1], 1);

	// remove the .fits extension to filename
	if( M.imass ) M.massfile[strlen(M.massfile)-5]=0;
	if( M.ishear ) M.shearfile[strlen(M.shearfile)-5]=0;
	if( M.iampli ) M.amplifile[strlen(M.amplifile)-5]=0;
	if( M.idpl ) 
	{
	       M.dplxfile[strlen(M.dplxfile)-5]=0;
	       M.dplyfile[strlen(M.dplyfile)-5]=0;
	}
    if( M.pixel ) M.pixelfile[strlen(M.pixelfile)-5]=0;
    if( M.iclean ) ps.pixfile[strlen(ps.pixfile)-5]=0;
	
	// Read catalog of multiple images 
	readConstraints();

	// Initialise the grid
	if( G.pol != 0 )
		gridp();
	else
		grid();
	
	// Switch to silent mode
	M.verbose = 0;

	//if( G.nmsgrid != G.nlens ) 
	//	prep_non_param();
	
	// Read the bayes.dat file
	array = readBayesModels(&nParam, &nVal);
	if( array == NULL )
	{
		fprintf(stderr, "ERROR: bayes.dat file not found\n");
		return -1; 
	}

	// Create the ./tmp directory
	i = system("mkdir -p tmp");
	
	// Prepare the index list
	index = (double *) malloc((unsigned) nVal*sizeof(double));
	for( i = 0 ; i < nVal ; i++ ) index[i]=i;
	seed = -2; 

	// Handle CTRL-C to interrupt the optimisation but not lenstool
	signal(SIGINT, signalReset);
	optInterrupt = 0;

	// Loop over each line
	for( i = 0; i < nVal && !optInterrupt && i < 2000; i++ )
	{
		// Randomly draw a line from index array 
		tmp = (int) floor(d_random(&seed) * (nVal - i));
		iVal = index[i+tmp];
		// and swap the indexes in the index list
		index[i+tmp] = index[i]; 		
		
		// Set the lens parameters from <array>
		setBayesModel( iVal, nVal, array );

  //   	if( G.nmsgrid != G.nlens ) 
  //          rhos2b0();


		if( M.imass != 0 )
		{
			sprintf( fname, "tmp/%s%ld.fits",M.massfile, iVal );
			printf("INFO: Compute file %d/%ld : %s [CTRL-C to interrupt]\n",i+1, nVal,fname);
			fflush(stdout);
			pFile = fopen( fname, "r" );
			if( pFile == NULL ) 
			{
				pFile = fopen( fname, "w");
				fprintf( pFile, "busy\n" );
				fclose(pFile);
				g_mass( M.imass, M.nmass, M.zmass, S.zs, fname );
			}
			else
				fclose(pFile);
		}

		if( M.ishear != 0 )
		{
			sprintf( fname, "tmp/%s%ld.fits",M.shearfile, iVal );
			printf("INFO: Compute file %d/%ld : %s [CTRL-C to interrupt]\n",i+1, nVal,fname);
			fflush(stdout);
			pFile = fopen( fname, "r" );
			if( pFile == NULL ) 
			{
				pFile = fopen( fname, "w");
				fprintf( pFile, "busy\n" );
				fclose(pFile);
				g_shear( M.ishear, M.nshear, M.zshear, fname );
			}
			else
				fclose(pFile);
		}

		if( M.idpl != 0 )
		{
			sprintf( fname, "tmp/%s%ld.fits",M.dplxfile, iVal );
			sprintf( fname2, "tmp/%s%ld.fits",M.dplyfile, iVal );
			printf("INFO: Compute file %d/%ld : %s %s [CTRL-C to interrupt]\n",i+1, nVal,fname,fname2);
			fflush(stdout);
			pFile = fopen( fname, "r" );
			if( pFile == NULL ) 
			{
				pFile = fopen( fname, "w");
				fprintf( pFile, "busy\n" );
				fclose(pFile);
				g_dpl( M.idpl, M.ndpl, M.zdpl,fname,fname2 );
			}
			else
				fclose(pFile);
		}

		if( M.iampli != 0 )
		{
			sprintf( fname, "tmp/%s%ld.fits",M.amplifile, iVal );
			printf("INFO: Compute file %d/%ld : %s [CTRL-C to interrupt]\n",i+1, nVal,fname);
			fflush(stdout);
			pFile = fopen( fname, "r" );
			if( pFile == NULL ) 
			{
				pFile = fopen( fname, "w");
				fprintf( pFile, "busy\n" );
				fclose(pFile);
				g_ampli( M.iampli, M.nampli, M.zampli, fname );
			}
			else
				fclose(pFile);
		}

		if( M.pixel != 0 )
		{
			sprintf( fname, "tmp/%s%ld.fits",M.pixelfile, iVal );
			sprintf( fname2, "tmp/%s_src%ld.fits",M.pixelfile, iVal );
			printf("INFO: Compute file %d/%ld : %s [CTRL-C to interrupt]\n",i+1, nVal,fname);
			fflush(stdout);
			pFile = fopen( fname, "r" );
			if( pFile == NULL ) 
			{
				pFile = fopen( fname, "w");
				fprintf( pFile, "busy\n" );
				fclose(pFile);
				e_pixel( M.npixel, fname, fname2, source);
			}
			else
				fclose(pFile);
		}
		if( M.iclean == 1 )
		{
			sprintf( fname, "tmp/%s%ld.fits",ps.pixfile, iVal );
			printf("INFO: Compute file %d/%ld : %s [CTRL-C to interrupt]\n",i+1, nVal,fname);
			fflush(stdout);
			pFile = fopen( fname, "r" );
			if( pFile == NULL ) 
			{
				pFile = fopen( fname, "w");
				fprintf( pFile, "busy\n" );
				fclose(pFile);
				imtosou( M.zclean, fname);
			}
			else
				fclose(pFile);
		}
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
