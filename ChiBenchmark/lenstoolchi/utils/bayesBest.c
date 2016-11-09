#include<stdio.h>
#include<signal.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include "dimension.h"
#include "structure.h"
#include "constant.h" 
#include "fonction.h"
#include "bayesChires.h"

/******************************************************************************
 * Name: bayesBest.c
 * Authors: EJ
 * Date: 31/08/07
 * 
 * Analyse a bayes.dat file and compute the best.par and bestopt.par files
 * from it. 
 *
 * syntax : bayesBest <.par> 
 *
 * The bayes.dat file in the current directory is used.
 * Output : best.par bestopt.par
 *
 ******************************************************************************/
typedef void (*sighandler_t)(int);
int optInterrupt;

void help_msg()
{
		fprintf(stderr, "Syntax : bayesBest [OPTION] <.par> <method>\n");
		fprintf(stderr, "Available methods : best, mean, median, mode\n");
		fprintf(stderr, "\t-u\n\t\tset uniform prior for the bestopt.par limits.\n");
		fprintf(stderr, "\t-g\n\t\tset Gaussian prior for the bestopt.par limits (Default).\n");
		fprintf(stderr, "\t-s[LIMIT]\n\t\tset the width of the prior to LIMIT x 1 sigma (Default 3).\n");
		exit(-1);
}

void writeSource()
{

}


int main( int argc, char** argv )
{
	int    method;
	int    ival, i;
	int    prior;
	double limit;
	char   *pch;
	FILE   *IN;
    extern int sblock[NFMAX][NPAMAX];
    extern struct g_source  S;
    extern struct galaxie  sources[NFMAX];
    int ipx;

	// Check the arguments
	if( argc < 2 )
		help_msg();

	// Analyse the options
	prior	= 3;
	limit = 3.;
	while( argv[1][0] == '-' )
	{
		pch = argv[1];
		if( pch[1] == 'u'  )
			prior = 1; // uniform prior
		else if( pch[1] ==  's' )
		{	
			pch+=2;
			if( sscanf(pch, "%lf", &limit) == 0) help_msg();
		}

		for( i = 1; i < argc-1 ; i++ )
			argv[i]=argv[i+1];
		argc--;
	}	

	if( strstr(argv[1], ".par") == NULL )
		help_msg();

	// Default method : best
	method = -4;
	if( argc == 3 )
	{
		if( !strcmp( argv[2], "mean" ) ) method = -1;
		else if( !strcmp( argv[2], "median" ) ) method = -2;
		else if( !strcmp( argv[2], "mode" ) ) method = -3;
		else if( sscanf( argv[2], "%d", &ival ) == 1 ) method = ival;
	}

	// Read the .par file
	init_grille( argv[1], 1);
	
	// Read constraints
	readConstraints();

	if( G.nmsgrid != G.nlens ) 
	{
		prep_non_param();
	}
	
	// Set the best model in memory
	o_set_lens_bayes(method, prior, limit);

	// Reset the lens parameters
	set_res_par();

	// Print the best.par and bestopt.par files
	// save the original best.par and bestopt.par
	if( IN =  fopen( "best.par" ,"r") )
	{
		fclose(IN);
		time_t rawtime;
		struct tm *timeinfo;
		char   buffer[15],cmd[50];
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		strftime(buffer,15,"%Y%m%d%H%M%S",timeinfo);
		sprintf(cmd,"cp best.par best.par.%s",buffer);
		system(cmd);
		sprintf(cmd,"cp bestopt.par bestopt.par.%s",buffer);
		system(cmd);
	}

    // Write a new source.best file if source optimization
    int sourcelimit_flag = 0;
    for( i = 0; i < S.ns; i++ )
        for( ipx = SCX; ipx <= SFLUX; ipx ++ )
            sourcelimit_flag += sblock[i][ipx];

    if( sourcelimit_flag > 0 )
        strcpy(M.sourfile, "source.best" );

    M.source = 2;
	o_print_res(o_chi(), 0.);  // we don't know the evidence

    // clean o_global() variables... put source coordinates back to absolute
    long int sns_sav = S.ns; 
    o_global_free();
    S.ns = sns_sav;

    if( sourcelimit_flag > 0 )
        ecrire_r(0, S.ns, source, "source.best", 6);

	return 0;
}

