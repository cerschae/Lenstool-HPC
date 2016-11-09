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
 * Name: bayesImage.c
 * Authors: EJ
 * Date: 12/09/07
 * 
 * Analyse a bayes.dat file and compute a image.all file for each line. 
 *
 * syntax : bayesImage <.par> 
 *
 * The bayes.dat file in the current directory is used.
 * Output : image??.all
 *
 ******************************************************************************/
typedef void (*sighandler_t)(int);
static void signalReset();
static int bc_mode();
static int image_mode();
int optInterrupt;
struct point gsource_global[NGGMAX][NGGMAX]; //(for e_lensing)

void help_msg()
{
		fprintf(stderr, "Syntax : bayesImage [OPTION] <.par>\n");
        fprintf(stderr, "Available OPTIONS:\n");
        fprintf(stderr, " -b : compute images of system barycenter [default]\n");
        fprintf(stderr, " -i : compute counter images, image per image\n");
		exit(-1);
}

int main( int argc, char** argv )
{
	double **array; // contains the bayes.dat data
	int nParam;  // size of array
	long int iVal, nVal;
	char fname[128]; // image<ival>.all
    int mode;  // lens barycenter (0) or lens image per image (1)
    char *pch;
	int ni, status, i;

	// Check the arguments
	if( argc < 2 )
        help_msg();

    mode = 0;  // image of barycenter
	while( argv[1][0] == '-' )
	{
		pch = argv[1];
		if( pch[1] ==  'b' )
            mode = 0;
		if( pch[1] == 'i'  )
            mode = 1;

		for( i = 1; i < argc-1 ; i++ )
			argv[i]=argv[i+1];
		argc--;
	}	

	if( strstr( argv[1], ".par" ) == NULL )
        help_msg();

	// Read the .par file
	init_grille( argv[1], 1);

    // image per image computation reads runmode/image keyword
    if( mode == 1 && M.imafile[0] == '\0' )
    {
        fprintf(stderr, "ERROR: keyword runmode/image is not defined\n");
        exit(-1);
    }
	
	// Read catalog of multiple images 
	readConstraints();

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
		return -1; 
	}

    // Create the ./tmp directory
    status = system("mkdir -p images");

	// Handle CTRL-C to interrupt the optimisation but not lenstool
	signal(SIGINT, signalReset);
	optInterrupt = 0;

	// Loop over each line
	for( iVal = 0; iVal < nVal && !optInterrupt; iVal++ )
	{
		// Set the lens parameters from <array>
		setBayesModel( iVal, nVal, array );
		printf( "INFO: Compute the image.all file %ld/%ld  [CTRL-C to interrupt]\r", iVal+1, nVal);
		fflush(stdout);

        if( mode == 0 )
           ni = bc_mode();
        else
           ni = image_mode();

		// write image array
		sprintf( fname, "images/image%ld.all", iVal );
		ecrire_r(0,ni,cimage,fname,2);
	}
	printf( "\n" );
	free( array );
	return 0;
}

static int bc_mode()
{
	struct point Bs;
	struct point Ps[NIMAX];
    int i, j;
    int ni;

	for( i = 0 ; i < I.n_mult ; i++ )
	{
		S.ns = i;

		// arclet --> source list
		e_unlens(I.mult[i],multi[i], &S.ns, source);

		for( j = 0; j < I.mult[i]; j++ )
		{
			multi[i][j].grad.x=multi[i][j].grad.y=0.;
			multi[i][j].grad2.a=multi[i][j].grad2.c=0.;
		}
		
		o_dpl(I.mult[i], multi[i], Ps, NULL);

		/*Find the barycenter position of the computed sources*/
		if (I.forme == 4 || I.forme == 6)
			Bs=weight_baryc(Ps,multi[i],I.mult[i] , i); //i==n_familly
		else
			Bs=bcentlist(Ps,I.mult[i]);  /* barycentre */

		source[i].C = Bs;
		source[i].n[strlen(source[i].n)-1] = 0;
	}
	S.ns = I.n_mult;

       
	// source list --> image array
	e_lensing(source,image);

	// Count the number of images
	ni = 0;
	for( i = 0 ; i < S.ns ; i++ )
		for( j = 0 ; j < NIMAX && strcmp(image[i][j].n,"") ; j++ )
			cimage[ni++]=image[i][j];

    return(ni);
}

static int image_mode()
{
    long int ni, ncistart;

    s_source();
    e_lensing(source, image);
    classer(image, cimage, &ni, &ncistart);
    return(ni);
}


static void signalReset()
{
	signal(SIGINT, SIG_DFL);
	optInterrupt = 1;
	printf( "\nINFO: Optimisation interrupted by CTRL-C\r");
}
