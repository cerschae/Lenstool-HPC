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

/******************************************************************************
 * Name: bayesGrad.c
 * Authors: EJ
 * Date: 03/01/08
 *
 * Analyse a bayes.dat file and for each model, compute the deflection angle
 * due to each lens for each image. Return the largest deflection produced
 * (generally for the closest image) in a bayesGrad.dat file, to be analysed
 * with the bayesResults.pl script.
 *
 * syntax : bayesGrad <.par>
 *
 * The bayes.dat file in the current directory is used.
 * Output : bayesGrad.dat
 *
 ******************************************************************************/

typedef void (*sighandler_t)(int);
static void signalReset();
int optInterrupt;

void help_msg()
{
	fprintf(stderr, "Syntax : bayesGrad [OPTION] <.par>\n");
    fprintf(stderr, "Available OPTIONS:\n");
    fprintf(stderr, " -i <potential id> : compute gradients due to a specific potential\n");
    exit(-1);
}

/* Print a line of bayesGrad.dat
 */
void printLine(FILE *bayes, long int ilens)
{
    int i, j;
    long int  k;
    double dx, dy;
    struct point grad;

    for ( i = 0 ; i < I.n_mult; i++ )
        for ( j = 0 ; j < I.mult[i] ; j++ )
        {
            dx = dy = 0.;
            if( ilens != -1 )
            {
                grad = e_grad_pot(&multi[i][j].C, ilens );
                dx = multi[i][j].dr * grad.x; 
                dy = multi[i][j].dr * grad.y; 
            }
            else
                for (k = 0; k < G.nlens; k++ )
                {
                    grad = e_grad_pot(&multi[i][j].C, k );
                    dx += multi[i][j].dr * grad.x; 
                    dy += multi[i][j].dr * grad.y; 
                }
            fprintf( bayes, "%lf ", sqrt(dx * dx + dy * dy));
        }
    fprintf( bayes, "\n" );
}

int main( int argc, char** argv )
{
    extern struct pot lens[];
    double **array; // contains the bayes.dat data
    int nParam, j;  // size of array
    long int    nVal, iVal, i;
    FILE   *bayes;
    double *index;  // list of bayes.dat lines
    int    seed;   // random seed
    int tmp;

    // Check the arguments
    if ( argc < 2 )
        help_msg();

    if ( strstr(argv[1], ".par") == NULL )
        help_msg;

    char id_lens[IDSIZE];
    id_lens[0] = 0;   // default empty
    if( !strcmp(argv[1], "-i" ) )
    {
        strcpy(id_lens, argv[2]);    
		for( i = 1; i < argc-1 ; i++ )
			argv[i]=argv[i+2];
        argc-=2;
    }

    // Read the .par file
    init_grille( argv[1], 1);

    // look for lens id in lens[] array
    long int ilens = -1;
    if( id_lens[0] != 0 )
    {
        i = 0;
        while( i < G.nlens && strcmp(id_lens, lens[i].n))  i++;
        ilens = i;
    }

    // Read constraints
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

    // Write the header of the bayesGrad.dat file
    bayes = fopen( "bayesGrad.dat", "w" );

    fprintf( bayes, "#Nsample\n");
    fprintf( bayes, "#Chi2\n");

    for ( i = 0 ; i < I.n_mult ; i++ )
        for ( j = 0; j < I.mult[i]; j++ )
            fprintf( bayes, "#Image %s\n", multi[i][j].n);

    // Prepare the index list
    index = (double *) malloc((unsigned) nVal * sizeof(double));
    for ( i = 0 ; i < nVal ; i++ ) index[i] = i;
    seed = -2;

    // Handle CTRL-C to interrupt the optimisation but not lenstool
    signal(SIGINT, signalReset);
    optInterrupt = 0;

    // Loop over each line
    for ( i = 0; i < nVal && !optInterrupt; i++ )
    {
        // Randomly draw a line from index array
//        tmp = (int) floor(d_random(&seed) * (nVal - i));
//        iVal = index[i + tmp];
        // and swap the indexes in the index list
//        index[tmp] = i;
        iVal = i;

        // Set the lens parameters from <array>
        setBayesModel( iVal, nVal, array );
        printf( "INFO: Compute the bayesGrad.dat file %ld/%ld  [CTRL-C to interrupt]\r", i + 1, nVal);
        fflush(stdout);

        fprintf( bayes, "%ld %lf ", iVal, array[0][iVal] );
        printLine(bayes, ilens);
        fflush(bayes);
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
