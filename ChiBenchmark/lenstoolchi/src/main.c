#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<omp.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        LENSTOOL            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

/*Global variables declaration*/
struct g_mode   M;
struct g_pot    P[NPOTFILE];  
struct g_pixel  imFrame, wFrame, ps, PSF;
struct g_cube   cubeFrame;
struct g_dyn	Dy;      //   //TV


struct g_source S;
struct g_image  I;
struct g_grille G;
struct g_msgrid H;  // multi-scale grid
struct g_frame  F;
struct g_large  L;
struct g_cosmo  C;
struct g_cline  CL;
struct g_observ O;
struct pot      lens[NLMAX];
struct pot      lmin[NLMAX], lmax[NLMAX], prec[NLMAX];
struct g_cosmo  clmin, clmax;       /*cosmological limits*/
struct galaxie  smin[NFMAX], smax[NFMAX];       // limits on source parameters
struct ipot     ip;
struct MCarlo   mc;
struct vfield   vf;
struct vfield   vfmin,vfmax; // limits on velocity field parameters
struct cline    cl[NIMAX];
lensdata *lens_table;

int  block[NLMAX][NPAMAX];      /*switch for the lens optimisation*/
int  cblock[NPAMAX];                /*switch for the cosmological optimisation*/
int  sblock[NFMAX][NPAMAX];                /*switch for the source parameters*/
int  vfblock[NPAMAX];                /*switch for the velocity field parameters*/
double excu[NLMAX][NPAMAX];
double excd[NLMAX][NPAMAX];

/* suppléments tableaux de valeurs pour fonctions g pour Einasto
 * Ce sont trois variables globales qu'on pourra utiliser dans toutes les fonctions du projet
*/

#define CMAX 20
#define LMAX 80

float Tab1[LMAX][CMAX];
float Tab2[LMAX][CMAX];
float Tab3[LMAX][CMAX];
 

int      nrline, ntline, flagr, flagt;
long int  narclet;

struct point    gimage[NGGMAX][NGGMAX], gsource_global[NGGMAX][NGGMAX]; 
struct biline   radial[NMAX], tangent[NMAX];
struct galaxie  arclet[NAMAX], source[NFMAX], image[NFMAX][NIMAX];
struct galaxie  cimage[NFMAX];
struct pointgal     gianti[NPMAX][NIMAX];

struct point    SC;
double elix;
double alpha_e;

double *v_xx;
double *v_yy;
double **map_p;
double **tmp_p;
double **map_axx;
double **map_ayy;

double myseconds()
{
        struct timeval  tp;
        struct timezone tzp;
        //
        int i = gettimeofday(&tp,&tzp);
        //
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

int  main(int argc, char *argv[])
{

    /*************  declaration de common et locale ****************/

    extern struct g_mode    M;
    extern struct g_grille  G;
    extern struct galaxie   source[NFMAX];
    extern struct galaxie   image[NFMAX][NIMAX];
    int noedit, init;
    double  chi2, lhood0;
    char    infile[80];
    long int ni, ncistart;

    #ifdef _OPENMP
    fprintf(stderr, "You are running openMP version of lenstool with %d threads\n", omp_get_max_threads());
    fprintf(stderr, "You can change number of threads by set environment variable OMP_NUM_THREADS\n");
//    fprintf(stderr, "ATTENTION!!! Make sure that you have at least %d FREE cores, otherwise you will have huge slow down\n", omp_get_max_threads());
    #endif 
   
    /*************  Verification du format de la commande  ****************/
    if ( argc != 2 && argc != 3 )
    {
        fprintf(stderr, "\nUnexpected number of arguments\n");
        fprintf(stderr, "\nUSAGE:\n");
        fprintf(stderr, "lenstool  input_file [-n]\n\n");
        exit(-1);
    }
    if ( argc == 2 || argc == 3 )
        strcpy (infile, argv[1]);

    if ( strstr( infile, ".par" ) == NULL )
        strcat( infile, ".par" );

    noedit = 0;
    if ( argc == 3 && !strcmp(argv[2], "-n") )
        noedit = 1;

    /************** Read the .par file and initialise potentials ***********/
    init = init_grille(infile, noedit);

    if ( noedit == 0 && M.verbose > 0 )
        copyright();


    /************** Print common information *******************************/
    if ( M.iref != 0 )
        NPRINTF(stderr, "Reference coordinates WCS  : %lf, %lf\n",
                M.ref_ra, M.ref_dec);
    NPRINTF(stderr, "Conversion Factor @ z = %.3lf, 1 arcsec == %.3lf kpc\n",
            lens[0].z, d0 / C.h*distcosmo1(lens[0].z));

    /*************** Creation d'un grille polaire ou cartesienne ***********/
    if ( G.pol != 0 )
        gridp();
    else
        grid();

    /*************** Optimisation globale **********************************/
    if ( M.inverse != 0 )
        o_global();

    tracepot();

    if ( M.ichi2 != 0 )
    {
        int error;
        NPRINTF(stdout, "INFO: compute chires.dat\n");
        readConstraints();
        o_chires("chires.dat");
        
        double t1;
	t1 = -myseconds();
	
        error = o_chi_lhood0(&chi2, &lhood0, NULL);
        	
        t1 += myseconds();
	printf("o_chi_lhood0 Time taken %f seconds \n", t1);
        
        printf("INFO: chi2 %lf  Lhood %lf\n", chi2, -0.5 * ( chi2 + lhood0 )  );
        o_global_free();
    }

    return 0;
}