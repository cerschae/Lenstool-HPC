
#include <setup.hpp>

#ifdef __WITH_LENSTOOL
void setup_lenstool(){
/*
	extern  struct  g_mode    M;

    extern  struct  g_grille    G;
    extern  struct  g_frame     F;
    extern  struct g_image   I;
    const extern struct pot       lens[NLMAX];

    ////////Initialisation/////////
    M.iclean = 0;
    //nbgridcells
    G.ngrid = 1000;
    //Frame
    F.xmin = -50.;
    F.xmax = 50.;
    F.ymin = -50.;
    F.ymax = 50.;
    //Image Control Structure
    I.shmap = 0;
    I.stat = 0;
    I.n_mult=1;
    I.srcfit = 0;
    I.npcl = 0;
    I.forme = -1;



    char *multfile = "../ConfigFiles/theoretical_images_time.txt";
    strcpy(I.multfile, multfile);
    */
	char    infile[80];
	strcpy (infile, "../ConfigFiles/90Pot_Lenstool.par");

	extern struct g_mode    M;
	extern struct g_grille  G;
	extern struct galaxie   source[NFMAX];
	extern struct galaxie   image[NFMAX][NIMAX];
	extern struct g_cosmo  C;
	extern struct pot      lens[NLMAX];
	int noedit = 1, init;
	//double  chi2, lhood0;
	//char    infile[80];
	long int ni, ncistart;

	#ifdef _OPENMP
	fprintf(stderr, "You are running openMP version of lenstool with %d threads\n", omp_get_max_threads());
	fprintf(stderr, "You can change number of threads by set environment variable OMP_NUM_THREADS\n");
//    fprintf(stderr, "ATTENTION!!! Make sure that you have at least %d FREE cores, otherwise you will have huge slow down\n", omp_get_max_threads());
	#endif

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
}
#endif
