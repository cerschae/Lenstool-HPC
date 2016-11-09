/*
* Dimension definition
*/

#ifndef DIMENSION_H
#define DIMENSION_H

#define NGGMAX 	128   /* maximum grid points in the I->S mapping */
#define NAMAX 	30000 /* maximum number of arclets */
#define NASMAX 	1000  /* maximum number of arclets for study*/
#define IDSIZE	10    /* size in characters of the id of clumps and images*/
#define ZMBOUND   10    /* maximum number of redshift bounded families*/
#define NGMAX 	400	// maximum number of point with the grille command 
#define NMAX 	3000  // maximum number of segments for the critical lines 
#define NPMAX 	5000
#define NPZMAX	9	/* maximum number of critical lines in g_cline struct*/
#define NLMAX 	2000   // maximum number of clumps in the lens[] array
#define NIMAX 	50    /* maximum images per family */
#define NFMAX 	1200    /* maximum number of families */
#define NPAMAX 	35     // number of free parameters (see #define in structure.h)
#define NTMAX 	1024
#define NPOINT 	1024  /* Number of contour points in cleanlens mode*/
#define NPARMAX 50	
#define NMCMAX  500	
#define NSRCFIT 30    // Number of points for source plane fitting
#define NPOTFILE 6    // Maximum number of potfiles

#define ARRAY_SIZE 2000000

/* zero pour les calculs de dichotomie, amplification infinie, pente nulle */
#define PREC_ZERO 	.00001
/* erreur sur dlsds pour le calcul inverse source->image */
#define PREC_DLSDS 	.00001
/* nombre maximal de points sur une ligne critique tangentielle ou radiale*/
#define NTLINEMAX     250
#define NRLINEMAX     250
#define DMIN	1e-4	// distance minimale de convergence dans le plan image (in arcsec)	
#define NITMAX 	100

#define IDPARAM1 1   // column of the 1st physical parameter in array (cf readBayesModel.c)
#define LINESIZE 16000  // size of a line in bayes.dat  
#define FILENAME_SIZE  50  // size of a filename in .par file
 
#endif // if DIMENSION_H
