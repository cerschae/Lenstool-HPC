/**
* @file   structure.h
* @Author Thomas Jalabert, EPFL (me@example.com) , Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch)
* @date   July 2015
* @version 0,1
* @brief  Header file to define the used structures (e.g. defined structs)
* 
* @param configuration file (parameters.par)
* @return Depends on choice in the configuration file, e.g. least chi2 model
*/


// Header guard
#ifndef STRUCTURE_H
#define STRUCTURE_H


#include <iostream>



/*****************************************************************/
/*                               */
/* Constants: Will be migrated to constants.h when there are too many of them*/
/*                               */
/*****************************************************************/


// GPU definitions
#define threadsPerBlock 512 // threads for each set of images
#define MAXIMPERSOURCE 20 // maximum number of multiple images for one source
#define MAXIM 200 // maximum  total number of images treated

// Dimensions definitions
#define NPZMAX	9	/* maximum number of critical lines in g_cline struct*/
//#define CLINESIZE	500	/* maximum number of critical and caustic points for Cline mode. Adjust depending on RAM*/
#define NPOTFILESIZE 2000 //maximum number of potential in potfiles
#define DMIN	1e-4	// distance minimale de convergence dans le plan image (in arcsec)	
#define NITMAX 	100


// gNFW definitions
#define gNFW_ARRAY_SIZE 1800 // Set the dimension of the gnfw table gNFW_ARRAY_SIZE, must be 1800 for the current table file

// Filename definition
#define FILENAME_SIZE 50  // size of a filename in .par file

//constants definition

#define pi_c2  7.209970e-06	/* pi en arcsecond/ c^2 =  648000/vol/vol */
#define cH2piG  0.11585881	/* cH0/2piG en g/cm^2 (H0=50) */
#define cH4piG  0.057929405	/* cH0/4piG en g/cm^2 (H0=50) */
#define cH0_4piG  2.7730112e-4	/* cH0/4piG en 10^12 M_Sol/kpc^2 (H0=50) */
#define d0	29.068701	/* vol/h0*1000/rta -- c/H0 en h-1.kpc/arcsecond  (h0=50)*/

#define MCRIT12	.2343165	/* c^3/4piGh0/RTA/RTA in 1e12 M_sol/arcsec^2 (h0=50) */

/*****************************************************************/
/*                               */
/* 			Types definition			*/
/*                               */
/*****************************************************************/

/** @brief Point: Structure of 2 coordinates
 * 
 * @param x: X coordinate
 * @param y: Y coordinate
 */
struct point    
{
	double x;	
	double y;
};

/** @brief Complex: Structure of 2 doubles
 * @param re: Real Part
 * @param im: Imaginary Part
 */
struct complex    
{
	double re;	
	double im;
};

/** @brief Segment: Structure of two points
 */
struct segment    
{
	point a;	
	point b;
};

/** @brief triplet: Structure of three points defining a triangle
 */
struct triplet
{
	struct point a;
	struct point b;
	struct point c;
};

/** @brief bitriplet: Defines two linked triangles (one in source plane and one in image plane)
 * @param i: Triangle on image plane
 * @param s: Triangle on source plane 
 */
struct bitriplet
{
	struct triplet i;	
	struct triplet s;
};

/** @brief contains the table needed to compute the potential derivatives of general NFW profile
 */
typedef struct
{
	double alpha_now, x_now, kappa, dpl;
} gNFWdata;

/** @brief Matrix: 2by2 doubles
 */
struct matrix
{
	double  a;
	double  b;
	double  c;
	double  d;
};

/** @brief ellipse: for shape computation
 * @param a: semimajor axis
 * @param b: semiminor axis
 * @param theta: shape ellipticity
 */
struct ellipse
{
	double  a;		
	double  b;		
	double  theta;	
};

/** @brief Storage type for sources, lenses and arclets
 * @param center: position of the center of galaxy
 * @param shape: shape of galaxy
 * @param mag: magnitude
 * @param redshift: redshift
 * @param dls: Distance lens-source
 * @param dos: Distance observer-source
 * @param dr: dls/dos
 */

struct galaxy
{
	//char    name[IDSIZE];
	struct point    center;		
	struct ellipse  shape;		
	double  mag;				
	double  redshift;			
	double  dls;         		
	double  dos;          		
	double  dr;           		
};

/** @brief Contains the information for optimising a parameter in the inverse mode
 * @param block: blockorfree variable (whether a parameter is blocked or free for the mcmc algorithm)
 * @param min: lower optimisation value
 * @param max: upper optimisation value
 * @param sigma: optimisation step (MIGHT NOT BE TAKEN INTO ACCOUNT)
 */
struct optimize_block
{
	int block; 			
	double min;			
	double max;			
	double sigma;	
};
/** @brief two optimize_block to simulate a point
 */
struct optimize_point    // blockorfree for the point structure
{
	struct optimize_block x;
	struct optimize_block y;
};

/** @brief Contains the information for optimising the potential in the inverse mode
 * @param position: position of the center of the halo
 * @param weight: weight of the clump (the projected mass sigma0 for PIEMD, the density rhoscale for NFW)
 * @param b0: Impact parameter
 * @param ellipticity_angle: orientation of the clump
 * @param ellipticity: Mass ellipticity
 * @param ellipticity_potential: Potential ellipticity
 * @param rcore: PIEMD specific value
 * @param rcut: PIEMD specific value
 * @param rscale: scale radius for NFW, Einasto
 * @param exponent: exponent for Einasto
 * @param vdisp: Dispersion velocity
 * @param alpha: exponent for general NFW
 * @param einasto_kappacritic: critical kappa for Einasto profile
 * @param z: redshift
 */
struct potentialoptimization  // block or free variable for the MCMC for the potential
{
	struct optimize_point position;					
	struct optimize_block weight;
	struct optimize_block  b0; 						
	struct optimize_block ellipticity_angle;
	struct optimize_block ellipticity;				
	struct optimize_block ellipticity_potential; 	
	struct optimize_block rcore;					
	struct optimize_block rcut;						
	struct optimize_block rscale;
	struct optimize_block exponent;
	struct optimize_block vdisp;					
	struct optimize_block alpha;
	struct optimize_block einasto_kappacritic;
	struct optimize_block z;						
	

};

/** @brief Contains the information of the lens potential
 * @param type: 1=PIEMD , 2=NFW, 3=SIES, 4=point mass, 5=SIS, 8=PIEMD
 * @param type_name: IEMD, NFW, SIES, point
 * @param name: name of the clump (e.g. name of the galaxy) : not compulsory
 * @param position: position of the center of the halo
 * @param weight: weight of the clump (the projected mass sigma0 for PIEMD, the density rhoscale for NFW)
 * @param b0: Impact parameter
 * @param ellipticity_angle:
 * @param ellipticity: Mass ellipticity
 * @param ellipticity_potential: Potential ellipticity
 * @param rcore: PIEMD specific value
 * @param rcut: PIEMD specific value
 * @param rscale: scale radius for NFW, Einasto
 * @param exponent: exponent for Einasto
 * @param vdisp: Dispersion velocity
 * @param alpha: exponent for general NFW
 * @param einasto_kappacritic: critical kappa for Einasto profile
 * @param z: redshift
 */

struct Potential_SOA
{
	int*    type; // 1=PIEMD ; 2=NFW; 3=SIES, 4=point mass
	char    type_name[10]; // PIEMD, NFW, SIES, point
	char    name[20]; // name of the clump (e.g. name of the galaxy) : not compulsory
	//struct point position; // position of the center of the halo
	double* position_x; // position of the center of the halo
	double* position_y; // position of the center of the halo
	double* weight; // weight of the clump (the projected mass sigma0 for PIEMD, the density rhoscale for NFW)
	double* b0; // Impact parameter
	double*	vdisp;	//Dispersion velocity
	double* ellipticity_angle; // orientation of the clump
	double* ellipticity; // ellipticity of the mass distribition
	double* ellipticity_potential; //ellipticity of the potential
	double* rcore;  // core radius
	double* rcut; // cut radius
	double* rscale; // scale radius for NFW, Einasto
	double*	exponent; // exponent for Einasto
	double* alpha; // exponent for general NFW
	double* einasto_kappacritic; // critical kappa for Einasto profile
	double* z; // redshift of the clump
	double* mag; //magnitude
	double* lum; //luminosity
	double* theta; //theta
	double* sigma; // sigma
};

struct Potential
{
       int    type; // 1=PIEMD ; 2=NFW; 3=SIES, 4=point mass
       char    type_name[10]; // PIEMD, NFW, SIES, point
       char    name[20]; // name of the clump (e.g. name of the galaxy) : not compulsory
       struct  point position; // position of the center of the halo
       double  weight; // weight of the clump (the projected mass sigma0 for PIEMD, the density rhoscale for NFW)
       double  b0; // Impact parameter
       double  vdisp;  //Dispersion velocity
       double  ellipticity_angle; // orientation of the clump
       double  ellipticity; // ellipticity of the mass distribition
       double  ellipticity_potential; //ellipticity of the potential
       double  rcore;  // core radius
       double  rcut; // cut radius
       double  rscale; // scale radius for NFW, Einasto
       double  exponent; // exponent for Einasto
       double  alpha; // exponent for general NFW
       double  einasto_kappacritic; // critical kappa for Einasto profile
       double  z; // redshift of the clump
       double  mag; //magnitude
       double  lum; //luminosity
       double  theta; //theta
       double  sigma; // sigma
};

/*****************************************************************/
/*                               */
/*      Control structure definition        */
/*                               */
/*****************************************************************/

/** @brief Control structure for runmode parameters
 * 
 * Default values are set in module_readParameters_readRunmode
 *
 * @param nbgridcells: Number of grid cells
 * @param source: flag for simple source to image conversion
 * @param sourfile: file name for source information
 * @param image: flag for simple image to source conversion
 * @param imafile: file name for image information
 * @param mass: flag for mass fitsfile
 * @param massgridcells: Nb of cell for fitsfile
 * @param z_mass: redshift for which to be computed
 * @param z_mass_s: redshift of source for which effect of mass will be computed
 * @param potential: flag for potential fitsfile
 * @param potgridcells: Nb of cell for fitsfile
 * @param z_pot: redshift for which to be computed
 * @param dpl: flag for displacement fitsfile
 * @param dplgridcells: Nb of cell for fitsfile
 * @param z_dpl: redshift for which to be computed
 * @param inverse: flag for inversion mode (MCMC etc,)
 * @param arclet: flag for arclet mode 
 * @param debug: flag for debug mode
 * @param nimagestotal: total number of lensed images in file
 * @param nsets: number of sources attributed to the lensed images
 * @param nhalo: Number of halos
 * @param grid: 0 for automatic grid (not implemented), 1 for grid infor applying on source plane, 2 for grid information applying on image plane
 * @param nbgridcells: Number of grid cells
 * @param zgrid: redshift of grid
 * @param cline: flag for critical and caustic line calculation
 */
 
struct runmode_param
{
	int 	nbgridcells;
	//Source Mode
	int     source;     
	std::string    sourfile;
	int 	nsets;
	//Image Mode
	int     image;
	std::string    imagefile;
	int		nimagestot;
	//Mass Mode
	int		mass;
	int		mass_gridcells;
	double	z_mass;
	double	z_mass_s;
	//Potential Mode
	int		potential;
	int		pot_gridcells;
	double	z_pot;
	int 	nhalos;
	//Potfile Mode
	int		potfile;
	int		npotfile;
	std::string		potfilename;
	//displacement Mode
	int		dpl;
	int		dpl_gridcells;
	double	z_dpl;
	//Inverse Mode
	int     inverse; 
	//Arclet Mode
	int     arclet; 
	//Debug Mode      
	int 	debug;	
	//Grid Mode
	int 	grid;
	int 	gridcells;
	double 	zgrid;
	//Critic and caustic mode
	int		cline;
	//Amplification Mode
	int 	amplif;
	int 	amplif_gridcells;
	double 	z_amplif;
	//Time/Benchmark mode
	int		time;
};

/** @brief Not used yet
 * 
 */
struct image_param
{

};

/** @brief Not used yet
 * 
 */
struct source_param
{

};

/** @brief Contains Grid information
 */

struct grid_param
{
	double  xmin;
	double  xmax;
	double  ymin;
	double  ymax;
	double  lmin;
	double  lmax;
	double  rmax;
};

/** @brief Control structure for cosmological parameters
 * 
 * @param model: Cosmological model
 * @param omegaM: 
 * @param omegaX: 
 * @param curvature: curvature parameter
 * @param wX: 
 * @param wa: 
 * @param H0: Hubble constant
 * @param h: H0/50
 */

struct cosmo_param  
{
    int     model;            	
    double  omegaM;				
    double  omegaX;				
    double  curvature;			
    double  wX;					
    double  wa;					
    double  H0;					
    double  h;					
};

/** @brief Control structure for potfile parameters
 *
 *  @param     	potid:  1: pot P, 2: pot Q
	@param     	ftype:
	@param  	potfile[FILENAME_SIZE];
	@param     	type;
	@param  	zlens;
	@param  	core;CCC
	@param  	corekpc;
	@param  	mag0;
	@param     	select;
	@param     	ircut;
	@param  	cut, cut1, cut2;
	@param  	cutkpc1, cutkpc2;
	@param     	isigma;
	@param  	sigma, sigma1, sigma2;
	@param     	islope;
	@param  	slope, slope1, slope2;
	@param     	ivdslope;
	@param  	vdslope, vdslope1, vdslope2;
	@param     	ivdscat;
	@param  	vdscat, vdscat1, vdscat2;
	@param     	ircutscat;
	@param  	rcutscat, rcutscat1, rcutscat2;
	@param     	ia;   // scaling relation of msm200
	@param  	a, a1, a2;
	@param     	ib;   // scaling relation of msm200
	@param  	b, b1, b2;

 */

struct  potfile_param
{
    int     potid;   // 1: pot P, 2: pot Q
	int     ftype;
	char    potfile[FILENAME_SIZE];
	int     type;
	double  zlens;
	double  core;
	double  corekpc;
	double  mag0;
	int     select;
	int     ircut;
	double  cut, cut1, cut2;
	double  cutkpc1, cutkpc2;
	int     isigma;
	double  sigma, sigma1, sigma2;
	int     islope;
	double  slope, slope1, slope2;
	int     ivdslope;
	double  vdslope, vdslope1, vdslope2;
	int     ivdscat;
	double  vdscat, vdscat1, vdscat2;
	int     ircutscat;
	double  rcutscat, rcutscat1, rcutscat2;
	int     ia;   // scaling relation of msm200
	double  a, a1, a2;
	int     ib;   // scaling relation of msm200
	double  b, b1, b2;
	int npotfile;
};

/** @brief Control structure for caustic and critical lines
 * 
 * @param nplan: number of sourceplanes for which the caustic lines will be computed
 * @param cz: redshift values array for respective sourceplanes
 * @param dos: distcosmo1 to redshift z
 * @param dls: distcosmo2 between lens[0] and z
 * @param dlsds:    ratio of dl0s/dos
 * @param limitLow: minimum size of the squares in MARCHINGSQUARES
 * @param dmax: Size of search area
 * @param xmin: 
 * @param xmax: 
 * @param ymin:    
 * @param ymax: 
 * @param limithigh: maximum size of the squares in MARCHINGSQUARES algorithm
 * @param nbgridcells: nbgridcells * nbgridcells = number of pixels for critical line computation
*/

struct cline_param
{
	int     nplan;
	double  cz[NPZMAX];
	double  dos[NPZMAX];	// distcosmo1 to redshift z
	double  dls[NPZMAX];	// distcosmo2 between lens[0] and z
	double  dlsds[NPZMAX];	// ratio of dl0s/dos
	double  limitLow;   // minimum size of the squares in MARCHINGSQUARES or initial step size in SNAKE
	double  dmax;
	double  xmin;
	double  xmax;
	double  ymin;
	double  ymax;
	double  limitHigh; // maximum size of the squares in MARCHINGSQUARES algorithm
	int nbgridcells; // nbgridcells * nbgridcells = number of pixels for critical line computation
};

#endif // if STRUCTURE_H
