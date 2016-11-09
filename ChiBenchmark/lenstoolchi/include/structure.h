#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "dimension.h"
#include<gsl/gsl_matrix.h>

// Parameter constants
#define CX       0
#define CY       1
#define EPOT     2
#define EMASS    3
#define THETA    4
#define PHI      5
#define RC       6
#define B0       7
#define ALPHA    8
#define BETA     9
#define RCUT    10 
#define MASSE   11 
#define ZLENS   12
#define RCSLOPE 13
#define PMASS   14
#define OMEGAM  15
#define OMEGAX  16
#define WX      17
#define WA      18
#define SCX     19
#define SCY     20
#define SA      21
#define SB      22
#define SEPS    23
#define STHETA  24
#define SINDEX  25
#define SFLUX   26
#define VFCX    27
#define VFCY    28
#define VFVT    29
#define VFRT    30
#define VFI     31
#define VFTHETA 32
#define VFLCENT 33
#define VFSIGMA 34

/*
* structure definition
*/

/*****************************************************************/
/*                               */
/*      Definition de type               */
/*                               */
/*****************************************************************/

typedef struct
{
	double re;
	double im;
} complex;

/*****************************************************************/
/*                               */
/*      structure  point                 */
/*                               */
/*****************************************************************/

struct point
{
	double x;
	double y;
};

/*****************************************************************/
/*                               */
/*      structure  lens-data                 */
/*                               */
/*****************************************************************/

typedef struct
{
	double alpha_now, x_now, kappa, dpl;
} lensdata;


/*****************************************************************/
/*                               */
/*      Definition des structures de controle        */
/*                               */
/*****************************************************************/

struct g_observ
{
	int     bruit;
	int     setbin;
	int     bin;
	int     setseeing;
	int     filtre;
	double  seeing;
	double  seeing_a;
	double  seeing_b;
	double  seeing_angle;
        char    psffile[FILENAME_SIZE];
	double  r0st;
	double  r0st_a;
	double  r0st_b;
	double  prec;
	double  SKY;
	double  gain;
	int     idum;
};

struct g_mode
{
	int     verbose;
	int     sort;
	int     grille;
	int     ngrille;
	double  zgrille;
	int     inverse;    /* mode inversion de la lentille */
	int     itmax;      /* iter. max lors de l'optimisation*/
	double  rate;       /* rate in the bayesian optimisation*/
	double  minchi0;        /* minimal chi2 for the optimisation*/
	int     ichi2;          /* compute the chi2 with the given model and constraints*/
	int     image;      /* flag pour fichier d'images */
	char    imafile[FILENAME_SIZE];
	int     source;     /* flag pour fichier de sources */
	char    sourfile[FILENAME_SIZE];
	/* from an image file compute only the source catalogue */
	int     sof;        /* flag pour fichier de sources */
	char    imsfile[FILENAME_SIZE];
	char    sfile[FILENAME_SIZE];
	/* correlate shear with model */
	int     icorshear;  /* flag pour fichier de correlation */
	char    corshfile[FILENAME_SIZE];
	/* etude pour le comportement local des image */
	int     local;  /* flag d'etude de e_i dans le repere local */
	char    localfile[FILENAME_SIZE];
	/* etude pour l'inversion des arclets en fonction du redshift */
	int     study;  /* flag d'etude de e_s en fonction de z_s */
	char    studyfile[FILENAME_SIZE];
	int     fake;
	int     mean;
	double  seeing;
	/* etude pour le calcul du amplification du plan image */
	int     iampli;
	int     nampli;
	double  zampli;
	char    amplifile[FILENAME_SIZE];
	/* etude pour le calcul du poten du plan image */
	int     ipoten;
	int     npoten;
	double  zpoten;
	char    potenfile[FILENAME_SIZE];
	/* etude pour le calcul du mass du plan image */
	int     imass;
	int     nmass;
	double  zmass;
	char    massfile[FILENAME_SIZE];
	/* etude pour le calcul du dpl du plan image */
	int     idpl;
	int     ndpl;
	double  zdpl;
	char    dplxfile[FILENAME_SIZE];
	char    dplyfile[FILENAME_SIZE];
	/* etude pour le calcul du curvature du plan image */
	int     icurv;
	int     ncurv;
	double  zcurv;
	char    cxxfile[FILENAME_SIZE];
	char    cxyfile[FILENAME_SIZE];
	char    cyyfile[FILENAME_SIZE];
	/* etude pour le calcul du shear-field du plan image */
	int     ishearf;
	double  zshearf;
	char    shearffile[FILENAME_SIZE];
	int     nshearf;
	/* etude pour le calcul du amplification-field du plan image */
	int     iamplif;
	double  zamplif;
	char    ampliffile[FILENAME_SIZE];
	/* etude pour le calcul du shear du plan image */
	int     ishear;
	int     nshear;
	double  zshear;
	char    shearfile[FILENAME_SIZE];
	/* etude pour le calcul du time_delay du plan image */
	int     itime;
	int     ntime;
	double  ztime;
	char    timefile[FILENAME_SIZE];
	/* etude pour le calcul des props du plan image */
	int     prop;
	int     nprop;
	double  zprop;
	char    propfile[FILENAME_SIZE];
	double  radius;
	double  masse;
	/* visualisation pixelise du champ */
	int     pixel;
	int     npixel;
	char    pixelfile[FILENAME_SIZE];
	/* datacube */
	int 	cube;
	int 	nslices;
	char 	cubefile[FILENAME_SIZE];
	/* Reference absolue */
	int     iref;
	double  ref_ra;
	double  ref_dec;
	/* calcul des markers sources de markers images donnes */
	int     marker;
	double  zmarker;
	char    markfile[FILENAME_SIZE];
	/* etude d'une coupe des props du plan image */
	int     radial;
	double  zradial;
	double  theta;
	/* visualisation pixelise du plan source */
	int     iclean;
	double  zclean;
	/* center of multiple images file */
	char    centerfile[FILENAME_SIZE];
};

/* potfile parameters */
struct  g_pot
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
};
/* dynfile parameters */
struct	g_dyn	{
	int		    dyntype;
	int	        dynnumber;
    double      dynvel;
    double      dynevel;
	double      indmass;
	double      indemass;
	double      refradius;
};

// parameters of an image in pixels with WCS coordinates
struct g_pixel
{
	int     column;
	int     ech;
	int     format;
	char    pixfile[FILENAME_SIZE];
	int     ncont;
	char    outfile[FILENAME_SIZE];
	char    contfile[10][FILENAME_SIZE];
	double  pixelx;
	double  pixely;
	int     header;
	double  xmin;
	double  xmax;
	double  ymin;
	double  ymax;
	int     nx;
	int     ny;
        double  meanFlux;
	struct  WorldCoor *wcsinfo;
};

struct g_cube
{
	int     format;
	char    pixfile[FILENAME_SIZE];
	double  xmin;
	double  xmax;
	double  ymin;
	double  ymax;
	double  lmin;
	double  lmax;
	double  pixelx;
	double  pixely;
	double  pixelz;
	int     nx;
	int     ny;
	int     nz;
        int header;
        double  meanFlux;
	struct  WorldCoor *wcsinfo;
};

struct g_image
{
	int     random;
	int     nzlim;
	int     npcl;
	/* shear map parameteres */
	int     shmap;
	double  zsh, dl0ssh, dossh, drsh;
	char    shfile[FILENAME_SIZE];
	/* arclets parameteres */
	int     stat;
	int     statmode;
	char    arclet[FILENAME_SIZE];
	char    nza[10];    /*index of the redshift known arclet in z_arclet*/
	double  zarclet;
	double  drarclet;
	double  sigell, dsigell; // ellipticity of sources and associated error
	double  sig2ell;
    // source plane fitting
    int     srcfit;    // boolean yes/no
    long int  nsrcfit;   // number of points in srcfit global variable
    char    srcfitFile[FILENAME_SIZE];  // name of the srcfit catalog of points
    char    srcfitMethod[50]; // source plane fitting method
	/* multiple arcs parameteres */
	int     forme;
	int     n_mult;
	int     mult_abs;
	int     mult[NFMAX];
	char    multfile[FILENAME_SIZE];
	double  sig2pos[NFMAX][NIMAX];  // position error per system
    double  **weight;       // inverse of covariance matrix
    double  detCov;      // determinant of the covariance matrix
	double  sig2amp;
	double  Dmag;
	int     adjust;
	char    Afile[FILENAME_SIZE];
	int     Anx;
	int     Any;
	int     Abin;
	int     Anfilt;
	int     Anpixseeing;
	double  Apixel;
	double  Aseeing;
	double  Axmin;
	double  Axmax;
	double  Aymin;
	double  Aymax;
	double  Amean;
	double  Adisp;
};

struct g_source
{
	int     grid;
	int     rand;
	long int  ns;
	int     distz;
	double  emax;
	double  zs;
	double  zsmin;
	double  zsmax;
	double  par1; // for Smail et al. distrib
	double  par2;
	double  par3;
	double  taille;
	double  lfalpha;
	double  lfm_star;
	double  lfm_min;
	double  lfm_max;
};

struct g_grille
{
	int     ngrid;       // number of cells in source-image plane inversion
	int     pol;         // boolean : 0 regular grid, 1 polar grid
	long int     nlens;       // size of lens[] list
	long int     nplens[NPOTFILE+1];      // array of potfile starting indexes in lens[]
    int     npot;    // number of potfiles
	long int     nlens_crit;  // number of lens for which the critical lines have to be computed
	long int     no_lens;     // individual clump to optimise index in lens[]
	long int     nmsgrid;     // multi-scale grid final index in lens[]
	char    splinefile[FILENAME_SIZE];
	double  xmin;
	double  xmax;
	double  ymin;
	double  ymax;
	double  dx;
	double  dy;
	int     nx;
	int     ny;
	int     echant;
	double  exc;
	double  excmin;
    double **invmat;      // msgrid transformation matrix between vector Sigma and sig2
};

// Multi-scale grid definition
struct g_msgrid
{
	double  threshold;   // splitting threshold
	int     levels;      // number of splitting
	double  param;       // for PIEMD grid : rcut/rcore
	char    gridfile[FILENAME_SIZE];// file containing the description of a grid
};

struct g_cline
{
	int     nplan;
	double  cz[NPZMAX];
	int     zone;
	int     npzone;
	char    zonefile[FILENAME_SIZE];
	double  cpas;   // minimum size of the squares in MARCHINGSQUARES or initial step size in SNAKE
	double  dmax;
	char    algorithm[20];  //SNAKE or MARCHINGSQUARES algorithm for critical lines search
	double  limitHigh; // maximum size of the squares in MARCHINGSQUARES algorithm
};

struct g_frame
{
	double  xmin;
	double  xmax;
	double  ymin;
	double  ymax;
	double  lmin;
	double  lmax;
	double  rmax;
};
struct g_large
{
	double  dlarge;
	char    iname[50];
	int     vitesse;
	int     iso;
	int     nmaxiso;
	double  scale;
	double  zonex;
	double  zoney;
	int     profil;
	int     pt;
	int     ncourbe;
	int     npt;
};

struct g_cosmo
{
   //function g_cosmo_equal make comparison of two g_cosmo structures
   //see distcosmo2_cash.c
    int     model;            //TV
    double  omegaM;
    double  omegaX;
    double  kcourb;
    double  wX;
    double  wa;
    double  H0;
    double  h;
};

/*****************************************************************/
/*                               */
/*              Definition des types             */
/*                               */
/*****************************************************************/

struct polar
{
	double  r;
	double  theta;
};

struct segment
{
	double  dist;
	int         k;
};

struct matrix
{
	double  a;
	double  b;
	double  c;
	double  d;
};

struct vecteur
{
	double  x;
	double  y;
};

struct biline
{
	int         i;
	struct point I;
	struct point S;
};

struct ellipse
{
	double  a;
	double  b;
	double  theta;
};

struct pointgal
{
	int         n;
	struct point C;
	double  z;
};

struct pointell
{
	struct point   C;
	struct ellipse E;
};

struct shear
{
	int n;
	struct point C;
	double  mx;
	double  my;
	double  dx;
	double  dy;
	double  err;
};

//"class" for cashing results of distcosmo2 (see cosmoratio.c )
struct distcosmo2_cash
{
   //you should initialise distcosmo2_cash before using
   //set z1 = 0 z2 = 0 dsl=0 (or run distcosmo2_cash_init)
   //or define it as global variables which will set zero automatically
   //Of course, in c++ life is easy and we simple could write constructor
   struct g_cosmo C;  //cosmologe for which dls = distcosmo2(z1,z2) was calculated 
   double z1, z2;  
   double dls;  //cashing value
};
//"members" of "class" distcosmo2_cash
//you should run distcosmo2_cash_init after creatig distcosmo2_hash
void  distcosmo2_cash_init(struct distcosmo2_cash* cash);
//calculate distcosmo2(z1,z2) and store result in the case
//or get value from the cash if we've already calculated it
double distcosmo2_cash_request(struct distcosmo2_cash* cash, double z1, double z2);
  


struct galaxie
{
	char    n[IDSIZE];
	struct point    C;
	struct point    Grad;   // total deflection with all clumps
	struct ellipse  E;
	char    c;
	int     type;
	double  magabs;
	double  mag;
	double  flux;
	double  mu;
	double  I0;
	double  z;
	double  dl0s;         // DLS distcosmo2 between lens[0] and z
	double  dos;          // DOS distance to z
	double  dr;           // ratio dl0s/dos
	double  q;
	double  eps;
	double  tau;
	double  dis;
	double  A;
	double  ep;
	double  tp;
	double  dp;
	double  thp;
	double  taux;
	double  tauy;
	double  time;
	double  kappa;
	double  gamma1;
	double  gamma2;
	double  var1;           // multipurpose variable1
	double  var2;           // multipurpose variable2
	struct point grad;      // deflection due to the not optimised clumps
	struct matrix grad2;    // absolute projected potential of the not optimised clumps
	struct point *np_grad; // deflection due to non parametric clumps with b0=1
    double *np_grad2a; // laplacian due to non parametric clumps with b0=1
    double *np_grad2b; // laplacian due to non parametric clumps with b0=1
    double *np_grad2c; // laplacian due to non parametric clumps with b0=1    

   struct  distcosmo2_cash dls_cash; //cash of distcosmo2 result

   
   struct  point   (*gsource)[NGGMAX][NGGMAX]; //we use it only for single image families!
                              //(see chi2SglImage function in o_chi.c)
      //gsource is a 2D map of G.ngrid^2 points in the source plane defined
      //by dlsds. Each point in gsource is linked to a single point in the
      //gimage global variable. Those 2 maps define a kind of bijection
      //between the source and the image planes. 
   double grid_dr; //dlsds ratio for the grid used from source to image plane
};

struct gal_pol
{
	int     n;
	struct polar    C;
	struct ellipse  E;
	double  z;
};

struct z_lim
{
	int     opt;
	char    n[128];
	int     bk;
	int     bk0;   // opt_mc stuff golse2002
	double   percent;  // opt_mc stuff golse2002
	double  min;
	double  max;
	double  err;
	double  ddmin;
	double  ddmax;
	double  dderr;
	double  excu;
	double  excd;
};

struct sigposStr
{
	int bk;
	double min;
	double max;
	double prec;
	double excu;
	double excd;
};

struct propertie
{
	int     n;
	struct point P;
	double  k;
	double  s;
	double  e;
	double  t;
	double  q;
	double  theta;
	double  d;
	double  A;
	double  g;
};

struct ligne
{
	double  i;
	double  e;
	double  theta;
};

struct line
{
	double  e;
	double  theta;
	double  i;
	double  phi;
};

struct cline
{
	int n;
	struct  point   C;
	double  phi;
	double  dl;
	double  z;
	double  dos;    // distcosmo1 to redshift z
	double  dl0s;   // distcosmo2 between lens[0] and z
	double  dlsds;  // ratio of dl0s/dos
};

struct pot
{
	int     type;
	char    n[IDSIZE];
	struct point C;
	double  epot;
	double  emass;
	double  theta;
	double  phi;  // triaxiality
	double  sigma;
	double  rc;
	double  rckpc;
	double  rcut;
	double  rcutkpc;
	double  alpha;
	double  beta;
	double  rcslope;
	double  z;
	double  dlsds;  /*ratio D(LS)/D(OS) angular distances*/
	double  psicut;
	double  psimcut;
	double  psiccut;
	double  b0;
	double  masse;
	double  pmass;
	double  cr;
	double  ct;
	double  mag;
	double  lum;
	double  mtol;
    double  effradius;
//      double  omegaM;
//      double  omegaX;
//      double  wX;
};

struct ipot
{
	int     pmax;
	int     map;
	int     map_z;
	int     lens[10];
	int     para[10];
	int     lens_z[10];
	int     para_z[10];
	int     zlim[2];
	int     extend;
	int masse;
};

struct triplet
{
	struct point a;
	struct point b;
	struct point c;
};

struct bitriplet
{
	struct triplet i;
	struct triplet s;
};

struct chaine
{
	struct triplet S;
	struct triplet I;
	struct chaine *F;
};

struct arbre
{
	struct galaxie *N;
	struct arbre *FG;
	struct arbre *FD;
};

struct variation
{
	struct point C;
	double elip;
	double theta;
};

struct pixlist
{
	double  i;
	double  j;
	double flux;
};

struct MCarlo
{
	int     optMC;
	int     n_MonteCarlo;
	int     iterations;
	int     squares_par;
	int     tosses_sq;
};

//velocity field
struct vfield
{
	double vt;
	double rt;
	struct point C;
	double i;
	double theta;
	double lcent;
	double sigma;
	int profile;
};

#endif // if STRUCTURE_H
