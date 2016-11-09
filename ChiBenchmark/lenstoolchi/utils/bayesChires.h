
/*Global variables declaration*/
struct g_mode	M;
struct g_pot	P[NPOTFILE];
struct g_pixel	imFrame,wFrame,ps,PSF;
struct g_cube   cubeFrame;
struct g_dyn	Dy;      //   //TV
struct g_source	S;
struct g_image	I;
struct g_grille	G;
struct g_frame	F;
struct g_msgrid H;
struct g_large	L;
struct g_cosmo	C;
struct g_cline	CL;
struct g_observ	O;
struct pot		lens[NLMAX];
struct pot		lmin[NLMAX],lmax[NLMAX],prec[NLMAX];
struct pot		clmin,clmax;		/*cosmological limits*/
struct galaxie  smin[NFMAX], smax[NFMAX];       // limits on source parameters
struct ipot		ip;
struct MCarlo	mc;
struct cline	cl[NIMAX];
struct galaxie *srcfit;  // points for source plane fitting
struct vfield vf;
struct vfield vfmin,vfmax;  // limits on velocity field parameters

#define CMAX 20
#define LMAX 80

float Tab1[LMAX][CMAX];
float Tab2[LMAX][CMAX];
float Tab3[LMAX][CMAX];


lensdata *lens_table;

int	 block[NLMAX][NPAMAX];		/*switch for the lens optimisation*/
int	 cblock[NPAMAX];				/*switch for the cosmological optimisation*/
int  sblock[NFMAX][NPAMAX];                /*switch for the source parameters*/
int   vfblock[NPAMAX];                     /*switch for the velocity field parameters*/
double excu[NLMAX][NPAMAX];
double excd[NLMAX][NPAMAX];

int 	 it;
int 	 nrline,ntline,flagr,flagt;
long int	 narclet;
struct point 	gimage[NGGMAX][NGGMAX],gsource[NGGMAX][NGGMAX];
struct biline	radial[NMAX],tangent[NMAX];
struct galaxie 	arclet[NAMAX],source[NAMAX],image[NFMAX][NIMAX];
struct galaxie 	cimage[NAMAX];
struct pointgal 	gianti[NPMAX][NIMAX];

struct point   	SC;
double elix;
double grid_dr;	// dlsds ratio for the grid used from source to image plane
double alpha_e;

double *v_xx;
double *v_yy;
double **map_p;
double **tmp_p;
double **map_axx;
double **map_ayy;

// Global variables defined in o_global.c
struct shear   shm[NAMAX];
struct galaxie multi[NFMAX][NIMAX];
struct z_lim   zlim[NIMAX];
struct z_lim   zalim;
struct z_lim   zlim_s[NIMAX];
struct pot     lmin_s[NLMAX],lmax_s[NLMAX],prec_s[NLMAX];
struct sigposStr sigposAs_s, sigposAs;
struct matrix  amplifi_mat[NIMAX][NIMAX],amplifi_matinv[NIMAX][NIMAX];
struct pixlist plo[100000]; /*list of points that compose the polygon around the image to inverse*/

double z_dlsds;
double chi_im,chip,chix,chiy,chia,chil,chis,chi_vel,chi_mass;   //I added chi_vel and chi_mass TV
double amplifi[NIMAX][NIMAX];
double **imo,**soo,**ero;
double drima;	/* distance ratio between the image to inverse and the main lens*/
double **sm_p;		/* optimization spline map : initial save */

int    nshmap; 
int    optim_z, n_famille;
int    block_s[NLMAX][NPAMAX];
int	 **imuo;
int    nplo;	/*number of points of the polygon around the image to inverse*/

double z_dlsds;
double chi_im,chip,chix,chiy,chia,chil,chis;
double amplifi[NIMAX][NIMAX];
double **imo,**soo,**ero;
double drima;	/* distance ratio between the image to inverse and the main lens*/
double **sm_p;		/* optimization spline map : initial save */

int    nshmap; 
int    optim_z, n_famille;
int    block_s[NLMAX][NPAMAX];
int	 **imuo;
int    nplo;	/*number of points of the polygon around the image to inverse*/

double distmin[NIMAX]; /* List of minimal euclidean distances between 
			  the images of a same familly i*/

int    nwarn;	// warning emitted during image plane chi2

double *np_b0;  // non parametric accelerator
