#ifndef FONCTION_H
#define FONCTION_H

#include <stdio.h>
#include "wcs.h"
#include "errors.h"
#include <gsl/gsl_randist.h>
#include <time.h>
#include <sys/time.h>

/*
*  fonction.h
*  kneib jean-paul
*  June 97
*  OMP Toulouse
*/

/*
*  useful function definition
*/
#define Min(A,B)    ((A)<(B)?(A):(B))
#define Max(A,B)    ((A)>(B)?(A):(B))
#define cube(A)     A*A*A

#define FPRINTF         if (M.verbose > 1) fprintf  /* Full verbose*/
#define NPRINTF         if (M.verbose == 1) fprintf /* Normal verbose */

#define CHECK_THIRD(A)   if(strlen(third) > A) { fprintf(stderr, "[ERROR] %s\n argument too long. Maximum argument size %d\n", third, A); exit(1); }

/*
 * structures in lenstool
 */

#include "structure.h"

/*
 * functions in lenstool
 */

//#ifdef __cplusplus
//"C" {
//#endif
double myseconds();
complex acpx(complex c1, complex c2);
complex acpxflt(complex c1, double f2);
void      add_pm(double **map, int nx, int ny, double x0, double y0, double m);
//void      amplif();
void      amplif_mat();
void      amplif_matinv();
int       bayesHeader();
void      bicubic_coef(double *y, double *y1, double *y2, double *y12, double d1, double d2, double **c);
void      e_im_prec(const struct bitriplet *E, const struct point *P, double dlsds, int *it, struct bitriplet *A);
void      e_inthere(const struct bitriplet *E, const struct bitriplet *I, const struct point *P, double dlsds, struct bitriplet *res);
double  brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
double  chi_invim(double **im, struct pixlist pl[], int npl, double dlsds, double **so, double **er, int **imu);
int   chsigne(struct point A, struct point B, double dl0s, double dos, double zs);
complex ci05(double x, double y, double eps, double rc);
complex ci10(double x, double y, double eps, double rc, double b0);
complex ci15(double x, double y, double eps, double rc, double b0);
double  chiz(double z);
void      classer(struct galaxie image[NFMAX][NIMAX], struct galaxie *cimage, long int *ni, long int *ncistart);
void      cleanlens(double zimage);
int   comparer_pos(struct galaxie *A, struct galaxie *B);
int   comparer_tau(struct galaxie *A, struct galaxie *B);
int   comparer_z(struct galaxie *A, struct galaxie *B);
int   comp_asc(double A, double B);
double  comp_chi_osv(double *vect);
double  comp_chi_osv(double *vect);
void      comp_dchi_osv(double *vect, double *dvect);
int   comp_desc(double A, double B);
void convertXY( double *x, double *y, int iref, double ra, double dec );
void      copyright();
void      copyTriplet(struct triplet *in, struct triplet *out);
void      cor_seeing(long int n, struct galaxie *image, double seeing);
void      cor_shear();
complex coscpx(complex c);
complex coshcpx(complex c);
void      cp_diffim(double **im1, double **im2, int ni, int nj, double **resim);
double  cp_errdiff(double **im1, double **im2, int ni, int nj, double sigim);
void      cp_im(double **im, int nx, int ny, double xmin, double xmax, double ymin, double ymax, struct galaxie *source, int nbs);
complex cpx(double re, double im);
void      crea_filtre(double seeing, double scale, double **filtre, int n);
void      critic_an();
void      criticinv(double xmin, double ymin, double xmax, double ymax);
void      criticnew(int verbose);
complex csiemd(double x, double y, double eps, double b0);
void      cv_cpsf(double **im, int nx, int ny, double xmin, double xmax, double ymin, double ymax, double seeing);
void      d_bruiter(double **z, int nx, int ny);
complex dcpx(complex c1, complex c2);
complex dcpxflt(complex c, double f);
double  determinant(const struct point *A, const struct point *B, const struct point *C);
double  d_gauss(double sig, int *idum);
struct ellipse diag(double a, double b, double c);
double  diff_mag(struct galaxie *arclet, struct point *guess);
double  d_integrer(struct galaxie A, double x, double y, double t, double f, int n);
double  dist2(struct point A, struct point B);
double  distcosmo1(double z);
double  distcosmo2(double z1, double z2);
void      dist_min();
void      distor(struct galaxie *image, long int ni);
double  distprime1(double z);
double  distprime2(double z1, double z2);
double  dist(struct point A, struct point B);
double  dlumcosmo1(double z);
void      do_itos(double **im, struct pixlist *pl, int npl, double dlsds, double **source, double **erreur, int **imult);
double  d_poisson(double xm, int *idum);
double  d_profil(double x, double y, const struct galaxie *gal);
double  d_random(int *idum);
double  dratio(double zl, double zs);
void dratio_gal(struct galaxie *arclet, double zl);
double  dratioprime(double zl, double zs);
double  d_rndschechter(int *idum);
int       d_rndtype(int *idum);
double  d_rndz(double zmax, int *idum);
gsl_ran_discrete_t* smailpreproc();
double d_rndzsmail(gsl_rng * r, gsl_ran_discrete_t * g);
void      d_seeing(double **im, int nx, int ny, double scale);
double  e_amp(const struct point *position, double dl0s, double dos, double zs);
double  e_amp_gal(struct galaxie *image, double *np_b0);
complex ecpx(complex c);
void      ecrire_r(long int nstart, long int nstop, struct galaxie *liste, char *name, int printShearFlag);
void      e_dpl(const struct point *gi, double dlsds, struct point *gs);
void      e_giant(int *ng, struct galaxie giants, struct galaxie *image );
struct point e_grad(const struct point *pi);
struct point e_grad_gal(struct galaxie *image, double *np_b0);
struct point e_grad_pot(const struct point *pi, long int ilens);
struct matrix e_grad2(const struct point *pi, double dl0s, double zs);
struct matrix e_grad2_gal(struct galaxie *image, double *np_b0);
struct matrix e_grad2_pot(const struct point *pi, long int ilens);
void   e_lensing(struct galaxie source[NFMAX], struct galaxie image[NFMAX][NIMAX]);
int    e_lens_P(struct point ps, struct point pim[NIMAX], double dlsds);
int    e_lens(struct galaxie source, struct galaxie image[NIMAX]);
struct ellipse e_mag(struct point *position, double dl0s, double dos, double zs);
struct ellipse e_mag_gal(struct galaxie *image);
double e_mass(long int icl, double radius);
struct ellipse e_unmag(const struct point *position, double dl0s, double dos, double zs);
struct ellipse e_unmag_gal(struct galaxie *image);
struct point e_zeroamp(struct point A, struct point B, double dl0s, double dos, double zs);
void   e_pixel(int np, char *iname, char *sname, struct galaxie *source);
double e_pot(struct point pi, double dlsds);
double sersic_dpl(double r, double re, double n, double kappae);
double sersic_kappa_eps(double r, double re, double n, double theta, double kappas, double eps);
struct point sersic_gamma_eps(double r, double re, double n, double theta, double kappas, double eps);
double err_invim(double **errim, int **imult);
void   e_tau(long int n, struct galaxie gal[NAMAX]);
void   e_testg(int i, int j, struct bitriplet *Tsol, int ni, struct point *P, double dlsds);
int    e_test_P(struct bitriplet *Tsol, int ni, struct point *ps, struct point image[NIMAX], double dlsds, double err);
int    e_test(struct bitriplet *Tsol, int ni, struct galaxie source, struct galaxie image[NIMAX]);
double e_time(struct point pi, double dlsds);
void   e_unlens_fast(int nima, struct galaxie *strima, struct galaxie *source);
void   e_unlensgrid(struct  point gsource[][NGGMAX], double dlsds);
void   e_unlens(long int na, struct galaxie *arclet, long int *ns, struct galaxie *source);
void      fftcc_im(double **r_im, double **i_im, double **tfr_im, double **tfi_im, int n, int flag);
void      fftc_im(double **im, double **tfr_im, double **tfi_im, int n, int flag);
void      fft(double *data, int *nn, int ndim, int isign);
double  fmin_ell(double dl0s, double dos, double zs);
void      followi(struct point A, struct point B, struct point O, double dl0s, double dos, double zs);
void      follow(struct point A, struct point B, struct point O, double dl0s, double dos, double zs);
struct ellipse formeli(double a, double b, double c);
void      frprmn( double *p, int n, double ftol, int *iter, double *fret, double (*func)(double*), void (*dfunc)(double*, double*) );
void      fr_sq_point(struct point **square, int nbr_lin);
void      f_shape2(long int *istart, struct galaxie *liste, char *name);
void      f_shape3(long int *istart, struct galaxie *liste, char *name);
void      f_shape4(long int *istart, struct galaxie *liste, char *name);
void      f_shape_abs(long int *istart, struct galaxie *liste, char *name);
void      f_shape(long int *istart, struct galaxie *liste, char *name, int flag);
void      f_shmap(int *istart, struct shear *liste, char *name);
void      f_source(char *name, struct galaxie *galaxy, int *n);
double  fz_dlsds(double z);
struct galaxie unlens1(struct galaxie arclet, double dlsds);
void    g_amplif(int iampf, double z, char *file);
void    g_ampli(int iamp, int np, double z, char *file);
void    g_curv(int icurv, int np, double z, char *file1, char *file2, char *file3);
void    g_dpl(int idpl, int np, double z, char *filex, char *filey);
void    g_grid(int igrid, int ngrid, double zgrid);
void    g_mass(int imass, int np, double zl, double zs, char *file);
void    g_poten(int ipoten, int np, double z, char *file);
double  g_profil(double x, double y, struct galaxie gal);
void    g_prop(int nprop, double z);
void    g_radial(int iradial, double zradial, double theta);
void    grid();
void    gridp();
void    g_shearf(int ishear, double z, char *file, int nshearf);
void    g_shear(int ishear, int np, double z, char *file);
void    g_time(int flag, int np, double z, char *file);
int     getNzmlimit();
char*   getParName(int ipx, char* name, int type);
int     getNConstraints();
int     getNParameters();
void    getRADEC(char *line, int *iref, double *ra, double *dec );
double  hern_dpl(double r, double rs, double kappas);
double  hern_kappa_eps(double r, double rs, double theta, double kappas, double eps);
struct point hern_gamma_eps(double r, double rs, double theta, double kappas, double eps);
void    ic_product(double **r_im, double **i_im, double **r_im2, double **i_im2, double **r_prod, double **i_prod, int nbr_lin, int nbr_col);
complex icpx(complex c);
void      i_marker(char markfile[], double z);
void      imtosou(double zimage, char *sname);
int     inconvexe(struct point P, int np, struct point I[NPOINT]);
int     init_grille(char *infile, int noedit);
int     insidebord(struct point *P, struct triplet *T);
int     inside(struct point *P, struct triplet *T);
void    d_binning(double **im, int *nx, int *ny, int bin);
double  interpol(double xx, const double *fxx, const double *fyy, const double *fy2, int imax);
int     inverse(const struct point gsource[][NGGMAX], struct point *P, struct bitriplet *Tsol);
void    isoima(struct ellipse *es, struct ellipse *ampli, struct ellipse *ei);
double  iter(double phi0, double coeur, double ct);
void    keep_cl(double **image, struct pixlist *pl, int *npl);
int     lire(struct galaxie *liste, char *name);
complex lncpx(complex c);
void    local(int type, char localfile[]);
void    mdci05(double x, double y, double eps, double rc, double b0, struct matrix *res);
struct matrix mdci10(double x, double y, double eps, double rc, double b0);
struct matrix mdci15(double x, double y, double eps, double rc, double b0);
void    mdcsiemd(double x, double y, double eps, double b0, struct matrix *res);
struct matrix rotmatrix(struct matrix *P, double theta);
struct matrix sp_grad2(struct point P);
double  mass2d_NFW(double velocity_disp,double reference_rad, double scale_rad);      //  TV
double  mass3d_NFW(double velocity_disp,double reference_rad, double scale_rad);      //  TV
double  vel_Mampost(double velocity_disp,double reference_rad, double scale_rad, double redshift_0, double betta_index); //  TV     
double  min_parabol(double x0, double y0, double x1, double y1, double x2, double y2);
double  min_slope_par(double x0, double y0, double x1, double y1, double x2, double y2);
double  mean(int n, double *array);
double  median(int n, double *array);
double  mode(int n, double *array);
void    multiscale_grid(long int *pilens);
double  ncpx(complex c);
double  nfw_dpl(double r, double rs, double kappas);
double  nfw_gamma(double r, double rs, double kappas);
struct point nfw_gamma_eps(double r, double rs, double theta, double kappas, double eps);
double  nfwg_dpl(double r, double rs, double kappas, double alpha);
double  nfwg_gamma(double r, double rs, double kappas, double alpha);
struct point  nfwg_gamma_eps(double r, double rs, double theta, double kappas, double eps, double alpha);
double  nfwg_kappa(double r, double rs, double kappas, double alpha);
double  nfwg_kappa_eps(double r, double rs, double theta, double kappas, double eps, double alpha);
double  nfw_kappa(double r, double rs, double kappas);
double  nfw_kappa_eps(double r, double rs, double theta, double kappas, double eps);
void    e_nfw_rs2c(double sigma_s, double r_s, double *rhos, double *c, double *M_vir, double z);
void    e_nfw_cm200_sigrs( double c, double mvir, double *sigma_s, double *r_s, double z );
void    e_nfw_cr200_sigrs( double c, double rvir, double *sigma_s, double *r_s, double z );
void    e_nfw_crs2sig(double c, double rs, double *sigma_s, double z);
void    e_nfw_rsm200_sigrs( double rs, double mvir, double *sigma_s, double z );
void    e_nfw_rsr200_sigrs( double rs, double rvir, double *sigma_s, double z );
double  elli_tri(const struct pot *ilens);
void    e_nfw_c3D2c2D(double c3D,double *c2D);
int     o_big_slope();
double  o_chi();
void    o_chi_flux(struct galaxie *arclet, double fluxS, double *da, double *sig2flux);
int     o_chi_lhood0(double *chi2, double *lhood0, double *np_b0);
double  o_chi_pos();
void    o_chires(char *filename);
void    o_dmag(int i, struct galaxie *gali, double *da);
void    o_dpl(int n_im, struct galaxie *gali, struct point *ps, double *np_b0);
void    o_flux(int n, double *fluxS, int n_famille, double *np_b0);
double  o_get_err(int i, int ipx);
double  o_get_lens(int i, int ipx);
double  o_get_lmax(int i, int ipx);
double  o_get_lmin(int i, int ipx);
double  o_get_source(struct galaxie *source, int ipx);
void    o_global();
void    o_global_free();
double  **o_invim(double zim, double *drim, struct pixlist *plo, int *nplo);
void    o_keep_min(double x1, double x2, double y1, double y2, int ils, int ipx);
void    o_keepz_min(double x1, double x2, double y1, double y2, int iz);
double  o_lhood(int *error);
void    o_mag(int n, struct galaxie *arclet);
void    o_mag_m(int n, struct galaxie *arclet);
double  o_min_loc(double y0);
double  o_min_slope(double y0);
double  o_prep();
void    o_prep_mult(int ntmult, struct galaxie *mult);
void    o_print(FILE *OUT, double chi0);
void    o_print_res(double chi0, double evidence);
void    opt_source();
void    o_run1();
void    o_run2();
double  o_run_bayes();
void    o_run();
void    o_run_mc();
void    o_run_mc0();
//double  o_run_nest();
void    o_runpot1(int flag);
void    o_runpot2();
void    o_runz_mc();
void    o_runz_mc0();
void    o_scale_pot();
void    o_set_err(int i, int ipx, double x);
void    o_set_exc(int i, double excu[NLMAX][NPAMAX], double excd[NLMAX][NPAMAX], int block[NLMAX][NPAMAX]);
void    o_set_ip();
void    o_set_lens_bayes(int method, int prior, double limit);
void    o_set_lens(int i, int ipx, double x);
void    o_set_lmax(int i, int ipx, double x);
void    o_set_lmin(int i, int ipx, double x);
void    o_set_limit_bayes(int nParam, long int nVal, double **array, int gauss, double limit);
void    o_set_map();
void    o_set_map_mc();
void    o_set_map_z();
void    o_set_start(int i, int block[NLMAX][NPAMAX]);
void    o_set_source(struct galaxie *source, int ipx, double val);
void    o_shape(int i, struct galaxie *gali, double *dx, double *sigx2, double *dy, double *sigy2, double *da);
int     o_slope_sp(double *y0);
void    o_stat(int na, struct galaxie arclet[NAMAX]);
double  o_step(double chi0);
int     o_swi_big(int i, int ipx);
void    packzmlimit(char limages[ZMBOUND][IDSIZE], int nimages, char name[IDSIZE]);
complex pcpx(complex c1, complex c2);
complex pcpxflt(complex c, double f);
double  pi05(double x, double y, double eps, double rc, double b0);
double  pi10(double x, double y, double eps, double rc, double b0);
double  pi15(double x, double y, double eps, double rc, double b0);
struct point **al_sq_point(int nbr_lin, int nbr_col);
struct point barycentre(struct triplet *A);
struct point bcentlist(struct point *P, int n);
double  max(int n, double *array);
double  mean(int n, double *array);
double  median(int n, double *array);
struct point milieu(struct point *A, struct point *B);
double  min(int n, double *array);
struct point next(struct point A, struct point B, double dpl, double dl0s, double dos, double zs);
struct point rotation(struct point P, double theta);
struct point sp_grad(struct point P);
struct point wcenter(struct point A, double wa, struct point B, double wb);
struct point weight_baryc(struct point *P, struct galaxie *multi, int n, int n_famille);
struct polar polxy(struct point xy);
void    prep_non_param();
void    pro_arclet(long int n, struct galaxie *gal);
void    pro_arclet_m(int n, struct galaxie gal[NAMAX]);
double  psiemd(double x, double y, double eps, double b0);
double  ptau(double x, double y);
void    r_cleanlens(FILE *IN, FILE *OUT);
void    r_cline(FILE *IN, FILE *OUT);
void    r_cosmolimit(FILE *IN, FILE *OUT);
void    r_cosmologie(FILE *IN, FILE *OUT);
void    r_frame(FILE *IN, FILE *OUT);
void    r_vfield(FILE *IN, FILE *OUT);
void    r_vfieldlimit(FILE *IN, FILE *OUT);
void    r_grille(FILE *IN, FILE *OUT);
void    r_image(FILE *IN, FILE *OUT);
void    r_large(FILE *IN, FILE *OUT);
void    r_limit(FILE *IN, FILE *OUT, int i);
void    r_msgrid(FILE *IN, FILE *OUT);
void    r_observ(FILE *IN, FILE *OUT);
void    r_potentiel(FILE *IN, FILE *OUT, int ilens);
void    r_potfile(FILE *IN, FILE *OUT, struct g_pot *pot);
void    r_dynfile(FILE *IN,FILE *OUT);   //TV Oct2011
void    r_runmode(FILE *IN, FILE *OUT);
void    r_shapelimit(FILE *IN, FILE *OUT, long int isrc);
void    r_shapemodel(FILE *IN, FILE *OUT, long int isrc);
void    r_source(FILE *IN, FILE *OUT);
double  **readBayesModels(int *nParam, long int *nVal);
void    readConstraints();
double  **readimage(struct g_pixel *P);
void    read_lenstable();
int     rescaleCube_1Atom(double *cube, int npar);
double  rho_cri(double z);
void    s_compmag(struct galaxie *gal, int *idum);
void    scale_pot(struct g_pot *pot);
complex scpx(complex c1, complex c2);
complex scpxflt(complex c1, double f2);
void    rhos2b0();
void      set_default();
void      set_dynamics(long int i);  // in set_lens.c
void      set_lens();
void      set_lens_par(FILE *OUT);
int     set_potfile(int nplens, struct g_pot *pot);
void      set_res_par();
void      setBayesModel( long int iVal, long int nVal, double **array);// in readBayesModels.c
void    setScalingRelations(struct g_pot *pot);
double  sgn_darg(complex z1, complex z2);
double  sgn(double x);
double  sig2posS4j(double sigmaArcsec, struct galaxie *multi, struct point ps, int j);
double  sig2posS(double sigmaArcsec, int n, int n_famille, double *np_b0);
double  sig2posSj(double sigmaArcsec, struct galaxie *multi, int n, int j, int n_famille);
void    sig2posSe(struct galaxie *multi, double *sigx2, double *sigy2);
double  signe(const struct triplet *T, const struct point *P);
complex sincpx(complex c);
complex sinhcpx(complex c);
double  slope(double x0, double y0, double x1, double y1);
void    sortf(int n, double *A, int (*comp)(double, double));
void    sort(long int n, struct galaxie *A, int (*comp)(struct galaxie *, struct galaxie *));
void    s_pixlist(double **image, struct pixlist *pl, int *npl);
void    spline(const double *x, const double *y, int n, double yp1, double ypn, double *y2);
int     split_image(double **image, int nbr_lin, int nbr_col);
double  sp_pot(struct point P);
void    sp_set(double **map, int nx, int ny, double **map_xx, double **map_yy);
complex sqcpx(complex c);
complex sqrtcpx(complex c);
void    s_sof();
void    s_source();
void    s_sourcebox(struct g_pixel *ps, const char *centerfile, double dlsds);
int     splitzmlimit(char name[IDSIZE], char limage[ZMBOUND][IDSIZE]);
void      st_opt(int mode, struct galaxie *imas, int nima, double *z, double ts[NASMAX][200], double pz[NASMAX][200], double ipz[NASMAX][200], double amp[NASMAX][200], double *kappa, double *gam, double *thetap, double *sumpz, double *tauix, double *tauiy, int jmax, double *zpm, double *nmag);
double  stddev(int n, double *array);
void    study_pg(int type, double seeing, char studyfile[], int fake);
complex tancpx(complex c);
void    tirer(struct galaxie *source);
double  tnfwg_dpl(double xwant, double alphawant);
double  tnfwg_kappa(double xwant, double alphawant);
void    tracepot();
void    e_transform(struct triplet *I, double dlsds, struct triplet *S);
void    interieur(const struct triplet *E, struct triplet *I);
int     unlens_bc(const struct point *Ps, struct point Bs, struct galaxie *multi, struct point *multib, int n, int n_famille, int *det_stop);
void    updatecut(int i);
void    update_emass(int i);
void    update_epot(int i, double *epot);
void    w_critic();
void    w_gianti(int ngs, char name[20]);
void    w_prop(int nprop, char propfile[]);
void    wr_mass(char *name, double **map_xx, double **map_yy);
void    wr_pot(char *name, double **map);
void    w_sicat(struct galaxie ima[NFMAX][NIMAX], struct galaxie ss[NFMAX]);
double  zero(double c1, double c2, double (*f)(double));
double  zero_t(double c1, double c2, double (*f)(double));
void    zonemult();

//partie einasto

float einasto_masse( double r, double rs, int n, double rhos);
float einasto_sigma( double r, double rs, int n, double rhos);
float einasto_kappa( double r, double rs, int n, double rhos, double kappa_crit);
float einasto_kappa_av( double r, double rs, int n, double rhos, double kappa_crit);
double einasto_gamma( double r, double rs, int n, double rhos, double kappa_crit);
struct point   einasto_gamma_eps(double r, double rs, double n, double theta, double kappa_crit, double rhos,double eps);
float einasto_phi(double r, double rs, int n, double rhos, double kappa_crit);
float einasto_alpha( double r, double rs, int n, double rhos, double kappa_crit);
void read_table_einasto();

void precompsource();
void testpixel();

//#ifdef __cplusplus
//};
//#endif

#endif // if FONCTION_H