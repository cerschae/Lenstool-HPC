#ifndef LT_H_
#define LT_H_
#include <stdio.h>

/* Functions declaration */

double 	**rdf_fits(char *filename,int *nx,int *ny,char **pheader);
double	***alloc_cubic_double(int nbr_lin,int nbr_col,int nbr_slice);
int	***alloc_cubic_int(int nbr_lin,int nbr_col,int nbr_slice);
double 	**alloc_square_double(int nbr_lin,int nbr_col);
int 	**alloc_square_int(int nbr_lin,int nbr_col);
double	*alloc_vector_double(int nbr_lin);
void	arraySort(double *a, int len);
void    cholesky(double **a, int n, double *p);
void    chop(char *str);
double 	**dmatrix(int nrl,int nrh,int ncl,int nch);
double 	*dvector(int nl,int nh);
void    free_square_double(double **square,int nbr_lin);
void    free_square_int(int **square,int nbr_lin);
void 	nrerror(char *error_text);
double 	*vector(int nl,int nh);
int 	*ivector(int nl,int nh);
int     indexCmp(const char *str1, const char *str2);
double 	**matrix(int nrl,int nrh,int ncl,int nch);
int 	**imatrix(int nrl,int nrh,int ncl,int nch);
double 	**submatrix(double **a,int oldrl,int oldrh,int oldcl,int newrl,int newcl);
void    fblanc(FILE *fichier, int n);
void    flire(FILE *fichier, char phrase[]);
void    fmot(FILE *fichier, char mot[]);
void    ftab(FILE *fichier, int n);
void 	free_vector(double *v,int nl);
void 	free_ivector(int *v,int nl);
void 	free_dvector(double *v,int nl);
void 	free_matrix(double **m,int nrl,int nrh,int ncl);
void 	free_dmatrix(double **m,int nrl,int nrh,int ncl);
void 	free_imatrix(int **m,int nrl,int nrh,int ncl);
void 	free_submatrix(double **b,int nrl);
double 	**convert_matrix(double *a,int nrl,int nrh,int ncl,int nch);
void 	free_convert_matrix(double *b,int nrl);
double	kth_smallest(double a[], int n, int k);
int     getWords(char *str);
void 	polint(double *xa,double *ya,int n,double x,double *y,double *dy);
double three_eighths(double a, double b, int n,double (*func)(double));       //THIS IS MINE
double three_eighths2(double a, double b,double c,double d, int n,double (*func)(double, double, double));       //THIS IS MINE
long int	rd_ipx(	char *file,int dim_r,int size[4],char *type_r,char mode[4],char *nature_r,
					char comments[1024],double *xmin,double *xmax,double *ymin,double *ymax );
long int	rd_ipxs(char *file,int dim_r,int size[4],char *type_r,char mode[4],char *nature_r,char comments[1024]);
double 	** rdf_fits(char *filename,int *nx,int *ny,char **pheader);
double 	**rdf_ipx(char *file,int *nx,int *ny,char *type,char *mode,char *nature,char comments[1024],
				double *xmin,double *xmax,double *ymin,double *ymax);
double 	**rdf_ipxs(char *file,int *nx,int *ny,char *type,char *mode,char *nature,char comments[1024]);
char*   upcase(char * str);
long int  wc(FILE * file);
int 	wrf_fits(char *filename,double **ima,int nx,int ny,double xmin,double xmax,double ymin,double ymax);
int 	wri_fits_abs( char *filename,int **ima,int nx,int ny,double xmin,double xmax,double ymin,double ymax, double ra,double dec);
int 	wrf_fits_abs(char *filename,double **ima,int nx,int ny,double xmin,double xmax,double ymin,double ymax, double ra,double dec);
int     wrf_tfits(char *filename, double **tbl, char* ttype[], int nx, long int ny);
int 	wri_fits(char *filename,int **ima,int nx,int ny,double xmin,double xmax,double ymin,double ymax);
int 	wrf_cube_fits(char *filename,double ***cube,int nx,int ny,int nz,double xmin,double xmax,double ymin,double ymax,double lmin, double lmax);
void	wr_ipx(	char *file,long int data,int dim,int *size,char *type,char *mode,char *nature,char *comments,
				double xmin,double xmax,double ymin,double ymax);
void	wr_ipxs(char *file,long int data,int dim,int *size,char *type,char *mode,char *nature,char *comments);
void	wrf_ipx(char *file,long int data,int dim,int *size,char *type,char *mode,char *nature,char *comments,
				double xmin,double xmax,double ymin,double ymax);
void 	printerror( int status);

#define median_WIRTH(n,a) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))

#endif /*LT_H_*/
