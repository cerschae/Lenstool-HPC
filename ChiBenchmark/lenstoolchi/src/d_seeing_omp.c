//d_seeing_omp.c
//
//  name of main function: d_seeing_omp    
//  author:     Sergey Rodionov  (astroseger@gmail.com)          
//  date:       01/2012                    
//  place:      Marseille                  
//                                             

//function for fast add seeing to images by means of fft 
//principal function is the last function
//ATTENTION!!! WE USE GLOBAL VARIABLES FOR OPTIMISATION REASONS
//so don't remove check_not_in_parallel();

//because historical reason we have nx, ny swapped in respect to lenstool
//here -- im[0..nx][0..ny] ... in lenstool im[0..ny][0..nx]


#include <gsl/gsl_fft_complex.h>
#include <stdio.h>
#include <stdlib.h>
#include "structure.h"
#include "lt.h"

struct cpa2d;  //predifiniton of cpa2d (see below)

//ATTENTION!!!! GLOBAL VARIABLES

//fourier transformation of psf
//psf_f->nx, psf_f->ny --- sizes of our image!
static struct cpa2d* psf_f = NULL;  

//preallocated "padded" image
//(we only us it directly from d_seeing_omp)
static struct cpa2d* pim_global = NULL; 

//structure which keep parameters for which  psf_f as constructed 
//(just to check that it haven't been changed)
static struct g_observ O_psf_f; 

//nx, ny and scale for first run (in addition to O_psf_f)

static int nx_init, ny_init;
static double scale_init;

//wavetables and workspaces for fourier trasform 
// *_x - for psf->nx
// *_y - for psf->ny
//for each row or column --- separeated workspace
static gsl_fft_complex_wavetable * wt_x;
static gsl_fft_complex_wavetable * wt_y;
static gsl_fft_complex_workspace **ws_x; //size psf_f->ny !!! (reverse!)
static gsl_fft_complex_workspace **ws_y; //size psf_f->nx !!! (reserse!)

//gsl "complex packed array" in 2D. For 2D fourier trasform
struct cpa2d
{
   double* v;
   int nx;
   int ny;
};
//                                                                            
//several basic function for the access to cpa2d, and allocate/free of it 
static struct cpa2d* cpa2d_alloc(int nx, int ny)
{
   struct cpa2d* c = (struct cpa2d*) malloc(sizeof(struct cpa2d*)); 
   c->v = (double*) malloc(sizeof(double) * nx * ny * 2);
   c->nx = nx;
   c->ny = ny;
   return c;
}
//                                                                            
static void cpa2d_free(struct cpa2d* c) 
{
   free(c->v);
   free(c);
}
//                                                                           
inline void cpa2d_set(struct cpa2d* c, int i, int j, double re, double im)
{
   c->v[(i * c->ny + j) * 2]         = re;
   c->v[(i * c->ny + j) * 2 + 1]     = im;
}
//                                                                           
inline void cpa2d_set_re(struct cpa2d* c, int i, int j, double re)
{
   c->v[(i * c->ny + j) * 2]         = re;
}
//                                                                           
inline double cpa2d_get_re(struct cpa2d* c, int i, int j)
{
   return  c->v[(i * c->ny + j) * 2];
}
//                                                                           
inline double cpa2d_get_im(struct cpa2d* c, int i, int j)
{
   return c->v[(i * c->ny + j) * 2 + 1];
}
//                                                                            

//from im (nx,ny) to pim (pim->nx,pim->ny) with zero padding
static void zero_padding(double **im,    int nx,     int ny, 
			 struct cpa2d* pim)
{
   if (pim->nx <= nx || pim->ny <= ny)
     {
	fprintf(stderr, "error | zero_paddinig | pim->nx <= nx || pim->ny <= ny\n");
	exit(EXIT_FAILURE);
     }
   int px = (pim->nx - nx) / 2;
   int py = (pim->ny - ny) / 2;
   
   int i;
   
#pragma omp parallel for schedule(static)  
   for (i = 0 ; i < pim->nx ; i++)
     {
	int j;
	for (j = 0 ; j < pim->ny ; j++)
	  {
	     if ((i >= px)  && (i < px + nx) && 
		 (j >= py)  && (j < py + ny))
	       {
		  cpa2d_set(pim, i, j, im[i - px][j - py], 0);
	       }
	     else
	       {
		  cpa2d_set(pim, i, j, 0, 0);
	       }
	  }   
     }
}


//from pim to im (unpadding and "unsplitting")
//We take into account "splitting" in pim
//if cntrx, cntry - is position of maximum in psf
//and it is center of "splitting"
static void split_unpadding(double **im,    int nx,     int ny, 
			    struct cpa2d* pim)
{
   int px = (pim->nx - nx) / 2;
   int py = (pim->ny - ny) / 2;   
   

   //calculate position of center in psf
   //for odd center in n/2
   //for even center in (n - 0.5)/2
   int cntrx = (int)floor((pim->nx - 0.5) / 2.0);
   int cntry = (int) floor((pim->ny - 0.5) / 2.0);
   
   int i;
   
#pragma omp parallel for schedule(static)  
   for (i = px ; i < px + nx ; i++)
     {
	int j;
	for (j = py ; j < py + ny ; j++)
	  {
	     int ip, jp;
	     
	     ip = i;
	     jp = j;
	     
	     //splitting
	     if (i < pim->nx - cntrx)
	       ip = i + cntrx;
	     else
	       ip = i - (pim->nx - cntrx);
	     
	     if (j < pim->ny - cntry)
	       jp = j + cntry;
	     else
	       jp = j - (pim->ny - cntry);	  
	     
	     //unpadding
	     im[i - px][j - py] = cpa2d_get_re(pim, ip, jp);
	  }
     }
}
//two simple function for 2D fast fourier trasform (direct and invert)
//ATENTIANON this two function use global variables 
static void fft_2D_omp(struct cpa2d* data)
{   
   int i;   

#pragma omp parallel for schedule(static)
   for (i = 0 ; i < data->nx ; i++)
     {		
	gsl_fft_complex_forward(data->v + data->ny * i * 2, 1, data->ny, 
				wt_y, ws_y[i]);
     }

#pragma omp parallel for schedule(static)
   for (i = 0 ; i < data->ny ; i++)
     {
	gsl_fft_complex_forward(data->v + i * 2, data->ny, 
				data->nx, wt_x, ws_x[i]);
     }
}
//                                                                            
static void fft_2D_inv_omp(struct cpa2d* data)
{
   int i;

#pragma omp parallel for schedule(static)
   for (i = 0 ; i < data->nx ; i++)
     {
	gsl_fft_complex_inverse(data->v + data->ny * i * 2, 1, 
				data->ny, wt_y, ws_y[i]);
     }
#pragma omp parallel for schedule(static)   
   for (i = 0 ; i < data->ny ; i++)
     {
	gsl_fft_complex_inverse(data->v + i * 2, data->ny, 
				data->nx, wt_x, ws_x[i]);
     }
}
//                                                                           Ñ

//we choose smallest number greater than n for which gsl fft is fast 
//(biggest factor should be equal to 7. see gsl documentation)
static int choose_fast_fft_n(int n)
{
   const int primes[4] = {2, 3, 5, 7};
   const int primes_n = 4;   
   
   int det = 0;
   do
     {
	int tmp = n;
	int i;
	for (i = 0 ; i < primes_n ; i++)
	  while (tmp % primes[i] == 0)
	    tmp /= primes[i];
	if (tmp == 1)
	  det = 1;
	else
	  n++;
     }
   while (!det);

   return n;   
}
//                                                                            
static void create_psf_psffile(double***spsf, int* nx, int* ny)
{
   char* header;
   const extern struct g_observ  O;
   
   //again nx and ny should be swapped here
   *spsf = rdf_fits((char*)O.psffile, ny, nx, &header);
   free(header);
   //check maximum (it should be at nx/2 ny/2)
   int max_i = *nx / 2;
   int max_j = *ny / 2;
   int i,j;
   for (i = 0 ; i < *nx ;  i++)
     for (j = 0 ; j < *ny; j++)
       {
	  if ((*spsf)[i][j] > (*spsf)[max_i][max_j])
	    {
	       fprintf(stderr, "Error | d_seeing_omp/create_psf_psffile | maximum is not located in center of psf\n");
	       exit(EXIT_FAILURE);
	    }   
       }
   
   //normalization
   long double sum = 0;
   for (i = 0 ; i < *nx ; i++)
     for (j = 0 ; j < *ny ; j++)
	  sum += (*spsf)[i][j];
   
   for (i = 0 ; i < *nx ; i++)
     for (j = 0 ; j < *ny ; j++)
       {
	  (*spsf)[i][j] /= sum;
       }
}
//                                                                            
static void create_psf_create_filtre(double scale, int big_n, double***spsf, int* nx, int*ny)
{
   int i,j;
   const double PSF_ZERO = 1e-20;
   const extern struct g_observ  O;
   
   //create_filtre work correctly only for even n
   if (big_n % 2 != 0)
     big_n++;      
   
   double** psf = (double **) alloc_square_double(big_n, big_n);
      
   if (O.filtre)
     {
	fprintf(stderr, "Error | d_seeing_omp.c/create_psf | O.filtre != 0 we don't know what it means");
	exit(EXIT_FAILURE);
     }
   if (O.setseeing == 1)
     {
	crea_filtre(O.seeing, scale, psf, big_n);  
     }
   else if (O.setseeing == 2) 
     {
	crea_filtre_e(O.seeing_a, O.seeing_b, O.seeing_angle, scale, psf, big_n); 
     }
   else
     {
	fprintf(stderr, "Error | d_seeing_omp.c/create_psf | unknown O.setseeing = %d", O.setseeing);
	exit(EXIT_FAILURE);
     }
      
   //just in case we check maximum value (it should be in big_n/2)
   int max_i = big_n / 2;
   for (i = 0 ; i < big_n ;  i++)
     for (j = 0 ; j < big_n; j++)
       if (psf[i][j] > psf[max_i][max_i])
	 {
	    fprintf(stderr, "Unexpected behavior of crea_filtre\n");
	    exit(EXIT_FAILURE);
	 }

   //now we can try to shrink psf
   int si = 1;
   int sj = 1;
   
   int det = 1; 

   do
     {
	for (j = 0 ; j < big_n ; j++)
	  {
	     if (psf[si][j] > PSF_ZERO || psf[2 * max_i - si][j] > PSF_ZERO)
	       det = 0;
	  }
	if (det)
	  si++;
	if (si == max_i)
	  {
	     fprintf(stderr, "Unexpected behavior 3 of crea_filtre\n");
	     exit(EXIT_FAILURE);
	  }
     }
   while (det);

   det = 1;
   do
     {
	for (i = 0 ; i < big_n ; i++)
	  {
	     if (psf[i][sj] > PSF_ZERO || psf[i][2 * max_i - sj] > PSF_ZERO)
	       det = 0;	
	  }
	if (det)
	  sj++;
	if (sj == max_i)
	  {
	     fprintf(stderr, "Unexpected behavior 3 of crea_filtre\n");
	     exit(EXIT_FAILURE);
	  }
     }
   while (det);


   *nx = (max_i - si) * 2 + 1;
   *ny = (max_i - sj) * 2 + 1;
   
   
   (*spsf) = (double **) alloc_square_double(*nx, *ny);
   
   long double sum = 0;
   for (i = 0 ; i < *nx ; i++)
     for (j = 0 ; j < *ny ; j++)
       {
	  (*spsf)[i][j] = psf[si + i][sj + j];
	  sum += (*spsf)[i][j];
       }
   for (i = 0 ; i < *nx ; i++)
     for (j = 0 ; j < *ny ; j++)
       {
	  (*spsf)[i][j] /= sum;
       }
   free_square_double(psf, big_n);   
}
//                                                                            
static void create_psf(double scale, int big_n, double***spsf, int* nx, int*ny)
{
   const extern struct g_observ  O;
      
   //we also save parameters for which psf was constructed
   O_psf_f = O; //we cannot use O_psf_f.psfile anymore!x
   
   if (O.setseeing == 1 || O.setseeing == 2)
     {
	create_psf_create_filtre(scale, big_n, spsf, nx, ny);
     }
   else if (O.setseeing == 3)
     {
	create_psf_psffile(spsf, nx, ny);   
     }
   else
     {
	fprintf(stderr, "Error | d_seeing_omp.c/create_psf | unknown O.setseeing = %d", O.setseeing);
	exit(EXIT_FAILURE);
     }
}

//                                                                            
//because historical reason we have nx, ny swapped in respect to lenstool
//here -- im[0..nx][0..ny] ... in lenstool im[0..ny][0..nx]
//ATTENTIAON!!! because of these nx ny was swapped here !!!
void d_seeing_omp(double **im, int ny, int nx, double scale)
{
   //we must not be in parallel now   
   check_not_in_parallel("d_seeing_omp");
   
   
   int i;
     
   if (psf_f == NULL)  //first run
     {	
	double** psf;
	int nx_psf,ny_psf;
	
	//create psf and save parameters for which it was constructed
	create_psf(scale, nx + 1, &psf, &nx_psf, &ny_psf);
	
	//save value of nx and ny of first run
	nx_init = nx;
	ny_init = ny;
	
	//save scale parameter of first run
	scale_init = scale;
	
	//we should add "zero" padding of size at least nx_psf/2 (ny_psf/2)
	//this value is exactly equal to size of wings of psf...
	//and after we search appropriative size of array for which fft 
	//could be calculated fast
	int nx1 = choose_fast_fft_n(nx + nx_psf / 2);
	int ny1 = choose_fast_fft_n(ny + ny_psf / 2);
//	printf("nx=%d ny=%d nx_psf=%d ny_psf=%d nx1=%d ny1=%d\n", 
//	       nx, ny, nx_psf, ny_psf, nx1, ny1);
	
	//now we can allocated wavetables and workspaces for fft
	wt_x = gsl_fft_complex_wavetable_alloc (nx1);
	wt_y = gsl_fft_complex_wavetable_alloc (ny1);
	
	ws_x = (gsl_fft_complex_workspace**) 
	  malloc(sizeof(gsl_fft_complex_workspace*) * ny1);
	ws_y = (gsl_fft_complex_workspace**) 
	  malloc(sizeof(gsl_fft_complex_workspace*) * nx1);
	
	for (i = 0 ; i < ny1 ; i++)
	  ws_x[i] = gsl_fft_complex_workspace_alloc (nx1);
	for (i = 0 ; i < nx1 ; i++)
	  ws_y[i] = gsl_fft_complex_workspace_alloc (ny1);		
	
	//calculate fft of psf
	psf_f = cpa2d_alloc(nx1, ny1);
	zero_padding(psf, nx_psf, ny_psf, psf_f);	
	fft_2D_omp(psf_f);
	
	//and finally allocated padded image with which we will work
	pim_global = cpa2d_alloc(psf_f->nx, psf_f->ny);
	
	//real psf we will not need anymore
	free_square_double(psf, nx_psf);	
     }
   else
     {
	//we just check that noting have changed after initialisation
	const extern struct g_observ  O;
	if (nx != nx_init || ny != ny_init || scale != scale_init)
	  {
	     fprintf(stderr, "Error | d_seeing_omp.c/d_seeing_omp | nx, ny or scale have been changed from first run");
	     exit(EXIT_FAILURE);
	  }
	if (O.setseeing   != O_psf_f.setseeing ||
	    O.seeing      != O_psf_f.seeing ||
	    O.seeing_a     != O_psf_f.seeing_a ||
	    O.seeing_b     != O_psf_f.seeing_b || 
	    O.seeing_angle != O_psf_f.seeing_angle )
	  {
	     fprintf(stderr, "Error | d_seeing_omp.c/d_seeing_omp | O structure  have been changed from first run");
	     exit(EXIT_FAILURE);
	  }
     }
   
   //Now we can start main part of algorim
   //psf_f->nx and psf_f->ny keep dimension of padded image   
      
   struct cpa2d* pim = pim_global;  //just for simplicity

   zero_padding(im, nx, ny, pim);
   
   fft_2D_omp(pim);
   
   //simple multiplication of pim and psf_f
   //(see basic theory of convolution)
#pragma omp parallel for schedule(static)
   for (i = 0 ; i < pim->nx * pim->ny * 2 ; i += 2)
     {
	double rr = pim->v[i] * psf_f->v[i]    - pim->v[i + 1] * psf_f->v[i + 1];
	double ii = pim->v[i] * psf_f->v[i + 1] + pim->v[i + 1] * psf_f->v[i];	
	pim->v[i] = rr;
	pim->v[i + 1] = ii;
     }
   
   fft_2D_inv_omp(pim);
   
   split_unpadding(im, nx, ny, pim);       
}
