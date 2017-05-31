#include <iostream>
#include <string.h>
//#include <cuda_runtime.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <immintrin.h>

//#include "simd_math.h"
#include "chi_GPU.cuh"
#include <chi_CPU.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void chi_bruteforce_SOA_GPU_grid_gradient(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images){

double dx,dy,x_pos,y_pos;        //pixelsize
dx = (frame->xmax - frame->xmin)/(runmode->nbgridcells-1);
dy = (frame->ymax - frame->ymin)/(runmode->nbgridcells-1);
int index=0;       // index of the image within the total image array
int grid_dim;
double im_dist[MAXIMPERSOURCE]; // distance to a real image for an "theoretical" image found from a theoretical source
int im_index;       // index of the closest real image for an image found from a theoretical source
int second_closest_id;   // index of the second closest real image for an image found from a theoretical source
int thread_found_image = 0;   // at each iteration in the lens plane, turn to 1 whether the thread find an image
struct point im_position, temp; // position of the image found from a theoretical source + temp variable for comparison
struct triplet Tsup, Tinf, Tsupsource, Tinfsource;// triangles for the computation of the images created by the theoretical sources
struct galaxy sources[runmode->nsets]; // theoretical sources (common for a set)
int nimagesfound[runmode->nimagestot][MAXIMPERSOURCE]; // number of images found from the theoretical sources
struct point tim[runmode->nimagestot][MAXIMPERSOURCE]; // theoretical images (computed from sources)


double *grid_gradient_x, *grid_gradient_y;

grid_gradient_x = (double *)malloc((int) (runmode->nbgridcells) * (runmode->nbgridcells) * sizeof(double));
grid_gradient_y = (double *)malloc((int) (runmode->nbgridcells) * (runmode->nbgridcells) * sizeof(double));
grid_dim = runmode->nbgridcells;
//Packaging the image to sourceplane conversion
gradient_grid_GPU_sorted(grid_gradient_x,grid_gradient_y,frame,lens,runmode->nhalos,grid_dim);

//printf ("%f %f \n", grid_srcplane_x[0],grid_srcplane_y[0]);

index = 0;
*chi = 0;

for( int  source_id = 0; source_id < runmode->nsets; source_id ++){

//////////////////////////////////////Initialisation//////////////////////////////////////
for (int i=0; i < nimages_strongLensing[source_id]; ++i){
  for (int j=0; j < nimages_strongLensing[source_id]; ++j){
	nimagesfound[i][j] = 0;
  }
}

for( int image_id = 0; image_id < nimages_strongLensing[source_id]; image_id++){

  //////////////////////////////////////computation of theoretical sources//////////////////////////////////////
  chi_transformImageToSourcePlane_SOA(runmode->nhalos, &images[index+image_id].center,images[index+image_id].dr,lens,&sources[source_id].center);

  //if (DEBUG ==1 )
	//printf("index %d image_id %d source_id %d %f \n",index, image_id, source_id,images[index+image_id].redshift);
  sources[source_id].redshift = images[index+image_id].redshift;
  sources[source_id].dr = images[index+image_id].dr;
  sources[source_id].dls = images[index+image_id].dls;
  sources[source_id].dos = images[index+image_id].dos;

  for (int x_id = 0; x_id < runmode->nbgridcells-1; ++x_id ){
	for (int y_id = 0; y_id < runmode->nbgridcells-1; ++y_id ){

	  x_pos = frame->xmin + x_id * dx;
	  y_pos = frame->ymin + y_id * dy;

	  // Define the upper + lower triangle, both together = square = pixel
	  Tsup.a.x=x_pos;
	  Tsup.b.x=x_pos;
	  Tsup.c.x=x_pos+dx;
	  Tinf.a.x=x_pos+dx;
	  Tinf.b.x=x_pos+dx;
	  Tinf.c.x=x_pos;

	  Tsup.a.y=y_pos;
	  Tsup.b.y=y_pos+dy;
	  Tsup.c.y=y_pos;
	  Tinf.a.y=y_pos+dy;
	  Tinf.b.y=y_pos;
	  Tinf.c.y=y_pos+dy;

	  // Lens to Sourceplane conversion of triangles

	  chi_transformtriangleImageToSourcePlane_SOA_grid_gradient_upper(&Tsup, sources[source_id].dr,&Tsupsource,grid_gradient_x,grid_gradient_y,x_id*grid_dim+y_id, grid_dim);
	  chi_transformtriangleImageToSourcePlane_SOA_grid_gradient_lower(&Tinf, sources[source_id].dr, &Tinfsource,grid_gradient_x,grid_gradient_y,x_id*grid_dim+y_id, grid_dim);

	  thread_found_image=0;
	  if(chi_insideborder(&sources[source_id].center,&Tsupsource)==1){
		thread_found_image=1; // thread has just found an image
		im_index=0;
		im_position=chi_barycenter(&Tsup);
		im_dist[im_index]=chi_dist(im_position,images[index+im_index].center);  // get the distance to the real image
		for(int i=1; i<nimages_strongLensing[source_id]; i++){  // get the distance to each real image and keep the index of the closest real image
		  im_dist[i]=chi_dist(im_position,images[index+i].center);
		  if(im_dist[i]<im_dist[im_index]){
			im_index=i;
		  }
		  //printf(" im_index %d im_dist actual %f im_dist %f \n",im_index, im_dist[im_index], im_dist[i]);
		}
	  }

	  if(chi_inside(&sources[source_id].center,&Tinfsource)==1){
		thread_found_image=1; // thread has just found an image
		im_index=0;
		im_position=chi_barycenter(&Tinf);  // get the barycenter of the triangle
		im_dist[im_index]=chi_dist(im_position,images[index+im_index].center);  // get the distance to the real image
		for(int i=1; i<nimages_strongLensing[source_id]; i++){  // get the distance to each real image and keep the index of the closest real image

		  im_dist[i]=chi_dist(im_position,images[index+i].center);
		  //printf("im_dist[i] %f, im_position %f %f , images[index+im_index].center %f %f\n",im_dist[i], im_position.x,im_position.y, images[index+i].center.x,images[index+i].center.y);
		  if(im_dist[i]<im_dist[im_index]){
			im_index=i;
		  }
		  //printf(" im_index %d im_dist actual %f im_dist %f \n",im_index, im_dist[im_index], im_dist[i]);
		}
	  }

	  int skip_image = 0;

	  if(thread_found_image == 1 ){
		skip_image = 0;
		// Sometimes due to the numerical errors at the centerpoint, for SIE potentials an additional image will appear at the center of the Potential.
		// This is due to the fact that it is not possible to simulate an infinity value at the center correctly, Check that sis correspond to Nlens[0]
		/*
		for (int iterator=0; iterator < runmode->Nlens[0]; ++iterator){
		  //printf("lens[i].type %d %f %f %f \n",lens[i].type, fabs(im_position.x - lens[i].position.x) ,  fabs(im_position.y - lens[i].position.y), dx/2.);
		  if ( fabs(im_position.x - lens[0].position_x[iterator]) <= dx/2. and fabs(im_position.y  - lens[0].position_y[iterator]) <= dx/2.){
			skip_image = 1;
			printf("WARNING: You are using SIE potentials. An image to close to one of the potential centers has been classified as numerical error and removed \n");
		  }
		}*/
		if(skip_image==0){
			  //checking whether a closest image has already been found
		  //printf("imagenumber %d im_index %d , im_position.x %f , im_position.y %f \n", image_id, im_index  , im_position.x  , im_position.y);
		  if(nimagesfound[image_id][im_index]==0){ // if no image found up to now
			tim[image_id][im_index]=im_position;  //image position is allocated to theoretical image
			nimagesfound[image_id][im_index]++;
		  }
		  else if(nimagesfound[image_id][im_index]>0){ // if we have already found an image
			// If the new image we found is closer than the previous image
			if(im_dist[im_index]<chi_dist(images[index+im_index].center,tim[image_id][im_index]))
			{
			  temp=tim[image_id][im_index]; // we store the position of the old image in temp
			  tim[image_id][im_index]=im_position; // we link the observed image with the image we just found
			}
			else
			{
			  temp=im_position; // we store the position of the image we just found in temp
			}
			// initialising second_closest_id to the highest value
			// Loop over all images in the set except the closest one
			// and initialize to the furthest away image
			second_closest_id=0;
			for(int i=1; i<nimages_strongLensing[source_id] && i!=im_index; i++)
			{
			  if(im_dist[i]>im_dist[second_closest_id]) second_closest_id=i;
			}
			///////////////////////////////////////////////////////////////
			// Loop over all images in the set that are not yet allocated to a theoretical image
			// and allocate the closest one
			for(int i=0; i<nimages_strongLensing[source_id] && nimagesfound[image_id][i]==0; i++) // we search for an observed image not already linked (nimagesfound=0)
			{
			  if(im_dist[i]<im_dist[second_closest_id])
			  {
				second_closest_id=i;
				im_index=i; // im_index value changes only if we found a not linked yet image
				tim[image_id][im_index]=temp; // if we found an observed and not already linked image, we allocate the theoretical image temp
			  }
			}
			nimagesfound[image_id][im_index]++; // increasing the total number of images found (If we find more than 1 theoretical image linked to 1 real image, these theoretical
							// images are included in this number)
		  }


		}
		thread_found_image=0; // for next iteration

	  }
	}
  }

}


//////////////////////////////////////computing the local chi square//////////////////////////////////////
double chiimage;

for( int iter = 0; iter < nimages_strongLensing[source_id]*nimages_strongLensing[source_id]; iter++){
  int i=iter/nimages_strongLensing[source_id];
  int j=iter % nimages_strongLensing[source_id];

  if(i!=j){
	// In the current method, we get the source in the source plane by ray tracing image in nimagesfound[i][i]. If we ray trace back,
	// we arrive again at the same position and thus the chi2 from it is 0. Thus we do not calculate the chi2 (-> if i!=j)
	if(nimagesfound[i][j]>0){
	  chiimage=pow(images[index+j].center.x-tim[i][j].x,2)+pow(images[index+j].center.y-tim[i][j].y,2);  // compute the chi2
	  *chi += chiimage;
	}
	else if(nimagesfound[i][j]==0){
	  // If we do not find a correpsonding image, we add a big value to the chi2 to disfavor the model
	  *chi += 100.*nimages_strongLensing[source_id];
	}
  }
}
/*
for (int i=0; i < nimages_strongLensing[source_id]; ++i){
  for (int j=0; j < nimages_strongLensing[source_id]; ++j){
	printf(" %d",nimagesfound[i][j]);
	}
  printf("\n");
  }*/

//Incrementing Index: Images already treated by previous source_id
index+=nimages_strongLensing[source_id];
}

free(grid_gradient_x);
free(grid_gradient_y);
}

void chi_bruteforce_GPU_CPU(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images){

double dx,dy,x_pos,y_pos;        //pixelsize
dx = (frame->xmax - frame->xmin)/(runmode->nbgridcells-1);
dy = (frame->ymax - frame->ymin)/(runmode->nbgridcells-1);
int index=0;       // index of the image within the total image array
int grid_dim;
double im_dist[MAXIMPERSOURCE]; // distance to a real image for an "theoretical" image found from a theoretical source
int im_index;       // index of the closest real image for an image found from a theoretical source
int second_closest_id;   // index of the second closest real image for an image found from a theoretical source
int thread_found_image = 0;   // at each iteration in the lens plane, turn to 1 whether the thread find an image
struct point im_position, temp; // position of the image found from a theoretical source + temp variable for comparison
struct triplet Tsup, Tinf, Tsupsource, Tinfsource;// triangles for the computation of the images created by the theoretical sources
struct galaxy sources[runmode->nsets]; // theoretical sources (common for a set)
int nimagesfound[runmode->nimagestot][MAXIMPERSOURCE]; // number of images found from the theoretical sources
struct point tim[runmode->nimagestot][MAXIMPERSOURCE]; // theoretical images (computed from sources)
int indexactual, ndim;

double *grid_gradient_x, *grid_gradient_y;

grid_gradient_x = (double *)malloc((int) (runmode->nbgridcells) * (runmode->nbgridcells) * sizeof(double));
grid_gradient_y = (double *)malloc((int) (runmode->nbgridcells) * (runmode->nbgridcells) * sizeof(double));

grid_dim = runmode->nbgridcells;
indexactual = 0;
ndim = grid_dim*grid_dim/2;
int ndivision = 2;

//Packaging the image to sourceplane conversion
gradient_grid_GPU_sub(grid_gradient_x,grid_gradient_y,frame,lens,runmode->nhalos,grid_dim,indexactual,ndim);
cudaStreamSynchronize(0);
indexactual += ndim;

for (int i = 1; i < ndivision; ++i ){
	gradient_grid_GPU_sub(grid_gradient_x,grid_gradient_y,frame,lens,runmode->nhalos,grid_dim,indexactual,ndim);
	indexactual += ndim;
}


index = 0;
*chi = 0;










free(grid_gradient_x);
free(grid_gradient_y);
}


