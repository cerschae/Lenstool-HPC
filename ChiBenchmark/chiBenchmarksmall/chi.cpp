#include <iostream>
#include <string.h>
//#include <cuda_runtime.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <immintrin.h>

#include "simd_math.h"
#include "chi.hpp"
#include "structure.h"

void chi_bruteforce(double *chi, int *error, runmode_param *runmode, const struct Potential *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images){

	int bid=0;
	int tid=0; 	// index of the thread within the block
	double dx,dy;				//pixelsize
	dx = (frame->xmax - frame->xmin)/(runmode->nbgridcells-1);
	dy = (frame->ymax - frame->ymin)/(runmode->nbgridcells-1);
	int imagenumber; 		// index of the image within a set
	int index=0; 			// index of the image within the total image array
	double im_dist[MAXIMPERSOURCE]; // distance to a real image for an "theoretical" image found from a theoretical source
	int im_index; 			// index of the closest real image for an image found from a theoretical source
	int second_closest_id; 	// index of the second closest real image for an image found from a theoretical source
	int thread_found_image = 0; 	// at each iteration in the lens plane, turn to 1 whether the thread find an image
	struct point im_position, temp; // position of the image found from a theoretical source + temp variable for comparison
	struct triplet Tsup, Tinf, Tsupsource, Tinfsource;// triangles for the computation of the images created by the theoretical sources
	int Nb_thread_row = (int) sqrtf(threadsPerBlock); // Nb_thread_row is the number of threads per line
	struct galaxy sources[MAXIMPERSOURCE]; // theoretical sources (common for a set)
	int nimagesfound[MAXIMPERSOURCE][MAXIMPERSOURCE]; // number of images found from the theoretical sources
	struct point tim[MAXIMPERSOURCE][MAXIMPERSOURCE]; // theoretical images (computed from sources)
	double temp_chi; // temporary chi square result (for a set of images)
	int mylock; // atomic lock variable for cuda atomic operations, used to avoid interference between threads when they write the images positions
	struct triplet Ttemp; //Temp triangle needed for rescaling information
	int tooManyImages;	//Overflow warning, If too many images overflow the memory allocated to the stacks, the code stops
	int DEBUG = 0;
	if(bid ==0)
		DEBUG = 1;

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
			module_chiClassic_transformImageToSourcePlane(runmode, &images[index+image_id].center,images[index+image_id].dr,lens,&sources[source_id].center);
			//if (DEBUG ==1 )
				//printf("index %d image_id %d source_id %d %f \n",index, image_id, source_id,images[index+image_id].redshift);
			sources[source_id].redshift = images[index+image_id].redshift;
			sources[source_id].dr = images[index+image_id].dr;
			sources[source_id].dls = images[index+image_id].dls;
			sources[source_id].dos = images[index+image_id].dos;

			for (double x_pos = frame->xmin; x_pos <= frame->xmax; x_pos += dx){
				for (double y_pos = frame->ymin; y_pos <= frame->ymax; y_pos += dx){
					//printf("%f %f %f %f %f %f \n", frame->xmin, frame->ymin, x_pos, y_pos,frame->xmax,frame->ymax);

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
					module_chiClassic_transformtriangleImageToSourcePlane(runmode,&Tsup,sources[source_id].dr, lens, &Tsupsource);
					module_chiClassic_transformtriangleImageToSourcePlane(runmode,&Tinf,sources[source_id].dr, lens, &Tinfsource);

					thread_found_image=0;
					if(module_chiClassic_insideborder(&sources[source_id].center,&Tsupsource)==1){
						thread_found_image=1; // thread has just found an image
						im_index=0;
						im_position=module_chiClassic_barycenter(&Tsup);
						im_dist[im_index]=module_chiClassic_dist(im_position,images[index+im_index].center);  // get the distance to the real image
						for(int i=1; i<nimages_strongLensing[source_id]; i++){  // get the distance to each real image and keep the index of the closest real image
							im_dist[i]=module_chiClassic_dist(im_position,images[index+i].center);
							if(im_dist[i]<im_dist[im_index]){
								im_index=i;
							}
							//printf(" im_index %d im_dist actual %f im_dist %f \n",im_index, im_dist[im_index], im_dist[i]);
						}
						if(DEBUG ==1){
							//printf(" %d Found something at  x %f y %f by bid %d and tid %d\n", image_id,x_pos,y_pos, bid , tid);
							//printf("T ax %f ay %f cx %f cy %f Tsour ax %f ay %f cx %f cy %fby bid %d and tid %d\n", Tinf.a.x, Tinf.a.y , Tinf.c.x, Tinf.c.y ,Tinfsource.a.x, Tinfsource.a.y , Tinfsource.c.x, Tinfsource.c.y , bid , tid);
							//printf(" %d Source x %f y %f by bid %d and tid %d\n", source_id,sources[source_id].center.x, sources[source_id].center.y , bid , tid);
							printf(" %d Img x %f y %f  im_dist %f and im_index %d\n",image_id, im_position.x, im_position.y , im_dist[im_index] , im_index);
						}
					}

					if(module_chiClassic_inside(&sources[source_id].center,&Tinfsource)==1){
						thread_found_image=1; // thread has just found an image
						im_index=0;
						im_position=module_chiClassic_barycenter(&Tinf);  // get the barycenter of the triangle
						im_dist[im_index]=module_chiClassic_dist(im_position,images[index+im_index].center);  // get the distance to the real image
						for(int i=1; i<nimages_strongLensing[source_id]; i++){  // get the distance to each real image and keep the index of the closest real image

							im_dist[i]=module_chiClassic_dist(im_position,images[index+i].center);
							//printf("im_dist[i] %f, im_position %f %f , images[index+im_index].center %f %f\n",im_dist[i], im_position.x,im_position.y, images[index+i].center.x,images[index+i].center.y);
							if(im_dist[i]<im_dist[im_index]){
								im_index=i;
							}
						}
						if(DEBUG ==1){
							//printf(" %d Found something at  x %f y %f by bid %d and tid %d\n",image_id, x_pos,y_pos, bid , tid);
							//printf("T ax %f ay %f cx %f cy %f Tsour ax %f ay %f cx %f cy %fby bid %d and tid %d\n", Tinf.a.x, Tinf.a.y , Tinf.c.x, Tinf.c.y ,Tinfsource.a.x, Tinfsource.a.y , Tinfsource.c.x, Tinfsource.c.y , bid , tid);
							//printf(" %d Source x %f y %f by bid %d and tid %d\n",source_id, sources[source_id].center.x, sources[source_id].center.y , bid , tid);
							printf(" %d Img x %f y %f  im_dist %f and im_index %d\n",image_id, im_position.x, im_position.y , im_dist[im_index] , im_index);
						}
					}

					int skip_image = 0;

					if(thread_found_image == 1 ){
						skip_image = 0;
						// Sometimes due to the numerical errors at the centerpoint, for SIE potentials an additional image will appear at the center of the Potential.
						// This is due to the fact that it is not possible to simulate an infinity value at the center correctly
						for (int i=0; i < runmode->nhalos; ++i){
							//printf("lens[i].type %d %f %f %f \n",lens[i].type, fabs(im_position.x - lens[i].position.x) ,  fabs(im_position.y - lens[i].position.y), dx/2.);
							if ( lens[i].type == 5 and fabs(im_position.x - lens[i].position.x) <= dx/2. and fabs(im_position.y  - lens[i].position.y) <= dx/2.){
								skip_image = 1;
								printf("WARNING: You are using SIE potentials. An image to close to one of the potential centers has been classified as numerical error and removed \n");
							}
						}
						if(skip_image==0){
									//checking whether a closest image has already been found
							//printf("imagenumber %d im_index %d , im_position.x %f , im_position.y %f \n", image_id, im_index  , im_position.x  , im_position.y);
							if(nimagesfound[image_id][im_index]==0){ // if no image found up to now
								tim[image_id][im_index]=im_position;  //image position is allocated to theoretical image
								nimagesfound[image_id][im_index]++;
							}
							else if(nimagesfound[image_id][im_index]>0){ // if we have already found an image
								// If the new image we found is closer than the previous image
								if(im_dist[im_index]<module_chiClassic_dist(images[index+im_index].center,tim[image_id][im_index]))
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

		for (int i=0; i < nimages_strongLensing[source_id]; ++i){
			for (int j=0; j < nimages_strongLensing[source_id]; ++j){
				printf(" %d",nimagesfound[i][j]);
				}
			printf("\n");
			}

		//Incrementing Index: Images already treated by previous source_id
		index+=nimages_strongLensing[source_id];
	}
}




/** @brief Tranform a point from image to source plane. Result stored in sourcepoint argument
 *
 * Tranform a point from image to source plane using lensequation
 *
 * @param image_point	 image position
 * @param dlsds		 dls/ds
 * @param nhalos	 number of halos
 * @param potential_param	 gravitational potential information
 * @param source_point	 address where source information will be stored
 *
*
*/
void module_chiClassic_transformImageToSourcePlane(const runmode_param *runmode, const struct point *image_point, double dlsds, const struct Potential *lens, struct point *source_point)
{   // dlsds is the distance between lens and source divided by the distance observer-source
    struct point Grad;  // gradient

    Grad = module_potentialDerivatives_totalGradient(runmode,image_point, lens);

    source_point->x = image_point->x - dlsds * Grad.x;
    source_point->y = image_point->y - dlsds * Grad.y;
    //printf("dlsds %f", dlsds);
}


/** @brief Tranform a triangle from image to source plane. Result stored in S triangle argument
 *
 * Return a triplet of points in the source plane corresponding to the triplet
 * of images. dlsds is the lens efficiency at the source redshift.
 * I is the triangle in the image plane (input), S is the same triangle in the source plane (output)
 *
 * @param I	 triangle in image plane
 * @param dlsds	 dls/ds
 * @param nhalos	 number of halos
 * @param potential_param	 gravitational potential information
 * @param S	 address where triangle source information will be stored
 *
*
*/
void module_chiClassic_transformtriangleImageToSourcePlane(const runmode_param *runmode, struct triplet *I, double dlsds, const struct Potential *lens, struct triplet *S)
{
  module_chiClassic_transformImageToSourcePlane(runmode, &I->a, dlsds,  lens, &S->a);
  module_chiClassic_transformImageToSourcePlane(runmode, &I->b, dlsds,  lens, &S->b);
  module_chiClassic_transformImageToSourcePlane(runmode, &I->c, dlsds,  lens, &S->c);
}

/** @brief Return the scalar triple product (a*b).c of the 3 vectors A[x,y,1], B[x,y,1], C[x,y,1].
 * If 2 of the 3 vectors are equal, colinear or form an orthogonal basis,
 * the triple product is 0.
 * This is also the determinant of the matrix
 *   | Ax  Bx  Cx |
 *   | Ay  By  Cy |
 *   |  1   1   1 |
 */
double module_chiClassic_determinant(const struct point *A,
                   const struct point *B,
                   const struct point *C)
{
    return( B->x * C->y - B->y * C->x +
            A->x * B->y - A->y * B->x +
            A->y * C->x - A->x * C->y );
}


/** @brief Return 1 if P is inside the triangle T, 0 otherwise.
*
* Return 1 if P is inside the triangle T, 0 otherwise.
* @param P  a point
* @param T  a triplet of points.
*
*
*/
int module_chiClassic_inside(const struct point *P, struct triplet *T)
{
  double  s, s1, s2, d;

  d = module_chiClassic_determinant(&T->a, &T->b, &T->c);
  s = module_chiClassic_determinant(&T->a, &T->b, P) * d;
  s1 = module_chiClassic_determinant(&T->b, &T->c, P) * d;
  s2 = module_chiClassic_determinant(&T->c, &T->a, P) * d;

  return((s > 0.) && (s1 > 0.) && (s2 > 0.));  // If all determinants are positive,
                                               // the point must be inside the triangle
}



/** @brief Return 1 if P is inside the triangle T or on its border, 0 otherwise.
*
* Return 1 if P is inside the triangle T or on its border, 0 otherwise.
* @param  P  a point
* @param  T  a triplet of points.
*
*
*/
int module_chiClassic_insideborder(const struct point *P, struct triplet *T)
{
  double  s, s1, s2, d;

  d = module_chiClassic_determinant(&T->a, &T->b, &T->c);
  s = module_chiClassic_determinant(&T->a, &T->b, P) * d;
  s1 = module_chiClassic_determinant(&T->b, &T->c, P) * d;
  s2 = module_chiClassic_determinant(&T->c, &T->a, P) * d;

  return((s >= 0.) && (s1 >= 0.) && (s2 >= 0.));  // If all determinants are positive or 0,
                                               // the point must be inside the triangle or on its border

}

/** @brief Barycentre of a triplet/triangle
 *
 * A is a structure triplet that contains 3 structures point a,b and c
 * Return value B is a point
 *
*
*/
struct  point   module_chiClassic_barycenter(struct triplet *A)
{
    struct  point   B;

    B.x = (A->a.x + A->b.x + A->c.x) / 3.;
    B.y = (A->a.y + A->b.y + A->c.y) / 3.;
    return(B);
}

/** @brief Euclidean distance between 2 points
*
* Euclidean distance between 2 points
*
*/
double module_chiClassic_dist(struct point A, struct point B)
{
   double  x, y;
   x = A.x - B.x;
   y = A.y - B.y;
   return(sqrt(x*x + y*y));
}

