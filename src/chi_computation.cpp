//
#include <string.h>
#include "structure_hpc.hpp"
#include "delense_CPU_utils.hpp"
#ifdef __WITH_MPI
#include <mpi.h>
#include "mpi_check.h"
#endif

void
chi_computation(double* chi, int* images_found, int* numimagesfound, struct point* imagesposition, const int *nimages_strongLensing, const galaxy* images, const int nsets)
{
	int    nimagesfound  [nsets][MAXIMPERSOURCE]; // number of images found from the theoretical sources
	struct point tim     [MAXIMPERSOURCE]; // theoretical images (computed from sources)
	int    index = 0;
	//
	for( int  source_id = 0; source_id < nsets; source_id ++)
	{
		unsigned short int nimages = nimages_strongLensing[source_id];
		//
		//______________________Initialisation______________________
		//
		//for (unsigned short int i = 0; i < nimages; ++i)
		//	nimagesfound[source_id][i] = 0;
		memset(&nimagesfound[source_id][0], 0, nimages*sizeof(int));  
		//
		//for( unsigned short int image_id = 0; image_id < nimages; image_id++)
		//{
		//#endif //collapsing
		//
		double image_dist[MAXIMPERSOURCE];
		//
		//printf("      Image %d, number of sources found %d\n", image_id, loc_images_found);
		//
		struct point image_position;
		int          image_index;
		//
		for (int ii = 0; ii < /*loc*/numimagesfound[source_id]; ++ii)
		{
			//
			//
			int image_index = 0;
			//
			image_position = imagesposition[source_id*MAXIMPERSOURCE + ii];
			image_dist[0] = mychi_dist(image_position, images[index + 0].center);  // get the distance to the real image
			for(int i = 1; i < nimages_strongLensing[source_id]; i++)
			{  // get the distance to each real image and keep the index of the closest real image

				image_dist[i] = mychi_dist(image_position, images[index + i].center);
				if (image_dist[i] < image_dist[image_index])
				{
					image_index = i;
				}
			}
			//
			// we should exit loops here
			//
			// p1_time += myseconds();
			//
			int skip_image = 0;
			// Sometimes due to the numerical errors at the centerpoint,
			// for SIE potentials an additional image will appear at the center of the Potential.
			// This is due to the fact that it is not possible to simulate an infinity value
			// at the center correctly, Check that sis correspond to Nlens[0]
			/*
			   for (int iterator = 0; iterator < runmode->Nlens[0]; ++iterator)
			   {
			   if ( fabs(image_position.x - lens[0].position_x[iterator]) <= dx/2. and fabs(image_position.y  - lens[0].position_y[iterator]) <= dx/2.)
			   {
			   skip_image = 1;
			   printf("WARNING: You are using SIE potentials. An image to close to one of the potential centers has been classified as numerical error and removed \n");
			   }
			   }*/
			//printf("%d %d %d %d\n", x_id, y_id,thread_found_image, skip_image);
			if (!skip_image)
			{
				//#pragma omp atomic
				(*images_found)++;
				struct point temp;
				//printf("                        source %d, image %d, index %d, Images found: %d\n", source_id, image_id, image_index, nimagesfound[source_id][image_id][image_index]);
				//checking whether a closest image has already been found
				if (nimagesfound[source_id][image_index] == 0)
				{ // if no image found up to now

					//image position is allocated to theoretical image
					//#pragma omp critical
					tim[image_index] = image_position;
					//#pragma omp atomic
					nimagesfound[source_id][image_index]++;
				}
				else if (nimagesfound[source_id][image_index] > 0)
				{ // if we have already found an image
					// If the new image we found is closer than the previous image
					//printf("      tim2: %f %f\n", image_dist[image_index], mychi_dist(images[index + image_index].center, tim[image_id][image_index]));
					if (image_dist[image_index] < mychi_dist(images[index + image_index].center, tim[image_index]))
					{
						temp = tim[image_index]; // we store the position of the old image in temp
						//#pragma omp critical
						tim[image_index] = image_position; // we link the observed image with the image we just found
						//printf("tim2 %d %d = %f %f\n", image_id, image_index, image_position.x, image_position.y);
					}
					else
					{
						temp = image_position; // we store the position of the image we just found in temp
					}
					// initialising second_closest_id to the highest value
					// Loop over all images in the set except the closest one
					// and initialize to the furthest away image
					int second_closest_id = 0;
					for (int i = 1; i < nimages_strongLensing[source_id] && (i != image_index); i++)
					{
						if(image_dist[i] > image_dist[second_closest_id]) second_closest_id=i;
					}
					///////////////////////////////////////////////////////////////
					// Loop over all images in the set that are not yet allocated to a theoretical image
					// and allocate the closest one
					// we search for an observed image not already linked (nimagesfound=0)
					for(int i = 0; i < nimages_strongLensing[source_id] && nimagesfound[source_id][i] == 0; i++)
					{
						if(image_dist[i] < image_dist[second_closest_id])
						{
							second_closest_id = i;
							// im_index value changes only if we found a not linked yet image
							image_index = i;
							//printf("tim3 %d %d = %f %f\n", image_id, image_index, temp.x, temp.y);
							//#pragma omp critical
							tim[image_index] = temp; // if we found an observed and not already linked image, we allocate the theoretical image temp
						}
					}
					// increasing the total number of images found (If we find more than 1 theoretical image linked to 1 real image,
					// these theoretical
					//#pragma omp atomic
					nimagesfound[source_id][image_index]++;
					// images are included in this number)
				}
			}
			//#pragma omp atomic
			//loc_images_found++;
			//thread_found_image  = 0; // for next iteration
		}
		//
		//}
		//#pragma omp barrier
		//____________________________ end of image loop
		//
		//____________________________ computing the local chi square
		//
		double chiimage;
		//
		int _nimages = nimages_strongLensing[source_id];
		//
		for( int iter = 0; iter < _nimages; iter++)
		{
			int i = iter;
			//int i=iter/nimages_strongLensing[source_id];
			//int j=iter%nimages_strongLensing[source_id];

			//printf("nimagesfound %d %d = %d\n", i, j, nimagesfound[i][j]);
			//if(i != j)
			//{
			if(nimagesfound[source_id][i] > 0)
			{
				double pow1 = images[index + i].center.x - tim[i].x;
				double pow2 = images[index + i].center.y - tim[i].y;
				//
				//chiimage = pow(images[index + j].center.x - tim[i][j].x, 2) + pow(images[index + j].center.y - tim[i][j].y, 2);  // compute the chi2
				chiimage = pow1*pow1 + pow2*pow2;  // compute the chi2
				//printf("%d %d = %.15f\n", i, j, chiimage);
				(*chi)    += chiimage;
			}
			else
				if(nimagesfound[source_id][i] == 0)
				{
					// If we do not find a correpsonding image, we add a big value to the chi2 to disfavor the model
					(*chi) += 1000.*nimages_strongLensing[source_id];
					printf("        source_id = %d, image number = %d has not found an image, chi = %f\n", source_id, i, *chi);
				}
			//}
		}
		//
		//____________________________ end of computing the local chi square
		//
		//printf("%d: chi = %.15f\n", source_id, *chi);
		/*
		   for (int i=0; i < nimages_strongLensing[source_id]; ++i){
		   for (int j=0; j < nimages_strongLensing[source_id]; ++j){
		   printf(" %d",nimagesfound[i][j]);
		   }
		   printf("\n");
		   }*/

		//Incrementing Index: Images already treated by previous source_id
		index += nimages_strongLensing[source_id];
	}
}
