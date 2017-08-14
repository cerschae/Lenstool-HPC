#include <iostream>
#include <string.h>
//#include <cuda_runtime.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
#ifndef __xlC__
#warning "gnu compilers"
#include <immintrin.h>
#endif

//#include "simd_math.h"
#include "chi_CPU.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


double myseconds()
{
        struct timeval  tp;
        struct timezone tzp;
        int i;

        i = gettimeofday(&tp,&tzp);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}




void mychi_bruteforce_SOA_CPU_grid_gradient_orig(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images){

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
//      struct triplet Tsup, Tinf, Tsupsource, Tinfsource;// triangles for the computation of the images created by the theoretical sources
        struct galaxy sources[runmode->nsets]; // theoretical sources (common for a set)
        int nimagesfound[runmode->nimagestot][MAXIMPERSOURCE]; // number of images found from the theoretical sources
        struct point tim[runmode->nimagestot][MAXIMPERSOURCE]; // theoretical images (computed from sources)


        double *grid_gradient_x, *grid_gradient_y;
	double tot_time = -myseconds();

        grid_gradient_x = (double *)malloc((int) (runmode->nbgridcells) * (runmode->nbgridcells) * sizeof(double));
        grid_gradient_y = (double *)malloc((int) (runmode->nbgridcells) * (runmode->nbgridcells) * sizeof(double));
        grid_dim = runmode->nbgridcells;
        //Packaging the image to sourceplane conversion
        double time = -myseconds();
        gradient_grid_CPU(grid_gradient_x,grid_gradient_y,frame,lens,runmode->nhalos,grid_dim);
        time += myseconds();
        printf("       Grid grad time = %f s.\n", time);

        //printf ("@@nsets = %d nx = %d ny = %d\n", runmode->nsets, runmode->nbgridcells, runmode->nbgridcells);

        index = 0;
        *chi  = 0;

        double trans_time = 0.;
        double chi2_time  = 0.;
	double inner_time = 0.;
	double p0_time    = 0.;
	double p1_time    = 0.;
	double p2_time    = 0.;
        time = -myseconds();
	//
	int images_found = 0;
	//
        for( int  source_id = 0; source_id < runmode->nsets; source_id ++)
	{

                //////////////////////////////////////Initialisation//////////////////////////////////////
                for (int i=0; i < nimages_strongLensing[source_id]; ++i){
                        for (int j=0; j < nimages_strongLensing[source_id]; ++j){
                                nimagesfound[i][j] = 0;
                        }
                }
		//
		//printf("Num sources = %d, num images = %d, Nx = %d, Ny = %d\n", runmode->nsets,  nimages_strongLensing[source_id], runmode->nbgridcells, runmode->nbgridcells);  
                for( int image_id = 0; image_id < nimages_strongLensing[source_id]; image_id++)
                {
			//printf("@@  nimages = %d\n",  nimages_strongLensing[source_id]);

                        //////////////////////////////////////computation of theoretical sources//////////////////////////////////////
                        trans_time -= myseconds();
                        mychi_transformImageToSourcePlane_SOA(runmode->nhalos, &images[index+image_id].center,images[index+image_id].dr,lens,&sources[source_id].center);
                        trans_time += myseconds();

                        //if (DEBUG ==1 )
                        //printf("index %d image_id %d source_id %d %f \n",index, image_id, source_id,images[index+image_id].redshift);
                        sources[source_id].redshift = images[index+image_id].redshift;
                        //
                        sources[source_id].dr       = images[index+image_id].dr;
                        sources[source_id].dls      = images[index+image_id].dls;
                        sources[source_id].dos      = images[index+image_id].dos;

			//inner_time -= myseconds();
                        for (int x_id = 0; x_id < runmode->nbgridcells - 1; ++x_id )
                        {
                                for (int y_id = 0; y_id < runmode->nbgridcells - 1; ++y_id )
                                {
					//@@p0_time -= myseconds();
                                        x_pos = frame->xmin + x_id * dx;
                                        y_pos = frame->ymin + y_id * dy;
                                        struct triplet Tsup, Tinf, Tsupsource, Tinfsource;
                                        // Define the upper + lower triangle, both together = square = pixel
                                        Tsup.a.x = x_pos;
                                        Tsup.b.x = x_pos;
                                        Tsup.c.x = x_pos + dx;
                                        Tinf.a.x = x_pos + dx;
                                        Tinf.b.x = x_pos + dx;
                                        Tinf.c.x = x_pos;
					//
                                        Tsup.a.y = y_pos;
                                        Tsup.b.y = y_pos + dy;
                                        Tsup.c.y = y_pos;
                                        Tinf.a.y = y_pos + dy;
                                        Tinf.b.y = y_pos;
                                        Tinf.c.y = y_pos + dy;
					//@@p0_time += myseconds();

                                        // Lens to Sourceplane conversion of triangles
                                        //@@trans_time -= myseconds();
                                        mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_upper(&Tsup, sources[source_id].dr, &Tsupsource, grid_gradient_x, grid_gradient_y, y_id*grid_dim + x_id, grid_dim);
                                        mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_lower(&Tinf, sources[source_id].dr, &Tinfsource, grid_gradient_x, grid_gradient_y, y_id*grid_dim + x_id, grid_dim);
                                        //@@trans_time += myseconds();

					//@@p1_time -= myseconds();
                                        thread_found_image=0;
					//if (x_id == 0) printf("@@ %d %d %d %d\n", source_id, y_id, mychi_insideborder(&sources[source_id].center, &Tsupsource), mychi_inside(&sources[source_id].center,&Tinfsource));
					//
                                        if(mychi_insideborder(&sources[source_id].center, &Tsupsource)==1)
                                        {
                                                thread_found_image = 1; // thread has just found an image
                                                im_index           = 0;
                                                im_position        = mychi_barycenter(&Tsup);
                                                im_dist[im_index]  = mychi_dist(im_position, images[index + im_index].center);  // get the distance to the real image
                                                for(int i=1; i<nimages_strongLensing[source_id]; i++)
                                                {  // get the distance to each real image and keep the index of the closest real image
                                                        im_dist[i] = mychi_dist(im_position,images[index+i].center);
                                                        if(im_dist[i] < im_dist[im_index])
                                                        {
                                                                im_index = i;
                                                        }
                                                }
                                        }

                                        if(mychi_inside(&sources[source_id].center,&Tinfsource)==1)
					{
                                                thread_found_image = 1; // thread has just found an image
                                                im_index           = 0;
                                                im_position        = mychi_barycenter(&Tinf);  // get the barycenter of the triangle
                                                im_dist[im_index]  = mychi_dist(im_position,images[index + im_index].center);  // get the distance to the real image
                                                for(int i=1; i<nimages_strongLensing[source_id]; i++){  // get the distance to each real image and keep the index of the closest real image

                                                        im_dist[i]=mychi_dist(im_position,images[index+i].center);
                                                        if(im_dist[i]<im_dist[im_index]){
                                                                im_index=i;
                                                        }
						//printf(" im_index %d im_dist actual %f im_dist %f \n",im_index, im_dist[im_index], im_dist[i]);
                                                }
                                        }
					//@@p1_time += myseconds();
					//	
					//continue;
                                        int skip_image = 0;
					//
					//@@p2_time -= myseconds();
                                        if (thread_found_image == 1)
					{
                                                skip_image = 0;
						images_found++;
                                                // Sometimes due to the numerical errors at the centerpoint, for SIE potentials an additional image will appear at the center of the Potential.
                                                // This is due to the fact that it is not possible to simulate an infinity value at the center correctly, Check that sis correspond to Nlens[0]
                                                for (int iterator=0; iterator < runmode->Nlens[0]; ++iterator){
                                                        if ( fabs(im_position.x - lens[0].position_x[iterator]) <= dx/2. and fabs(im_position.y  - lens[0].position_y[iterator]) <= dx/2.){
                                                                skip_image = 1;
                                                                printf("WARNING: You are using SIE potentials. An image to close to one of the potential centers has been classified as numerical error and removed \n");
                                                        }
                                                }
                                                if(!skip_image)
                                                {
                                                        //checking whether a closest image has already been found
                                                        //printf("imagenumber %d im_index %d , im_position.x %f , im_position.y %f \n", image_id, im_index  , im_position.x  , im_position.y);
							//printf("%d %d = %d\n", image_id, im_index, nimagesfound[image_id][im_index]);
                                                        if(nimagesfound[image_id][im_index]==0)
                                                        { // if no image found up to now
                                                                tim[image_id][im_index]=im_position;  //image position is allocated to theoretical image
								//printf("tim1 %d %d = %f %f\n", image_id, im_index, im_position.x, im_position.y);
                                                                nimagesfound[image_id][im_index]++;
                                                        }
                                                        else if(nimagesfound[image_id][im_index]>0)
                                                        { // if we have already found an image
                                                                // If the new image we found is closer than the previous image
                                                                if(im_dist[im_index]<mychi_dist(images[index+im_index].center,tim[image_id][im_index]))
                                                                {
                                                                        temp=tim[image_id][im_index]; // we store the position of the old image in temp
                                                                        tim[image_id][im_index]=im_position; // we link the observed image with the image we just found
									//printf("tim2 %d %d = %f %f\n", image_id, im_index, im_position.x, im_position.y);
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
                                                                // Loop over all images in the set that are not yet allocated to a theoretical image
                                                                // and allocate the closest one
                                                                for(int i=0; i<nimages_strongLensing[source_id] && nimagesfound[image_id][i]==0; i++) // we search for an observed image not already linked (nimagesfound=0)
                                                                {
                                                                        if(im_dist[i]<im_dist[second_closest_id])
                                                                        {
                                                                                second_closest_id=i;
                                                                                im_index=i; // im_index value changes only if we found a not linked yet image
                                                                                tim[image_id][im_index]=temp; // if we found an observed and not already linked image, we allocate the theoretical image temp
										//printf("tim3 %d %d = %f %f\n", image_id, im_index, temp.x, temp.y);
                                                                        }
                                                                }
                                                                nimagesfound[image_id][im_index]++; // increasing the total number of images found (If we find more than 1 theoretical image linked to 1 real image, these theoretical
                                                                // images are included in this number)
                                                        }


                                                }
                                                thread_found_image=0; // for next iteration
                                        }
					//@@p2_time += myseconds();
                                }
                        }
			//inner_time += myseconds();

                }


                //////////////////////////////////////computing the local chi square//////////////////////////////////////
                double chiimage;

                chi2_time -= myseconds();
                for( int iter = 0; iter < nimages_strongLensing[source_id]*nimages_strongLensing[source_id]; iter++){
                        int i=iter/nimages_strongLensing[source_id];
                        int j=iter % nimages_strongLensing[source_id];

                        if(i!=j){
                                // In the current method, we get the source in the source plane by ray tracing image in nimagesfound[i][i]. If we ray trace back,
                                // we arrive again at the same position and thus the chi2 from it is 0. Thus we do not calculate the chi2 (-> if i!=j)
                                if(nimagesfound[i][j]>0)
                                {
                                        chiimage = pow(images[index + j].center.x-tim[i][j].x, 2) + pow(images[index+j].center.y - tim[i][j].y, 2);  // compute the chi2
					printf("%d %d = %.15f\n", i, j, chiimage);
                                        *chi    += chiimage;
                                }
                                else if(nimagesfound[i][j]==0){
                                        // If we do not find a correpsonding image, we add a big value to the chi2 to disfavor the model
                                        *chi += 100.*nimages_strongLensing[source_id];
                                }
                        }
                }
		printf("%d: chi = %.15f\n", source_id, *chi);
                chi2_time += myseconds();
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
	//
        time += myseconds();
	//
        printf("        chi2 time = %f s.\n", time);
        printf("        - trans time = %f s.\n", trans_time);
	printf("	- inner time = %f s.\n", inner_time);
	printf("		- p0 time = %f s.\n", p0_time);
	printf("		- p1 time = %f s.\n", p1_time);
	printf("		- p2 time = %f s.\n", p2_time);
        printf("        - chi2 time = %f s.\n", chi2_time);
	//
	printf("	images found: %d\n", images_found);
        free(grid_gradient_x);
        free(grid_gradient_y);
	//
	tot_time += myseconds();
	printf("	Total time = %f\n", tot_time);  
}


void mychi_bruteforce_SOA_CPU_grid_gradient(double *chi, int *error, runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, galaxy *images){

	//double dx, dy; //, x_pos, y_pos;        //pixelsize
	const double dx = (frame->xmax - frame->xmin)/(runmode->nbgridcells-1);
	const double dy = (frame->ymax - frame->ymin)/(runmode->nbgridcells-1);
//	double im_dist[MAXIMPERSOURCE]; // distance to a real image for an "theoretical" image found from a theoretical source
//	int im_index;       		// index of the closest real image for an image found from a theoretical source
//	int second_closest_id;   	// index of the second closest real image for an image found from a theoretical source
//	int thread_found_image = 0;   	// at each iteration in the lens plane, turn to 1 whether the thread find an image
// 	struct point im_position, temp; // position of the image found from a theoretical source + temp variable for comparison
//	struct triplet Tsup, Tinf, Tsupsource, Tinfsource;// triangles for the computation of the images created by the theoretical sources
	struct galaxy sources[runmode->nsets]; // theoretical sources (common for a set)
	int nimagesfound[runmode->nimagestot][MAXIMPERSOURCE]; // number of images found from the theoretical sources
	struct point tim[runmode->nimagestot][MAXIMPERSOURCE]; // theoretical images (computed from sources)
	//
	double *grid_gradient_x, *grid_gradient_y;
	//
	grid_gradient_x = (double *)malloc((int) (runmode->nbgridcells) * (runmode->nbgridcells) * sizeof(double));
	grid_gradient_y = (double *)malloc((int) (runmode->nbgridcells) * (runmode->nbgridcells) * sizeof(double));
	//
	const int grid_dim = runmode->nbgridcells;
	//Packaging the image to sourceplane conversion
	double time = -myseconds(); 
	gradient_grid_CPU(grid_gradient_x, grid_gradient_y, frame, lens, runmode->nhalos, grid_dim);
	time += myseconds();
	printf("	Grid grad time = %f s.\n", time);
	//
	int index       = 0;       // index of the image within the total image array
	*chi  = 0;
	//
	double trans_time  = 0.;
	double trans2_time = 0.;
	double chi2_time   = 0.;
	double p1_time     = 0.;
	double p2_time	   = 0.;
	double loop_time   = 0.;
	//
	int images_found = 0;
	int images_total = 0;
        printf ("@@nsets = %d nx = %d ny = %d\n", runmode->nsets, runmode->nbgridcells, runmode->nbgridcells);
	//
	time = -myseconds();
	//
	for( int  source_id = 0; source_id < runmode->nsets; source_id ++)
	{
		unsigned short int nimages = nimages_strongLensing[source_id]; 
		//memset(&nimagesfound[0], 0, nimages_strongLensing[source_id]*nimages_strongLensing[source_id]*sizeof(nimagesfound[0][0]));
		//////////////////////////////////////Initialisation//////////////////////////////////////
		//
		for (unsigned short int i = 0; i < nimages; ++i)
		{
			for (unsigned short int j = 0; j < nimages; ++j)
		{
				nimagesfound[i][j] = 0;
			}
		}
		//
		//@@printf("nx = %d, ny = %d\n",  nimages_strongLensing[source_id], nimages_strongLensing[source_id]);
		//____________________________ image loop ________________________________
		for( unsigned short int image_id = 0; image_id < nimages; image_id++)
		{
			//printf("@@  nimages = %d\n",  nimages_strongLensing[source_id]);

			//////////////////////////////////////computation of theoretical sources//////////////////////////////////////
			double time = -myseconds();
			mychi_transformImageToSourcePlane_SOA(runmode->nhalos, &images[index+image_id].center,images[index + image_id].dr, lens, &sources[source_id].center);
			time += myseconds(); 
			//printf("%d: trans_time = %f\n", omp_get_thread_num(), time);
			trans_time += time;

			//if (DEBUG ==1 )
			sources[source_id].redshift = images[index+image_id].redshift;
			//
			sources[source_id].dr       = images[index+image_id].dr;
			sources[source_id].dls      = images[index+image_id].dls;
			sources[source_id].dos      = images[index+image_id].dos;

			loop_time -= myseconds();
#pragma omp parallel for reduction(+: trans2_time, p1_time, p2_time, images_found, images_total) 
			for (int x_id = 0; x_id < runmode->nbgridcells-1; ++x_id )
			{
				for (int y_id = 0; y_id < runmode->nbgridcells-1; ++y_id )
				{
					
					images_total++;
					double x_pos = frame->xmin + x_id*dx;
					double y_pos = frame->ymin + y_id*dy;
					//
					// Define the upper + lower triangle, both together = square = pixel
					// 
 					struct triplet Tsup, Tinf, Tsupsource, Tinfsource;
					//
					Tsup.a.x = x_pos;
					Tsup.b.x = x_pos;
					Tsup.c.x = x_pos + dx;
					Tinf.a.x = x_pos + dx;
					Tinf.b.x = x_pos + dx;
					Tinf.c.x = x_pos;
					//
					Tsup.a.y = y_pos;
					Tsup.b.y = y_pos + dy;
					Tsup.c.y = y_pos;
					Tinf.a.y = y_pos + dy;
					Tinf.b.y = y_pos;
					Tinf.c.y = y_pos + dy;
					//
					// Lens to Sourceplane conversion of triangles
					// 
					//double time = -myseconds();
					mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_upper(&Tsup, sources[source_id].dr, &Tsupsource,grid_gradient_x, grid_gradient_y, y_id*grid_dim + x_id, grid_dim);
					mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_lower(&Tinf, sources[source_id].dr, &Tinfsource,grid_gradient_x, grid_gradient_y, y_id*grid_dim + x_id, grid_dim);
					//time += myseconds();
					//printf("%d: trans_time2 = %f\n", omp_get_thread_num(), time);
					//trans2_time += time;
					//
					int    thread_found_image = 0;
					int    im_index;
					struct point im_position, temp;
					double im_dist[MAXIMPERSOURCE];
					//int    nimagesfound[runmode->nimagestot][MAXIMPERSOURCE]; // number of images found from the theoretical sources
					//
					//p1_time -= myseconds();
					if(mychi_insideborder(&sources[source_id].center, &Tsupsource) == 1)
					{
						thread_found_image = 1; // thread has just found an image
						im_index  	   = 0;
						im_position 	   = mychi_barycenter(&Tsup);
						// get the distance to the real image
						im_dist[im_index]  = mychi_dist(im_position, images[index + im_index].center);  
						for(int i = 1; i < nimages_strongLensing[source_id]; i++)
						{  
							// get the distance to each real image and keep the index of the closest real image
							im_dist[i] = mychi_dist(im_position, images[index + i].center);
							if(im_dist[i] < im_dist[im_index])
							{
								im_index = i;
							}
						}
					}
					//
					if (mychi_inside(&sources[source_id].center, &Tinfsource) == 1)
					{
						thread_found_image = 1; // thread has just found an image
						im_index           = 0;
						im_position        = mychi_barycenter(&Tinf);  // get the barycenter of the triangle
						im_dist[im_index]  = mychi_dist(im_position, images[index + im_index].center);  // get the distance to the real image
						for(int i = 1; i < nimages_strongLensing[source_id]; i++)
						{  // get the distance to each real image and keep the index of the closest real image

							im_dist[i] = mychi_dist(im_position, images[index + i].center);
							if(im_dist[i] < im_dist[im_index])
							{
								im_index = i;
							}
						}
					}
					//p1_time += myseconds();
					//
					int skip_image = 0;
					//p2_time -= myseconds();
					//
					if (thread_found_image == 1)
					{
						images_found++;
						//printf("%d %d %d: found image\n", x_id, y_id, image_id); 
						skip_image = 0;
						// Sometimes due to the numerical errors at the centerpoint, 
						// for SIE potentials an additional image will appear at the center of the Potential.
						// This is due to the fact that it is not possible to simulate an infinity value 
						// at the center correctly, Check that sis correspond to Nlens[0]
						for (int iterator=0; iterator < runmode->Nlens[0]; ++iterator)
						{
							if ( fabs(im_position.x - lens[0].position_x[iterator]) <= dx/2. and fabs(im_position.y  - lens[0].position_y[iterator]) <= dx/2.)
							{
								skip_image = 1;
								printf("WARNING: You are using SIE potentials. An image to close to one of the potential centers has been classified as numerical error and removed \n");
							}
						}
						//printf("%d %d %d %d\n", x_id, y_id,thread_found_image, skip_image);
						if (!skip_image)
						{
							//checking whether a closest image has already been found
							if(nimagesfound[image_id][im_index] == 0)
							{ // if no image found up to now
								
								//image position is allocated to theoretical image
								#pragma omp critical
								tim[image_id][im_index] = im_position;  
								#pragma omp atomic
								nimagesfound[image_id][im_index]++;
							}
							else 
							if (nimagesfound[image_id][im_index] > 0)
							{ // if we have already found an image
								// If the new image we found is closer than the previous image
								if(im_dist[im_index] < mychi_dist(images[index+im_index].center,tim[image_id][im_index]))
								{
									temp = tim[image_id][im_index]; // we store the position of the old image in temp
									#pragma omp critical
									tim[image_id][im_index]=im_position; // we link the observed image with the image we just found
									//printf("tim2 %d %d = %f %f\n", image_id, im_index, im_position.x, im_position.y);
								}
								else
								{
									temp = im_position; // we store the position of the image we just found in temp
								}
								// initialising second_closest_id to the highest value
								// Loop over all images in the set except the closest one
								// and initialize to the furthest away image
								int second_closest_id=0;
								for(int i=1; i<nimages_strongLensing[source_id] && i!=im_index; i++)
								{
									if(im_dist[i]>im_dist[second_closest_id]) second_closest_id=i;
								}
								///////////////////////////////////////////////////////////////
								// Loop over all images in the set that are not yet allocated to a theoretical image
								// and allocate the closest one
								for(int i=0; i < nimages_strongLensing[source_id] && nimagesfound[image_id][i]==0; i++) // we search for an observed image not already linked (nimagesfound=0)
								{
									if(im_dist[i]<im_dist[second_closest_id])
									{
										second_closest_id=i;
										im_index=i; // im_index value changes only if we found a not linked yet image
										//printf("tim3 %d %d = %f %f\n", image_id, im_index, temp.x, temp.y);
										#pragma omp critical
										tim[image_id][im_index] = temp; // if we found an observed and not already linked image, we allocate the theoretical image temp
									}
								}
								#pragma omp atomic
								nimagesfound[image_id][im_index]++; // increasing the total number of images found (If we find more than 1 theoretical image linked to 1 real image, these theoretical
								// images are included in this number)
							}
						}
						thread_found_image=0; // for next iteration
					}
					//p2_time += myseconds();
				}
			}
//#pragma omp barrier
			loop_time += myseconds();

		}
		//____________________________ end of image loop
		//
		//____________________________ computing the local chi square
		//
		double chiimage;
		//
		chi2_time -= myseconds();
		int _nimages = nimages_strongLensing[source_id];
		//
		for( int iter = 0; iter < _nimages*_nimages; iter++)
		{
			int i=iter/nimages_strongLensing[source_id];
			int j=iter%nimages_strongLensing[source_id];

			//printf("nimagesfound %d %d = %d\n", i, j, nimagesfound[i][j]);
			if(i != j)
			{
				// In the current method, we get the source in the source plane by ray tracing image in nimagesfound[i][i]. If we ray trace back,
				// we arrive again at the same position and thus the chi2 from it is 0. Thus we do not calculate the chi2 (-> if i!=j)
				if(nimagesfound[i][j] > 0)
				{
					chiimage = pow(images[index + j].center.x - tim[i][j].x, 2) + pow(images[index+j].center.y - tim[i][j].y, 2);  // compute the chi2
					//printf("%d %d = %.15f\n", i, j, chiimage);
					*chi    += chiimage;
				}
				else if(nimagesfound[i][j]==0)
				{
					// If we do not find a correpsonding image, we add a big value to the chi2 to disfavor the model
					*chi += 100.*nimages_strongLensing[source_id];
				}
			}
		}
		//
		//____________________________ end of computing the local chi square
		//
		//printf("%d: chi = %.15f\n", source_id, *chi);
		chi2_time += myseconds();
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
	time += myseconds();
	int nthreads = 1;
#pragma omp parallel
	nthreads = omp_get_num_threads();

	printf("	chi2 time = %f s. using %d threads\n", time, nthreads);
	printf("	- trans  time         = %f s.\n", trans_time);
	printf("	- loop   time 	      = %f s.\n", loop_time/nthreads);
	printf("		- trans2 time = %f s.\n", trans2_time/nthreads);
	printf("		- p1     time = %f s.\n", p1_time/nthreads);
	printf("		- p2     time = %f s.\n", p2_time/nthreads);
	printf("	- chi2 update    time = %f s.\n", chi2_time);
	
	printf("	images found: %d out of %d\n", images_found, images_total);

	free(grid_gradient_x);
	free(grid_gradient_y);
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
void mychi_transformImageToSourcePlane(const runmode_param *runmode, const struct point *image_point, double dlsds, const struct Potential *lens, struct point *source_point)
{   // dlsds is the distance between lens and source divided by the distance observer-source
	struct point Grad;  // gradient

	Grad = module_potentialDerivatives_totalGradient(runmode->nhalos, image_point, lens);
	//Grad = module_potentialDerivatives_totalGradient_SOA(image_point, lens, runmode->Nlens);

	source_point->x = image_point->x - dlsds*Grad.x;
	source_point->y = image_point->y - dlsds*Grad.y;
	//printf("dlsds %f", dlsds);
}


void mychi_transformImageToSourcePlane_SOA(const int Nlens, const struct point *image_point, double dlsds, const struct Potential_SOA *lens, struct point *source_point)
{   // dlsds is the distance between lens and source divided by the distance observer-source
	struct point Grad;  // gradient
	Grad = module_potentialDerivatives_totalGradient_SOA(image_point, lens, Nlens);

	source_point->x = image_point->x - dlsds*Grad.x;
	source_point->y = image_point->y - dlsds*Grad.y;
	//printf("dlsds %f", dlsds);
}

inline
void mychi_transformImageToSourcePlane_SOA_Packed( const struct point *image_point, double dlsds, struct point *source_point, double *grad_x, double * grad_y, int grad_id)
{

	source_point->x = image_point->x - dlsds*grad_x[grad_id];
	source_point->y = image_point->y - dlsds*grad_y[grad_id];
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

inline
void mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_upper( struct triplet *I, double dlsds, struct triplet *S, double *grad_x, double * grad_y, int grad_id, int nbgridcell)
{
	mychi_transformImageToSourcePlane_SOA_Packed( &I->a, dlsds,   &S->a, grad_x, grad_y, grad_id             );
	mychi_transformImageToSourcePlane_SOA_Packed( &I->b, dlsds,   &S->b, grad_x, grad_y, grad_id + nbgridcell);
	mychi_transformImageToSourcePlane_SOA_Packed( &I->c, dlsds,   &S->c, grad_x, grad_y, grad_id + 1);
}

inline
void mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_lower( struct triplet *I, double dlsds, struct triplet *S, double *grad_x, double * grad_y, int grad_id, int nbgridcell)
{
	mychi_transformImageToSourcePlane_SOA_Packed( &I->a, dlsds,   &S->a, grad_x, grad_y, grad_id + nbgridcell + 1);
	mychi_transformImageToSourcePlane_SOA_Packed( &I->b, dlsds,   &S->b, grad_x, grad_y, grad_id + 1    );
	mychi_transformImageToSourcePlane_SOA_Packed( &I->c, dlsds,   &S->c, grad_x, grad_y, grad_id + nbgridcell             );
}



/** @brief Return the scalar triple product (a*b).c of the 3 vectors A[x,y,1], B[x,y,1], C[x,y,1].
 * If 2 of the 3 vectors are equal, colinear or form an orthogonal basis,
 * the triple product is 0.
 * This is also the determinant of the matrix
 *   | Ax  Bx  Cx |
 *   | Ay  By  Cy |
 *   |  1   1   1 |
 */
//inline
double mychi_determinant(const struct point *A,
		const struct point *B,
		const struct point *C)
{
	return( B->x*C->y - B->y*C->x +
			A->x*B->y - A->y*B->x +
			A->y*C->x - A->x*C->y );
}


/** @brief Return 1 if P is inside the triangle T, 0 otherwise.
 *
 * Return 1 if P is inside the triangle T, 0 otherwise.
 * @param P  a point
 * @param T  a triplet of points.
 *
 *
 */
inline
int mychi_inside(const struct point *P, struct triplet *T)
{
	double  s, s1, s2, d;

	d  = mychi_determinant(&T->a, &T->b, &T->c);
	s  = mychi_determinant(&T->a, &T->b, P)*d;
	s1 = mychi_determinant(&T->b, &T->c, P)*d;
	s2 = mychi_determinant(&T->c, &T->a, P)*d;

	return((s > 0.) && (s1 > 0.) && (s2 > 0.));  // If all determinants are positive,
	// the point must be inside the triangle
}


/*
int
mychi_inside2(const struct point *A, const struct point *B, const struct point *C)
{

	// Compute vectors        
	v0 = C - A;
	v1 = B - A;
	v2 = P - A;

	// Compute dot products
	dot00 = dot(v0, v0);
	dot01 = dot(v0, v1);
	dot02 = dot(v0, v2);
	dot11 = dot(v1, v1);
	dot12 = dot(v1, v2);

	// Compute barycentric coordinates
	invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point is in triangle
	return (u >= 0) && (v >= 0) && (u + v < 1);
}
*/

/** @brief Return 1 if P is inside the triangle T or on its border, 0 otherwise.
 *
 * Return 1 if P is inside the triangle T or on its border, 0 otherwise.
 * @param  P  a point
 * @param  T  a triplet of points.
 *
 *
 */
	inline
int mychi_insideborder(const struct point *P, struct triplet *T)
{
	double  s, s1, s2, d;

	d  = mychi_determinant(&T->a, &T->b, &T->c);
	s  = mychi_determinant(&T->a, &T->b, P)*d;
	s1 = mychi_determinant(&T->b, &T->c, P)*d;
	s2 = mychi_determinant(&T->c, &T->a, P)*d;

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
struct  point   mychi_barycenter(struct triplet *A)
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
double mychi_dist(struct point A, struct point B)
{
	double  x, y;
	x = A.x - B.x;
	y = A.y - B.y;
	return(sqrt(x*x + y*y));
}

