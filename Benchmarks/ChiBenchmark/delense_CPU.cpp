// This file is part of lenstoolHPC
// authors: gilles.fourestey@epfl.ch

#include "delense_CPU.hpp"
#if 1
void delense_barycenter(struct point* image_pos, int* locimagesfound, int* numimg, const runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, const struct galaxy* sources, double* grid_gradient_x, double* grid_gradient_y)
{
#define INDEX2D_BAR(y, x)         (MAXIMPERSOURCE*y + x)
	//const unsigned int nimagestot  = runmode->nimagestot;
	const unsigned int nsets       = runmode->nsets;
	const unsigned int nbgridcells = runmode->nbgridcells;
	int world_size = 1;
	int world_rank = 0;
#ifdef __WITH_MPI
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif
	unsigned int verbose = (world_rank == 0);
	//
	const int grid_dim      = runmode->nbgridcells;
	const int grid_size     = nbgridcells;
	const int loc_grid_size = nbgridcells/world_size;
	//
	double y_len      = fabs(frame->ymax - frame->ymin);
	int    y_len_loc  = nbgridcells/world_size;
	int    y_pos_loc  = (int) world_rank*y_len_loc;
	int    y_bound    = y_len_loc;
	//
	if ((world_rank + 1) != world_size) y_bound++;
	//
	const double dx   = (frame->xmax - frame->xmin)/(nbgridcells - 1);
	const double dy   = (frame->ymax - frame->ymin)/(nbgridcells - 1);
	//
	int images_total  = 0;
	int index         = 0;
	//
	int lif[nsets][16];
	memset(&lif,0,  nsets*16*sizeof(int));	
	//
		//locimagesfound[source_id] = 0;
		//lif[source_id] = 0;
		// number of images in the image plane for the specific image (1,3,5...)
		//unsigned short int nimages = nimages_strongLensing[source_id];
		//printf("@@ source_id = %d, nimages = %d\n",  source_id, nimages_strongLensing[source_id]);
		//____________________________ image (constrains) loop ________________________________
		//
		//struct point image_pos [MAXIMPERSOURCE];
		//
		//MPI_Barrier(MPI_COMM_WORLD);
		//if (verbose) printf("source = %d, image = %d\n", source_id, image_id);
		//if (verbose) fflush(stdout);
	for( int  source_id = 0; source_id < nsets; source_id ++)
	{
		int loc_images_found = 0;
		
		//                      MPI_Barrier(MPI_COMM_WORLD);
		//#endif
#pragma omp parallel
#pragma omp for reduction(+: images_total) 
		for (int y_id = 0; y_id < (y_bound - 1); ++y_id)
		{
			for (int x_id = 0; x_id < runmode->nbgridcells - 1 ; ++x_id)
			{
				images_total++;
				//
				double x_pos = frame->xmin + (            x_id)*dx;
				double y_pos = frame->ymin + (y_pos_loc + y_id)*dy;
				//
				// Define the upper + lower triangle, both together = square = pixel
				//
				struct triplet Tsup, Tinf;
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
				// double time = -myseconds();
				//
				struct triplet Timage;
				struct triplet Tsource;
				//
				mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_upper(&Tsup, sources[source_id].dr, &Tsource, grid_gradient_x, grid_gradient_y, y_id*grid_dim + x_id, grid_dim);
				//
				int thread_found_image = 0;
				//
				if (mychi_insideborder(&sources[source_id].center, &Tsource) == 1)
				{
					thread_found_image = 1;
					Timage = Tsup;
				}
				else
				{
					mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_lower(&Tinf, sources[source_id].dr, &Tsource, grid_gradient_x, grid_gradient_y, y_id*grid_dim + x_id, grid_dim);
					if (mychi_inside(&sources[source_id].center, &Tsource) == 1)
					{
						thread_found_image = 1;
						Timage = Tinf;
					}
				}
#if 1
				if (thread_found_image)
				{
#pragma omp critical
					{
						// get the barycenter of the triangle
						image_pos[INDEX2D_BAR(source_id, loc_images_found)] = mychi_barycenter(&Timage);
                                                //locimagesfound[source_id]++;
                                                loc_images_found++;
                                                *numimg = *numimg + 1;
						//locimagesfound[source_id]++;
						lif[source_id][0]++;
						//loc_images_found++;
						//*numimg = *numimg + 1;
					}
				}
#endif
			}
		}
		index += nimages_strongLensing[source_id];


	}
#pragma omp parallel for
	for (int ii = 0; ii < nsets; ++ii)
		locimagesfound[ii] = lif[ii][0];
}
#endif

void delense(struct point* image_pos, int* locimagesfound, int* numimg, const runmode_param *runmode, const struct Potential_SOA *lens, const struct grid_param *frame, const int *nimages_strongLensing, const struct galaxy* sources, double* grid_gradient_x, double* grid_gradient_y)
{
	const unsigned int nsets       = runmode->nsets;
	const unsigned int nimagestot  = runmode->nimagestot;
#define INDEX2D(y, x)         (nimagestot*y + x)
#define INDEX3D(y, x, z) (nimagestot*MAXIMPERSOURCE*y + MAXIMPERSOURCE*x + z)

	//const unsigned int nsets       = runmode->nsets;
	//const unsigned int nimagestot  = runmode->nimagestot;
	const unsigned int nbgridcells = runmode->nbgridcells;
	//printf("nsets = %d, nimagestot = %d\n", nsets, nimagestot);
	int world_size = 1;
	int world_rank = 0;
#ifdef __WITH_MPI
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	unsigned int verbose = (world_rank == 0);
	//
	const int grid_dim 	= runmode->nbgridcells;
	const int grid_size     = nbgridcells;
	const int loc_grid_size = nbgridcells/world_size;
	//
	double y_len      = fabs(frame->ymax - frame->ymin);
	int    y_len_loc  = nbgridcells/world_size;
	int    y_pos_loc  = (int) world_rank*y_len_loc;
	int    y_bound    = y_len_loc;
	//
	if ((world_rank + 1) != world_size) y_bound++;
	//
	const double dx   = (frame->xmax - frame->xmin)/(nbgridcells - 1);
	const double dy   = (frame->ymax - frame->ymin)/(nbgridcells - 1);
	//
	int images_total  = 0;
	int index         = 0;
	//
	for( int source_id = 0; source_id < nsets; source_id ++)
	{
		// number of images in the image plane for the specific image (1,3,5...)
		const unsigned short int nimages = nimages_strongLensing[source_id];
		//____________________________ image (constrains) loop ________________________________
		for(unsigned short int image_id = 0; image_id < nimages; image_id++)
		{
			//
			//struct point image_pos [MAXIMPERSOURCE];
			//
			int loc_images_found = 0;
			//
#pragma omp parallel
#pragma omp for reduction(+: images_total)
			for (int y_id = 0; y_id < (y_bound - 1); ++y_id)
			{
				for (int x_id = 0; x_id < runmode->nbgridcells - 1 ; ++x_id)
				{
					//int yy_pos = MIN(y_pos_loc + y_id, runmode->nbgridcells - 1);
					images_total++;
					//
					double x_pos = frame->xmin + (            x_id)*dx;
					double y_pos = frame->ymin + (y_pos_loc + y_id)*dy;
					//
					// Define the upper + lower triangle, both together = square = pixel
					//
					struct triplet Tsup, Tinf;
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
					//
					Tsup.c.y = y_pos;
					Tinf.a.y = y_pos + dy;
					Tinf.b.y = y_pos;
					Tinf.c.y = y_pos + dy;
					//
					// Lens to Sourceplane conversion of triangles
					//
					struct triplet Timage;
					struct triplet Tsource;
					//
					double dr = sources[INDEX2D(source_id, image_id)].dr;
                                        mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_upper(&Tsup, dr, &Tsource, grid_gradient_x, grid_gradient_y, y_id*grid_dim + x_id, grid_dim);
                                        //
                                        int thread_found_image = 0;
                                        //
                                        if (mychi_insideborder(&sources[INDEX2D(source_id, image_id)].center, &Tsource) == 1)
                                        {
                                                thread_found_image = 1;
                                                Timage = Tsup;
                                        }
                                        else
                                        {
                                                mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_lower(&Tinf, sources[INDEX2D(source_id, image_id)].dr, &Tsource, grid_gradient_x, grid_gradient_y, y_id*grid_dim + x_id, grid_dim);
                                                if (mychi_inside(&sources[INDEX2D(source_id, image_id)].center, &Tsource) == 1)
                                                {
                                                        thread_found_image = 1;
                                                        Timage = Tinf;
                                                }
                                        }
#if 1
                                        if (thread_found_image)
                                        {
#pragma omp critical
                                                {
							// get the barycenter of the triangle
                                                        image_pos[INDEX3D(source_id, image_id, loc_images_found)] = mychi_barycenter(&Timage);
                                                        locimagesfound[INDEX2D(source_id, image_id)]++;
                                                        loc_images_found++;
                                                        *numimg = *numimg + 1;
                                                }
                                        }
#endif
                                }
                        }
//#if 1
                }
                index += nimages_strongLensing[source_id];
       } 
}
