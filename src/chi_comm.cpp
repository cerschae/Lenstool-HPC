//
#include <string.h>
#include "structure_hpc.hpp"
#ifdef __WITH_MPI
#include <mpi.h>
#include "mpi_check.h"
#endif

//
void
delense_comm(int* numimagesfound, struct point* imagesposition, int* numimg, const int *nimages_strongLensing, const int* locimagesfound, struct point* image_pos, const int nsets, const int world_rank, const int world_size)
{
        //
        int          numimagesfound_tmp[nsets];
        struct point imagesposition_tmp[nsets][MAXIMPERSOURCE];
        //
        memset(numimagesfound_tmp, 0, nsets*sizeof(int));
        memset(imagesposition_tmp, 0, nsets*sizeof(int));
	//
	int verbose = (world_rank == 0);
        //
        int total = 0;
        //
        if (!verbose)
        {
                MPI_CHECK(MPI_Send( numimg, 1, MPI_INT   , 0, 666 + world_rank, MPI_COMM_WORLD ));
                if (*numimg != 0)
                {
                        MPI_CHECK(MPI_Send( locimagesfound, nsets                 , MPI_INT   , 0, 666 + world_rank, MPI_COMM_WORLD ));
                        MPI_CHECK(MPI_Send( image_pos,      nsets*MAXIMPERSOURCE*2, MPI_DOUBLE, 0, 667 + world_rank, MPI_COMM_WORLD ));
                }
        }
	else
	{
		int image_sum = 0;
		//
		for (int ipe = 0; ipe < world_size; ++ipe)
		{
			MPI_Status status;
			//
			if (ipe == 0)
			{
				memcpy(&numimagesfound_tmp, locimagesfound, nsets*sizeof(int));
				memcpy(&imagesposition_tmp, image_pos,      nsets*MAXIMPERSOURCE*sizeof(point));
			}
			else
			{
				MPI_CHECK(MPI_Recv(numimg, 1, MPI_INT, ipe, 666 + ipe, MPI_COMM_WORLD, &status));
				if (*numimg != 0)
				{
					MPI_CHECK(MPI_Recv(&numimagesfound_tmp, nsets                 , MPI_INT   , ipe, 666 + ipe, MPI_COMM_WORLD, &status));
					MPI_CHECK(MPI_Recv(&imagesposition_tmp, nsets*MAXIMPERSOURCE*2, MPI_DOUBLE, ipe, 667 + ipe, MPI_COMM_WORLD, &status));
				}
			}
			//
			//MPI_Reduce(&imagesfound_tmp, &total, ipe, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			//
			if (*numimg != 0)
			{
				//for (int jj = 0; jj < runmode->nimagestot; ++jj)
				//{
				int loc_numimg   = 0;
				int loc_numimg_s = 0;
				for (int ii = 0; ii < nsets; ++ii)
				{
					//int img_len = numimagesfound[ii][jj];
					int img_len = numimagesfound_tmp[ii];
					image_sum  += img_len;
					if (img_len != 0)
					{
						int loc_length = numimagesfound[ii];
						memcpy(&imagesposition[ii*MAXIMPERSOURCE + loc_length], &imagesposition_tmp[ii], img_len*sizeof(point));
						numimagesfound[ii] += img_len;
						loc_numimg   += img_len;
						loc_numimg_s += nimages_strongLensing[ii];
						//printf("%d: %d , img_len = %d, source = %d, sum = %d %d\n", ipe, ii, numimagesfound_tmp[ii], nimages_strongLensing[ii], loc_numimg, loc_numimg_s);

					}
				}
			}
			printf("        rank %d: Total number of images found = %d\n", ipe, *numimg);
		}
	}
}
