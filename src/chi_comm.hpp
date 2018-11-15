#include <mpi.h>
#include <string.h>
#include "structure_hpc.hpp"
#include "mpi_check.h"
//
void
delense_comm(int* numimagesfound, struct point* imagesposition, int* numimg, const int *nimages_strongLensing, const int* locimagesfound, struct point* image_pos, const int nsets, const int world_rank, const int world_size);
