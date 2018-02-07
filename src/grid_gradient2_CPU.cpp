#include "grid_gradient2_CPU.hpp"
//
#ifdef __WITH_MPI
#include<mpi.h>
#endif
//
//

void gradient2_grid_CPU(matrix *grid_grad2, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells, int istart, int jstart)
{
	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells - 1);
	type_t dy = (frame->ymax - frame->ymin)/(nbgridcells - 1);
	//
	gradient2_grid_general_CPU(grid_grad2, frame, lens, nhalos, dx, dy, nbgridcells, nbgridcells, istart, jstart);
}

void gradient2_grid_CPU(matrix *grid_grad2,  const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells_x - 1);
	type_t dy = (frame->ymax - frame->ymin)/(nbgridcells_y - 1);
	//
	printf("istart = %d, jstart = %d, nbgridcells_x = %d, ngridcells_y = %d\n", istart, jstart, nbgridcells_x, nbgridcells_y);
	gradient2_grid_general_CPU(grid_grad2, frame, lens, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
}
//
//
//
void gradient2_grid_CPU(matrix *grid_grad2, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, double dx, double dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
	gradient2_grid_general_CPU(grid_grad2, frame, lens, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
}
//
//
//
void gradient2_grid_general_CPU(matrix *grid_grad2, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
	int bid = 0; // index of the block (and of the set of images)
	int tid = 0; // index of the thread within the block
	//
	type_t /*dx, dy,*/ x_pos, y_pos;        //pixelsize
	int    grid_dim;
	point   image_point;
	matrix Grad;
	//
	grid_dim = nbgridcells_x*nbgridcells_y;
	//
//#pragma omp parallel
//#pragma omp for private(Grad, image_point)
	for (int jj = 0; jj < nbgridcells_y; ++jj) 
		for (int ii = 0; ii < nbgridcells_x; ++ii)
		{
			int index = jj*nbgridcells_x + ii;

			image_point.x = frame->xmin + (ii + istart)*dx;
			image_point.y = frame->ymin + (jj + jstart)*dy;
			//
			Grad = module_potentialDerivatives_totalGradient2_SOA(&image_point, lens, Nlens);
			//
			grid_grad2[index] = Grad;
			//
		}
}
