/**
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)
* @date   July 2017
* @version 0,1
*
*/

#include "grid_gradient_CPU.hpp"
//
#ifdef __WITH_MPI
#include<mpi.h>
#endif
//
//
//
void gradient_grid_general_CPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart);

void gradient_grid_CPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells, int istart, int jstart)
{
	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells - 1);
	type_t dy = (frame->ymax - frame->ymin)/(nbgridcells - 1);
	//
	gradient_grid_general_CPU(grid_grad_x, grid_grad_y, frame, lens, nhalos, dx, dy, nbgridcells, nbgridcells, 0, 0);
}

void gradient_grid_CPU_print(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells, int istart, int jstart)
{
	gradient_grid_print_CPU(grid_grad_x, grid_grad_y, frame, nhalos, nbgridcells, lens, istart, jstart);
}


void gradient_grid_CPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells_x - 1);
	type_t dy = (frame->ymax - frame->ymin)/(nbgridcells_y - 1);
	//
	printf("istart = %d, jstart = %d, nbgridcells_x = %d, ngridcells_y = %d\n", istart, jstart, nbgridcells_x, nbgridcells_y);
	gradient_grid_general_CPU(grid_grad_x, grid_grad_y, frame, lens, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
}
//
//
//
void gradient_grid_CPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
        //@@printf("dx = %f, dy = %f, istart = %d, jstart = %d, nbgridcells_x = %d, ngridcells_y = %d\n", dx, dy, istart, jstart, nbgridcells_x, nbgridcells_y);
	gradient_grid_general_CPU(grid_grad_x, grid_grad_y, frame, lens, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
}
//
//
//
void gradient_grid_general_CPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
	int bid = 0; // index of the block (and of the set of images)
	int tid = 0; // index of the thread within the block
	//
	type_t /*dx, dy,*/ x_pos, y_pos;        //pixelsize
	int    grid_dim;
	point  Grad, image_point, true_coord_rotation;
	//
	grid_dim = nbgridcells_x*nbgridcells_y;
	//
//#pragma omp parallel
#pragma omp for private(Grad, image_point)
	for (int jj = 0; jj < nbgridcells_y; ++jj)
		for (int ii = 0; ii < nbgridcells_x; ++ii)
		{
			//  (index < grid_dim*grid_dim)

			int index = jj*nbgridcells_x + ii;
			//grid_grad_x[index] = 0.;
			//grid_grad_y[index] = 0.;
			//
			image_point.x = frame->xmin + (ii + istart)*dx;
			image_point.y = frame->ymin + (jj + jstart)*dy;

			//
			Grad = module_potentialDerivatives_totalGradient_SOA(&image_point, lens, Nlens);
			//

			grid_grad_x[index] = Grad.x;
			grid_grad_y[index] = Grad.y;
			//if(ii == 0 and jj == 0)
				//std::cout << grid_grad_x[index] << image_point.x <<image_point.y << std::endl;

			//
		}
	//std::cout << std::endl;
}

void gradient_grid_print_CPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, int Nlens, int nbgridcells, const struct Potential_SOA *lens, int istart, int jstart)
{
	int bid = 0; // index of the block (and of the set of images)
	int tid = 0; // index of the thread within the block
	//
	type_t dx,dy,x_pos,y_pos;        //pixelsize
	int    grid_dim, index;
	point  Grad, image_point, true_coord_rotation;
	type_t R;
	//
	dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
	dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
	grid_dim = (nbgridcells);
	//
	index = bid ;
	//
	for (int jj = jstart; jj < nbgridcells; ++jj)
		for (int ii = istart; ii < nbgridcells; ++ii)
		{
			//  (index < grid_dim*grid_dim)

			int index = jj*nbgridcells + ii;
			//grid_grad_x[index] = 0.;
			//grid_grad_y[index] = 0.;

			image_point.x = frame->xmin + ii*dx;
			image_point.y = frame->ymin + jj*dy;

			Grad = module_potentialDerivatives_totalGradient_5_SOA_print(&image_point, lens, 0, 1, index);
			//
			grid_grad_x[index] = Grad.x;
			grid_grad_y[index] = Grad.y;
			//
			//std::cout << image_point.x << " " << image_point.y << ": " << Grad.x << " " << Grad.y << std::endl;

			//bid += 1;
			//index = bid * 1 + tid;
		}
}
