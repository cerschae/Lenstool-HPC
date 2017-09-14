#include "grid_gradient_CPU.hpp"
//
#ifdef __WITH_MPI
#include<mpi.h>
#endif
//
//
//
void gradient_grid_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells)
{
	double dx = (frame->xmax - frame->xmin)/(nbgridcells - 1);
        double dy = (frame->ymax - frame->ymin)/(nbgridcells - 1);
	//
	gradient_grid_general_CPU(grid_grad_x, grid_grad_y, frame, lens, nhalos, dx, dy, nbgridcells, nbgridcells, 0, 0);
}
//
//
//
void gradient_grid_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
	double dx = (frame->xmax - frame->xmin)/(nbgridcells_x - 1);
        double dy = (frame->ymax - frame->ymin)/(nbgridcells_y - 1);
	//
	printf("istart = %d, jstart = %d, nbgridcells_x = %d, ngridcells_y = %d\n", istart, jstart, nbgridcells_x, nbgridcells_y);
	gradient_grid_general_CPU(grid_grad_x, grid_grad_y, frame, lens, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
}
//
//
//
void gradient_grid_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, double dx, double dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
        //@@printf("dx = %f, dy = %f, istart = %d, jstart = %d, nbgridcells_x = %d, ngridcells_y = %d\n", dx, dy, istart, jstart, nbgridcells_x, nbgridcells_y);
	gradient_grid_general_CPU(grid_grad_x, grid_grad_y, frame, lens, nhalos, dx, dy, nbgridcells_x, nbgridcells_y, istart, jstart);
}
//
//
//
void gradient_grid_general_CPU_old(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcells, const struct Potential_SOA *lens)
{
	int bid=0; // index of the block (and of the set of images)
	int tid=0; // index of the thread within the block

	double dx,dy,x_pos,y_pos;        //pixelsize
	int grid_dim, index;
	point Grad, image_point, true_coord_rotation;
	double      R;
	dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
	dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
	grid_dim = (nbgridcells);

	index = bid ;
	//std::cout << "xmin = " << frame->xmin << " ymin = " << frame->ymin << std::endl; 

	while(index < grid_dim*grid_dim)
	{
		grid_grad_x[index] = 0.;
		grid_grad_y[index] = 0.;

		image_point.x = frame->xmin + (index%grid_dim)*dx;
		image_point.y = frame->ymin + (index/grid_dim)*dy;

		Grad = module_potentialDerivatives_totalGradient_SOA(&image_point, lens, Nlens);
		//
		//std::cout << "col = " << index%grid_dim << " row = " << index/grid_dim << " " << image_point.x << " " << image_point.y << ": " << Grad.x << " " << Grad.y << std::endl;

		grid_grad_x[index] = Grad.x;
		grid_grad_y[index] = Grad.y;
		//return;

		bid += 1;
		index = bid * 1 + tid;
	}
}


void gradient_grid_general_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, double dx, double dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
	int bid = 0; // index of the block (and of the set of images)
	int tid = 0; // index of the thread within the block
	//
	double /*dx, dy,*/ x_pos, y_pos;        //pixelsize
	int    grid_dim;
	point  Grad, image_point, true_coord_rotation;
	//
	//dx = (frame->xmax - frame->xmin)/(nbgridcells_x - 1);
	//dy = (frame->ymax - frame->ymin)/(nbgridcells_y - 1);
	//
	grid_dim = nbgridcells_x*nbgridcells_y;
	//@@printf("-> dx = %f, dy = %f, istart = %d, jstart = %d, nbgridcells_x = %d, nbgridcells_y = %d\n", dx, dy, istart, jstart, nbgridcells_x, nbgridcells_y);
	//
#pragma omp parallel
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
			/*
#if __WITH_MPI
			//
			int world_rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
			int world_size;
			MPI_Comm_rank(MPI_COMM_WORLD, &world_size);
			//
			if (world_rank == 0) printf("%d x = %f y = %f -> %f %f\n", index, image_point.x, image_point.y, Grad.x, Grad.y); 
#endif
			*/
			grid_grad_x[index] = Grad.x;
			grid_grad_y[index] = Grad.y;
			//
		}
}


/*
void gradient_grid_general_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, int nbgridcells_x, int nbgridc
ells_y, const struct Potential_SOA *lens, int istart, int jstart)
{
        int bid = 0; // index of the block (and of the set of images)
        int tid = 0; // index of the thread within the block
        //
        double dx, dy, x_pos, y_pos;        //pixelsize
        int    grid_dim;
        point  Grad, image_point, true_coord_rotation;
        //
        dx = (frame->xmax - frame->xmin)/(nbgridcells_x - 1);
        dy = (frame->ymax - frame->ymin)/(nbgridcells_y - 1);
        //
        grid_dim = nbgridcells_x*nbgridcells_y;
        //
#pragma omp parallel
#pragma omp for private(Grad, image_point)
        for (int jj = jstart; jj < nbgridcells_y; ++jj)
                for (int ii = istart; ii < nbgridcells_x; ++ii)
                {
                        //  (index < grid_dim*grid_dim)
                        int index = jj*nbgridcells_x + ii;
                        //grid_grad_x[index] = 0.;
                        //grid_grad_y[index] = 0.;
                        //
                        image_point.x = frame->xmin + ii*dx;
                        image_point.y = frame->ymin + jj*dy;
                        //
                        Grad = module_potentialDerivatives_totalGradient_SOA(&image_point, lens, Nlens);
                        //
                        grid_grad_x[index] = Grad.x;
                        grid_grad_y[index] = Grad.y;
                        //
                        //std::cout << image_point.x << " " << image_point.y << ": " << Grad.x << " " << Grad.y << std::endl;
                }
}
*/




