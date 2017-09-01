#include "grid_gradient_CPU.hpp"
//
//
//
void gradient_grid_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells, int istart, int jstart)
{
	gradient_grid_general_CPU(grid_grad_x, grid_grad_y, frame, nhalos, nbgridcells, lens, istart, jstart);
}

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


void gradient_grid_general_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, int nbgridcells, const struct Potential_SOA *lens, int istart, int jstart)
{
	int bid = 0; // index of the block (and of the set of images)
	int tid = 0; // index of the thread within the block
	//
	double dx,dy,x_pos,y_pos;        //pixelsize
	int    grid_dim, index;
	point  Grad, image_point, true_coord_rotation;
	double R;
	//
	dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
	dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
	grid_dim = (nbgridcells);
	//
	index = bid ;
	//
#pragma omp parallel
#pragma omp for private(Grad, image_point)
	for (int jj = jstart; jj < nbgridcells; ++jj) 
		for (int ii = istart; ii < nbgridcells; ++ii)
		{
			//  (index < grid_dim*grid_dim)

			int index = jj*nbgridcells + ii;
			//grid_grad_x[index] = 0.;
			//grid_grad_y[index] = 0.;

			image_point.x = frame->xmin + ii*dx;
			image_point.y = frame->ymin + jj*dy;

			Grad = module_potentialDerivatives_totalGradient_SOA(&image_point, lens, Nlens);
			//
			grid_grad_x[index] = Grad.x;
			grid_grad_y[index] = Grad.y;
			//
			//std::cout << image_point.x << " " << image_point.y << ": " << Grad.x << " " << Grad.y << std::endl;

			//bid += 1;
			//index = bid * 1 + tid;
		}
}
