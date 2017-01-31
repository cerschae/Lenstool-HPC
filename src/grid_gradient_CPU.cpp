#include "grid_gradient_CPU.hpp"



void gradient_grid_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells){

  gradient_grid_general_CPU(grid_grad_x, grid_grad_y,frame,nhalos,nbgridcells, lens);

}

void gradient_grid_general_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens,int nbgridcells, const struct Potential_SOA *lens){
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

  while(index < grid_dim*grid_dim){

    grid_grad_x[index] = 0.;
    grid_grad_y[index] = 0.;

    image_point.x = frame->xmin + (index/grid_dim)*dx;
    image_point.y = frame->ymin + (index % grid_dim)*dy;

    Grad = module_potentialDerivatives_totalGradient_SOA(&image_point, lens, Nlens);

    grid_grad_x[index] = Grad.x;
    grid_grad_y[index] = Grad.y;

    bid += 1;
    index = bid * 1 + tid;
  }
}
