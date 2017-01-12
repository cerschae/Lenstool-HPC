#include "grid_srcplane_conversion.hpp"
#include <gradient_avx.hpp>


void grid_srcplane_conversion_CPU(double *grid_srcplane_x, double *grid_srcplane_y, const struct grid_param *frame, const struct Potential_SOA *lens, const double dlsds, int *Nlens, int nbgridcells){

  //SIS
  grid_srcplane_conversion_sis_CPU(grid_srcplane_x, grid_srcplane_y,frame,Nlens[0],nbgridcells,dlsds, &lens[0]);
  //PIEMD8
  grid_srcplane_conversion_piemd_CPU(grid_srcplane_x, grid_srcplane_y,frame,Nlens[1],nbgridcells,dlsds, &lens[1]);
  //PIEMD81


}

void grid_srcplane_conversion_sis_CPU(double *grid_srcplane_x, double *grid_srcplane_y, const struct grid_param *frame, int Nlens, int nbgridcells,const double dlsds, const struct Potential_SOA *lens){
  int bid=0; // index of the block (and of the set of images)
  int tid=0; // index of the thread within the block

  double dx,dy,x_pos,y_pos;        //pixelsize
  point Grad, image_point;
  int grid_dim, index;
  struct point true_coord, true_coord_rotation;
  double      R;
  dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
  dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
  grid_dim = (nbgridcells);

  index = bid ;

  while(index < grid_dim*grid_dim){
    grid_srcplane_x[index] = 0.;
    grid_srcplane_y[index] = 0.;

    image_point.x = frame->xmin + (index/grid_dim)*dx;
    image_point.y = frame->ymin + (index % grid_dim)*dy;


    Grad = module_potentialDerivatives_totalGradient_5_SOA(&image_point, lens,0, Nlens);

    grid_srcplane_x[index] = image_point.x - dlsds * Grad.x;  //position in srcplane
    grid_srcplane_y[index] = image_point.y - dlsds * Grad.y;

    bid += 1;
    index = bid * 1 + tid;
  }
}


//void grid_srcplane_conversion_piemd_CPU(double *grid_srcplane_x, double *grid_srcplane_y, const struct grid_param *frame, int Nlens, int nbgridcells,const double dlsds, double * lens_x, double *lens_y, double * lens_b0, double *lens_angle, double *lens_epot, double *lens_rcore, double *lens_rcut){
void grid_srcplane_conversion_piemd_CPU(double *grid_srcplane_x, double *grid_srcplane_y, const struct grid_param *frame, int Nlens, int nbgridcells,const double dlsds, const struct Potential_SOA *lens){
  int bid=0; // index of the block (and of the set of images)
  int tid=0; // index of the thread within the block

  double dx,dy,x_pos,y_pos;        //pixelsize
  point Grad, image_point;
  int grid_dim, index;
  struct point true_coord, true_coord_rotation;
  complex      zis;
  dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
  dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
  grid_dim = (nbgridcells);

  index = bid ;



  //while(index < grid_dim*grid_dim){
  while(index < grid_dim*grid_dim){

    image_point.x = frame->xmin + (index/grid_dim)*dx;
    image_point.y = frame->ymin + (index % grid_dim)*dy;


    Grad = module_potentialDerivatives_totalGradient_8_SOA(&image_point, lens,0, Nlens);

    grid_srcplane_x[index] = grid_srcplane_x[index] - dlsds * Grad.x;  //position in srcplane
    grid_srcplane_y[index] = grid_srcplane_y[index] - dlsds * Grad.y;

    bid += 1;
    index = bid * 1 + tid;
  }
}
