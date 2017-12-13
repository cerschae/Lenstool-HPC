#include <stdlib.h>
#include "grid_gradient_mixed_CPU.hpp"
//
//
//
#define SIS_THRESHOLD 0.0092
//
//
struct point_double module_potentialDerivatives_totalGradient_5_SOA_DP_v2(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos)
{
        asm volatile("# module_potentialDerivatives_totalGradient_SIS_SOA_v2 begins");
        //printf("# module_potentialDerivatives_totalGradient_SIS_SOA_v2 begins\n");
        //
        struct point_double grad, result;
        grad.x = 0;
        grad.y = 0;
        for(int i = shalos; i < shalos + nhalos; i++)
        {
                //
                struct point true_coord, true_coord_rotation;
                //
                true_coord.x = pImage->x - lens->position_x[i];
                true_coord.y = pImage->y - lens->position_y[i];
                //
                //true_coord_rotation = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[i]);
                double cose = lens->anglecos[i];
                double sine = lens->anglesin[i];
                //
                double x = true_coord.x*cose + true_coord.y*sine;
                double y = true_coord.y*cose - true_coord.x*sine;
                //
                double ell_pot = lens->ellipticity_potential[i];
                //
                double R = sqrt(x*x*(1 - ell_pot) + y*y*(1 + ell_pot));
                //
                result.x = (1 - ell_pot)*lens->b0[i]*x/R;
                result.y = (1 + ell_pot)*lens->b0[i]*y/R;
                //
                grad.x += result.x*cose - result.y*sine;
                grad.y += result.y*cose + result.x*sine;
        }
        return grad;
}


void gradient_grid_general_double_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, int nbgridcells, const struct Potential_SOA *lens, int istart, int jstart)
{
        int bid = 0; // index of the block (and of the set of images)
        int tid = 0; // index of the thread within the block
        //
        type_t dx,dy,x_pos,y_pos;        //pixelsize
        int    grid_dim;
        point_double  Grad_DP;
        point  Grad, image_point, true_coord_rotation;
        type_t R;
        //
        dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
        dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
        grid_dim = (nbgridcells);
        //
        //
#pragma omp parallel
#pragma omp for private(Grad, image_point)
        for (int jj = jstart; jj < nbgridcells; ++jj)
                for (int ii = istart; ii < nbgridcells; ++ii)
                {
                        int index = jj*nbgridcells + ii;
                        //
                        image_point.x = frame->xmin + ii*dx;
                        image_point.y = frame->ymin + jj*dy;

                        //Grad = module_potentialDerivatives_totalGradient_SOA(&image_point, lens, 0, Nlens);
                        Grad_DP = module_potentialDerivatives_totalGradient_5_SOA_DP_v2(&image_point, lens, 0, Nlens);
                        //
                        grid_grad_x[index] = (double) Grad_DP.x;
                        grid_grad_y[index] = (double) Grad_DP.y;
                        //
                }
}


//
//
/*
void gradient_grid_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells, int istart, int jstart)
{
	gradient_grid_general_CPU(grid_grad_x, grid_grad_y, frame, nhalos, nbgridcells, lens, istart, jstart);
}
*/
//
//
void gradient_grid_general_mixed_CPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, int Nlens, int nbgridcells, const struct Potential_SOA *lens, int istart, int jstart)
{
	int bid = 0; // index of the block (and of the set of images)
	int tid = 0; // index of the thread within the block
	//
	type_t dx,dy,x_pos,y_pos;        //pixelsize
	int    grid_dim;
	point_double  Grad_DP; 
	point  Grad, image_point, true_coord_rotation;
	type_t R;
	//
	dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
	dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
	grid_dim = (nbgridcells);
	//
	//
	double* grid_grad_temp_x = (double*) malloc(nbgridcells*nbgridcells*sizeof(double));
	double* grid_grad_temp_y = (double*) malloc(nbgridcells*nbgridcells*sizeof(double));

	//
	//
#pragma omp parallel
#pragma omp for private(Grad, image_point)
	for (int jj = jstart; jj < nbgridcells; ++jj) 
		for (int ii = istart; ii < nbgridcells; ++ii)
		{
			int index = jj*nbgridcells + ii;
			//
			image_point.x = frame->xmin + ii*dx;
			image_point.y = frame->ymin + jj*dy;

			//Grad = module_potentialDerivatives_totalGradient_SOA(&image_point, lens, 0, Nlens);
			Grad = module_potentialDerivatives_totalGradient_SOA(&image_point, lens, Nlens);
			//
			grid_grad_temp_x[index] = (double) Grad.x;
			grid_grad_temp_y[index] = (double) Grad.y;
			//
		}

	memcpy(grid_grad_x, grid_grad_temp_x, nbgridcells*nbgridcells*sizeof(double));
	memcpy(grid_grad_y, grid_grad_temp_y, nbgridcells*nbgridcells*sizeof(double));

#if 1
#pragma omp parallel
#pragma omp for private(Grad, image_point)
        for (int jj = jstart + 1; jj < nbgridcells; ++jj)
                for (int ii = istart + 1; ii < nbgridcells; ++ii)
                {
			int index       = (jj - 0)*nbgridcells + (ii - 0);
                        int index_north = (jj - 1)*nbgridcells + (ii - 0);
                        int index_west  = (jj - 0)*nbgridcells + (ii - 1);
                        //
			type_t grad_north_x = grid_grad_temp_x[index] - grid_grad_temp_x[index_north];
			type_t grad_north_y = grid_grad_temp_y[index] - grid_grad_temp_y[index_north];
			//
			type_t grad_west_x  = grid_grad_temp_x[index] - grid_grad_temp_x[index_west];
			type_t grad_west_y  = grid_grad_temp_y[index] - grid_grad_temp_y[index_west];
			//
			int c1 = fabs(dx - grad_north_x) < SIS_THRESHOLD;
			int c2 = fabs(dx -  grad_west_x) < SIS_THRESHOLD; 
			int c3 = fabs(dy - grad_north_y) < SIS_THRESHOLD;
			int c4 = fabs(dy -  grad_west_y) < SIS_THRESHOLD; 
			//
			if (c1 || c2 || c3 || c4)
			{
				struct point pImage;
				pImage.x = frame->xmin + ii*dx; 
                                pImage.y = frame->ymin + jj*dx; 
				//
				Grad_DP = module_potentialDerivatives_totalGradient_5_SOA_DP_v2(&pImage, lens, 0, Nlens);
				//
				grid_grad_x[ii*nbgridcells + jj] = Grad_DP.x;
                                grid_grad_y[ii*nbgridcells + jj] = Grad_DP.y;

			}
			//
                        //Grad = module_potentialDerivatives_totalGradient_SOA(&image_point, lens, 0, Nlens);
                        //Grad = module_potentialDerivatives_totalGradient_SOA(&image_point, lens, Nlens);
                        //
                        //grid_grad_x[index] = (double) Grad_DP.x;
                        //grid_grad_y[index] = (double) Grad_DP.y;
                        //
                }

#endif

#if 1
	for (int pot = 0; pot < Nlens; ++pot)
	{
		struct point center;
		int patch_size = 20;
		if (lens->vdisp[pot] > 900.)
		{
			int patch_size = 200;
		}
		struct point pImage, lens_int;
		for (int ii = 0; ii < patch_size; ++ii)
			for (int jj = 0; jj < patch_size; ++jj)
			{
				int px = (int) (lens->position_x[pot] - frame->xmin)/dx; 
				int py = (int) (lens->position_y[pot] - frame->ymin)/dy; 
				//
				//pImage.x = lens->position_x[pot] - (patch_size/2 + ii)*dx;
				//pImage.y = lens->position_y[pot] - (patch_size/2 + jj)*dy;
				int jc = px - patch_size/2 + ii;
                                int ic = py - patch_size/2 + jj; 
				//
				pImage.x = jc*dx + frame->xmin; 
				pImage.y = ic*dy + frame->ymin; 
				//	
				Grad_DP = module_potentialDerivatives_totalGradient_5_SOA_DP_v2(&pImage, lens, 0, Nlens);
				//
				//printf("%d %d\n", lens_int.x - 
				//int jc = (pImage.x - frame->xmin)/dx;
                                //int ic = (pImage.y - frame->ymin)/dy;	
				// needs checking
				//
				std::cout << "*** " << ic*nbgridcells + jc << " " << pImage.x << " " << pImage.y << " ic = " << ic << " jc = " << jc << " Grad.x = " << Grad_DP.x << " Grad.y = " << Grad_DP.y << std::endl;
				grid_grad_x[ic*nbgridcells + jc] = Grad_DP.x;  
				grid_grad_y[ic*nbgridcells + jc] = Grad_DP.y;  
			}

	}
#endif
}
