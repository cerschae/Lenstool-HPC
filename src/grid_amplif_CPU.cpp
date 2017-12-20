#include "grid_amplif_CPU.hpp"
#include "module_cosmodistances.hpp"
//
#ifdef __WITH_MPI
#include<mpi.h>
#endif
//
//This natrix handles the calling of the amplif functions for the different types without losing time to switch or if condition
//
typedef void (*amplif_func_t)(type_t *ampli, const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y,type_t z, int istart, int jstart);
amplif_func_t amplif_func[10] =
{
0, 0, 0, 0, 0,
amplif_5_grid_general_CPU, 0, 0, 0, 0,
};
//
void amplif_grid_CPU(type_t *amplif, const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int nhalos, int nbgridcells,int mode_amp,type_t z, int istart, int jstart)
{
	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells - 1);
	type_t dy = (frame->ymax - frame->ymin)/(nbgridcells - 1);
	//

	(*amplif_func[mode_amp])(amplif, cosmo, frame, lens, nhalos, dx, dy, nbgridcells, nbgridcells,z, istart, jstart);
}
//
//table for ampli methods
//
void amplif_5_grid_general_CPU(type_t *ampli, const struct cosmo_param *cosmo, const struct grid_param *frame, const struct Potential_SOA *lens, int Nlens, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, type_t z, int istart, int jstart)
{
	int bid = 0; // index of the block (and of the set of images)
	int tid = 0; // index of the thread within the block
	//
	type_t /*dx, dy,*/ x_pos, y_pos, kappa;        //pixelsize
	int    grid_dim;
	point   image_point;
	matrix Grad2;
	type_t dl0s = module_cosmodistances_objectObject(lens->z[0], z, *cosmo);
	//std::cerr  << lens->z[0] <<" "  << z << " "<< dl0s << std::endl;
	//
	grid_dim = nbgridcells_x*nbgridcells_y;
	//@@printf("-> dx = %f, dy = %f, istart = %d, jstart = %d, nbgridcells_x = %d, nbgridcells_y = %d\n", dx, dy, istart, jstart, nbgridcells_x, nbgridcells_y);
	//
#pragma omp parallel
#pragma omp for private(Grad2, image_point, kappa)
	for (int jj = 0; jj < nbgridcells_y; ++jj) 
		for (int ii = 0; ii < nbgridcells_x; ++ii)
		{
			int index = jj*nbgridcells_x + ii;

			image_point.x = frame->xmin + (ii + istart)*dx;
			image_point.y = frame->ymin + (jj + jstart)*dy;
			//
			Grad2 = module_potentialDerivatives_totalGradient2_SOA(&image_point, lens, Nlens);
			//
    		//Grad2.a /= dl0s;
    		//Grad2.b /= dl0s;
    		//Grad2.c /= dl0s;
    		//
            kappa = (Grad2.a + Grad2.c) / 2.;
            ampli[index] = kappa;
			//
		}
}
