/**
 * @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)
 * @date   July 2017
 * @version 0,1
 *
 */
#include <fstream>
#include "potential_GPU.cuh"
#include "structure_hpc.hpp"

#define BLOCK_SIZE_X 8
#define BLOCK_SIZE_Y 16

//#define ROT

#define _SHARED_MEM

#ifdef _SHARED_MEM
#define SHARED __shared__
#warning "shared memory"
extern __shared__ type_t shared[];
#else
#define SHARED 
#endif

#define Nx 1
#define Ny 0

extern "C" {
double myseconds();
}

__device__ double  pi05_gpu(double x, double y, double eps, double rc, double b0)
{
    double  sqe, ci, cxro, cyro, rem2, e1, e2, z;
    complex eta, zeta, b1, b2, a1, a2, c1, c2, ckk;


    sqe = sqrt(eps);
    ci = .5 * (1. - eps * eps) / sqe;
    cxro = (1. + eps) * (1. + eps);
    cyro = (1. - eps) * (1. - eps);
    rem2 = x * x / cxro + y * y / cyro;
    e1 = 2.*sqe / (1 - eps);
    e2 = 2.*sqe / (1 + eps);
    z = sqrt(x * x + y * y);
    eta = cpx_GPU(-.5 * asinh(e1 * y / z), .5 * asin(e2 * x / z));
    zeta = cpx_GPU( 0.5 * log( (sqrt(rem2) + sqrt(rc * rc + rem2)) / rc), 0. );
    b1 = coshcpx_GPU(acpx_GPU(eta, zeta));
    b2 = coshcpx_GPU(scpx_GPU(eta, zeta));
    a1 = lncpx_GPU( dcpx_GPU(sqcpx_GPU(coshcpx_GPU(eta)), pcpx_GPU(b1, b2)) );
    a2 = lncpx_GPU(dcpx_GPU(b1, b2));
    c1 = pcpx_GPU(sinhcpx_GPU(pcpxflt_GPU(eta, 2.)), a1);
    c2 = pcpx_GPU(sinhcpx_GPU(pcpxflt_GPU(zeta, 2.)), a2);
    ckk = acpx_GPU(c1, c2);

    return( b0*ci*rc / sqrt(rem2)*(ckk.im*x - ckk.re*y) );
}


__device__ type_t module_potential_81_SOA_GPU(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos){
	//asm volatile("# module_potentialDerivatives_totalGradient_81_SOA begins");
#if 0
	int col = blockIdx.x*blockDim.x + threadIdx.x;
	int row = blockIdx.y*blockDim.y + threadIdx.y;
	if(row == 0 && col == 0) printf("# module_potentialDerivatives_totalGradient_81_SOA begins\n");
#endif
	//std::cout << "# module_potentialDerivatives_totalGradient_81_SOA begins" << std::endl;
	// 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
	//
	type_t t05;
	type_t potential, pa, ps;
	potential =  0;

for(int i = shalos; i < shalos + nhalos; i++)
{

	struct point true_coord;
	//True coord
	true_coord.x = pImage->x - lens->position_x[i];
	true_coord.y = pImage->y - lens->position_y[i];
	//Rotation
	type_t cose = lens->anglecos[i];
	type_t sine = lens->anglesin[i];
	//
	type_t x = true_coord.x*cose + true_coord.y*sine;
	type_t y = true_coord.y*cose - true_coord.x*sine;
	// 81 comput
	t05 = lens->rcut[i] / (lens->rcut[i] - lens->rcore[i]);

	//if(blockIdx.x*blockDim.x + threadIdx.x == 0 && blockIdx.y*blockDim.y + threadIdx.y == 0)printf("t05 %f rcut: %f, rcore: %f, b0:%f \n",t05, lens->rcut[i],lens->rcore[i], lens->b0[i]);
	pa = pi05_gpu(x, y, lens->ellipticity_potential[i], lens->rcore[i], lens->b0[i]);
	////////////////////////////////////////////////
	ps = pi05_gpu(x, y, lens->ellipticity_potential[i], lens->rcut[i], lens->b0[i]);
	/////////////////////////////////////
	potential += t05 * (pa - ps);
	///////////
}

return(potential);

}

#if 1
typedef type_t (*potential_func_GPU_t) (const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);

__constant__ potential_func_GPU_t potential_func_GPU[100] =
{
		0, 0, 0, 0, 0, 0, 0, 0, 0,  0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0,  module_potential_81_SOA_GPU, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};
#endif


__global__
void
module_potential_totalPotential_SOA_GPU(type_t *potential_GPU, const struct Potential_SOA *lens, const struct grid_param *frame, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
	struct point image_point;
	type_t potential, potential_temp;
	//
	int col = blockIdx.x*blockDim.x + threadIdx.x;
	int row = blockIdx.y*blockDim.y + threadIdx.y;
	//if(row == 0 && col == 0) printf("Start GPU \n");

	//
	if ((row + jstart < nbgridcells_y) && (col + istart < nbgridcells_x))
	{
		int index = row*nbgridcells_x + col;
		// Create temp pot variable to minimise writing to global memory potential
		potential = 0;
		potential_temp = 0;
		//
		image_point.x = frame->xmin + (col + istart)*dx;
		image_point.y = frame->ymin + (row + jstart)*dy;
		//
		int shalos = 0;
		//if(row == 0 && col == 0) printf("Start 2 GPU \n");
		//if(row == 0 && col == 0) std::cout << std::endl;;


		while (shalos < nhalos)
		{
			int lens_type = lens->type[shalos];
			int count     = 1;
			while (lens->type[shalos + count] == lens_type and shalos + count < nhalos) count++;
			//
			//clumpgrad = (*halo_func2_GPU[lens_type])(&image_point, lens, shalos, count);
			//clumpgrad = module_potentialDerivatives_totalGradient2_81_SOA_GPU(&image_point, lens, shalos, count);

			if(lens_type == 81) potential_temp = module_potential_81_SOA_GPU(&image_point, lens, shalos, count);
			//else if(lens_type == 14) clumpgrad = module_potentialDerivatives_totalGradient2_14_SOA_GPU(&image_point, lens, shalos, count);
			else if(row == 0 && col == 0) printf("No kernel selected \n");

			//
			potential += potential_temp;
			shalos += count;
		}
		//if(row == 0 && col == 0) printf(" %f %f %f %f \n",grad.a,grad.b,grad.c,grad.d);
		// Write to global memory
		potential_GPU[index] = potential;
		//if(row == 0 && col == 0) printf("point = %lf \n", grid_grad2_a[index] );

	}
}

