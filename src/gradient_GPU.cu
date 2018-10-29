/**
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)
* @date   July 2017
* @version 0,1
*
*/
#include <fstream>
#include "grid_gradient_GPU.cuh"
#include "gradient.hpp"
#include "gradient_GPU.cuh"
#include "gradient.hpp"

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
/*
void
module_potentialDerivatives_totalGradient_SOA_CPU_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_cpu, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos);

void
module_potentialDerivatives_totalGradient_SOA_CPU_GPU_v2(type_t *grid_grad_x, type_t *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_cpu, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos);

void calculate_cossin_values(type_t *theta_cos, type_t *theta_sin, type_t *angles, int nhalos ){
	for(int i = 0 ; i < nhalos; i++)
	{
		theta_cos[i] = cos(angles[i]);
		theta_sin[i] = sin(angles[i]);
	}
}
*/



#if 0

__global__
void
module_potentialDerivatives_totalGradient_8_SOA_GPU_cur(type_t *grid_grad_x, type_t *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int shalos, int nhalos)
{
        //asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        //
        struct point grad, image_point;
	//
        grad.x = 0;
        grad.y = 0;
        //
        int col = blockIdx.x*blockDim.x + threadIdx.x;
        int row = blockIdx.y*blockDim.y + threadIdx.y;
        //
        if ((row < nbgridcells) && (col < nbgridcells))
        {
                //
                int index = row*nbgridcells + col;
                //
                //grid_grad_x[index] = 0.;
                //grid_grad_y[index] = 0.;
                //
                type_t dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
                type_t dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
                //
                image_point.x = frame->xmin + col*dx;
                image_point.y = frame->ymin + row*dy;
		//
                for(int i = shalos; i < shalos + nhalos; i++)
                {
                        //type_t       R, angular_deviation;
                        complex      zis;
                        //
                        // positionning at the potential center
                        // Change the origin of the coordinate system to the center of the clump
                        //
			//@@if ((row == Ny) && (col == Nx)) printf("image_x = %f, %f image_y = %f, %f\n",  image_point.x, frame->xmin, image_point.y,frame->ymin);
			type_t true_coord_x = image_point.x - __ldg(&lens->position_x[i]);
                        type_t true_coord_y = image_point.y - __ldg(&lens->position_y[i]);
			//if ((row == Ny) && (col == Nx)) printf("x = %f y = %f\n",  true_coord_x, true_coord_y);	
			//
                        type_t cosi = __ldg(&lens->anglecos[i]);
                        type_t sinu = __ldg(&lens->anglesin[i]);
			//
                        type_t x = true_coord_x*cosi + true_coord_y*sinu;
                        type_t y = true_coord_y*cosi - true_coord_x*sinu;
			//
			//if ((row == Ny) && (col == Nx)) printf("x = %f y = %f\n",  x, y);	
                        //
                        type_t eps = __ldg(&lens->ellipticity_potential[i]);
                        //
                        type_t sqe  = sqrt(eps);
                        //
                        type_t rem2 = x*x/((1. + eps)*(1. + eps)) + y*y/((1. - eps)*(1. - eps));
                        //
                        complex zci;
                        complex znum, zden, zres;
                        type_t norm;
                        //
			zci.re  = 0;
                        zci.im  = -0.5*(1. - eps*eps)/sqe;
			//@@if ((col == Nx) && (row == Ny)) printf("%d %d, zis: %f %f\n", row, col, zci.re, zci.im);
                        //
                        type_t rc  = __ldg(&lens->rcore[i]);
                        type_t cx1  = (1. - eps)/(1. + eps);
                        znum.re = cx1*x;
                        znum.im = 2.*sqe*sqrt(rc*rc + rem2) - y/cx1;
                        //
                        zden.re = x;
                        zden.im = 2.*rc*sqe - y;
                        norm    = (zden.re*zden.re + zden.im*zden.im);     // zis = znum/zden
			//@@if ((col == Nx) && (row == Ny)) printf("norm = %f\n", norm);
                        //
                        zis.re  = (znum.re*zden.re + znum.im*zden.im)/norm;
                        zis.im  = (znum.im*zden.re - znum.re*zden.im)/norm;
			//
			//@@if ((col == Nx) && (row == Ny)) printf("%d %d, zis: %f %f\n", row, col, zis.re, zis.im);
                        //
                        norm    = zis.re;
                        //
                        zis.re  = log(sqrt(norm*norm + zis.im*zis.im));  // ln(zis) = ln(|zis|)+i.Arg(zis)
                        zis.im  = atan2(zis.im, norm);
			//
			//@@if ((col == Nx) && (row == Ny)) printf("%d %d, zis: %f %f\n", row, col, zis.re, zis.im);
                        //
                        zres.re = zci.re*zis.re - zci.im*zis.im;   // Re( zci*ln(zis) )
                        zres.im = zci.im*zis.re + zis.im*zci.re;   // Im( zci*ln(zis) )
			//
			//@@if ((col == Nx) && (row == Ny)) printf("%d %d, zres: %f %f\n", row, col, zres.re, zres.im);
                        //
                        type_t b0  = __ldg(&lens->b0[i]);
                        grad.x += b0*(zres.re*cosi - zres.im*sinu);
                        grad.y += b0*(zres.im*cosi + zres.re*sinu);
			//@@if ((col == Nx) && (row == Ny)) printf("grad: %f %f\n", grad.x, grad.y);
                }
                //IACA_END;
                //
                grid_grad_x[index] = grad.x;
                grid_grad_y[index] = grad.y;
		//if ((row == 0) && (col == 9)) 
		//printf("%f %f: %f %f\n",  image_point.x, image_point.y, grid_grad_x[index], grid_grad_y[index]);
        }
}

__global__
void
module_potentialDerivatives_totalGradient_8_SOA_GPU_SM2(type_t *grid_grad_x, type_t *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int shalos, int nhalos)
{
	//
        //asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        //
        type_t grad_x, grad_y;
	type_t clumpgrad_x, clumpgrad_y;
	type_t image_point_x, image_point_y;
	//
	SHARED type_t cosi	[200];
	SHARED type_t sinu	[200];
	SHARED type_t rc	[200];
	SHARED type_t b0	[200];
	SHARED type_t epsi	[200];
	SHARED type_t position_x[200];
	SHARED type_t position_y[200];
	SHARED type_t rsqe	[200];
	SHARED type_t sonepeps	[200];
	SHARED type_t sonemeps	[200];
        //
        grad_x = 0;
        grad_y = 0;
        //
        int col = blockIdx.x*blockDim.x + threadIdx.x;
        int row = blockIdx.y*blockDim.y + threadIdx.y;
	int ithread  = threadIdx.y*blockDim.x + threadIdx.x;
	//
	int index = row*nbgridcells + col;
	//
	//grid_grad_x[index] = 0.;
	//grid_grad_y[index] = 0.;
	//
	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
	type_t dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
	//
	image_point_x = frame->xmin + col*dx;
	image_point_y = frame->ymin + row*dy;
	//
	int i = ithread;
	if (i < nhalos)
	{
		cosi[i]       = __ldg(&lens->anglecos		  [shalos + i]);
		sinu[i]       = __ldg(&lens->anglesin		  [shalos + i]);
		position_x[i] = __ldg(&lens->position_x		  [shalos + i]);
		position_y[i] = __ldg(&lens->position_y		  [shalos + i]);
		rc[i]         = __ldg(&lens->rcore		  [shalos + i]);
		b0[i]         = __ldg(&lens->b0		          [shalos + i]);
		epsi[i]       = __ldg(&lens->ellipticity_potential[shalos + i]);
		//sonemeps[i]   = 1 - epsi[i];
		//sonepeps[i]   = 1 + epsi[i];
		rsqe[i]	      = sqrt(epsi[i]);
	}
	__syncthreads();
	//
	if ((row < nbgridcells) && (col < nbgridcells))
	{
		for(int i = 0; i < nhalos; i++)
		{
			//
			type_t true_coord_x = image_point_x - position_x[i];
			type_t true_coord_y = image_point_y - position_y[i];
			//
			type_t x = true_coord_x*cosi[i] + true_coord_y*sinu[i];
			type_t y = true_coord_y*cosi[i] - true_coord_x*sinu[i];
			//
			type_t eps     = epsi[i];
			//type_t onemeps = 1 - eps;
			//type_t onepeps = 1 + eps;
			//
			//type_t eps     = epsi[i];
			type_t onemeps = sonemeps[i];
			type_t onepeps = sonepeps[i];
			//
			//type_t sqe  = sqrt(eps);
			type_t sqe  = rsqe[i];
			type_t rem2 = x*x/(onepeps*onepeps) + y*y/(onemeps*onemeps);
			//
			complex      zis;
			//
			type_t znum_re, znum_im;
			type_t zres_re, zres_im;
			type_t norm;
			type_t zden_re, zden_im;
			type_t  zis_re,  zis_im;
			//
			type_t zci_im  = -0.5*(1. - eps*eps)/sqe;
			//
			type_t cx1  = onemeps/onepeps;
			//
			znum_re = cx1*x;
			znum_im = 2.*sqe*sqrt(rc[i]*rc[i] + rem2) - y/cx1;
			//
			zden_re = x;
			zden_im = 2.*rc[i]*sqe - y;
			//
			norm    = (x*x + zden_im*zden_im);     // zis = znum/zden
			zis.re  = (znum_re*x + znum_im*zden_im)/norm;
			zis.im  = (znum_im*x - znum_re*zden_im)/norm;
			//
			norm    = zis.re;
			//
			zis.re  = log(sqrt(norm*norm + zis.im*zis.im));  // ln(zis) = ln(|zis|)+i.Arg(zis)
			zis.im  = atan2(zis.im, norm);
			//
			zres_re = - zci_im*zis.im;   // Re( zci*ln(zis) )
			zres_im =   zci_im*zis.re;   // Im( zci*ln(zis) )
			//
			grad_x += b0[i]*(zres_re*cosi[i] - zres_im*sinu[i]);
			grad_y += b0[i]*(zres_im*cosi[i] + zres_re*sinu[i]);
		}
		//
		grid_grad_x[index] = grad_x;
		grid_grad_y[index] = grad_y;
		//__syncthreads();
	}
}
//
//
//
__global__
void
module_potentialDerivatives_totalGradient_8_SOA_GPU_SM3(type_t *grid_grad_x, type_t *grid_grad_y, const struct Potential_SOA
 *lens, const struct grid_param *frame, int nbgridcells, int shalos, int nhalos)
{
        //
        //asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        //
        type_t grad_x, grad_y;
        type_t clumpgrad_x, clumpgrad_y;
        type_t image_point_x, image_point_y;
        //
        SHARED type_t cosi      [200];
        SHARED type_t sinu      [200];
        SHARED type_t rci       [200];
        SHARED type_t b0        [200];
        SHARED type_t epsi      [200];
        SHARED type_t position_x[200];
        SHARED type_t position_y[200];
        SHARED type_t rsqe      [200];
	//SHARED type_t sgrad_x   [(BLOCK_SIZE_X + 1)*BLOCK_SIZE_Y];
	//SHARED type_t sgrad_y   [(BLOCK_SIZE_X + 1)*BLOCK_SIZE_Y];

        //SHARED type_t sonepeps  [200];
        //SHARED type_t sonemeps  [200];
        //
        grad_x         = 0;
        grad_y 	       = 0;
        //
        int row        = blockIdx.y*blockDim.y + threadIdx.y;
        int col        = blockIdx.x*blockDim.x + threadIdx.x;
	//
	//int loc_row    = threadIdx.x;
	//int loc_col    = threadIdx.y*blockDim.x + threadIdx.x;
	//
        //int grid_size  = nbgridcells/blockDim.y; 
	//
	//if (threadIdx.x == 0) printf("%d %d %d: row = %d, col = %d, grid_size = %d\n", blockIdx.y, gridDim.y, threadIdx.y, row, col, grid_size);
        //
	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
	type_t dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
	//if (threadIdx.x == 0) printf("dx = %f, dy = %f\n", dx, dy);
	//
	image_point_x = frame->xmin + col*dx;
	image_point_y = frame->ymin + row*dy;
	//
	//int iloc  = threadIdx.x*blockDim.y + threadIdx.y;
	int iglob = row*nbgridcells + col;
	int numThreads = blockDim.x*blockDim.y;
	//
	for (int i = 0; i < (nhalos + numThreads - 1)/numThreads; ++i)
	{ 
		int iloc  = threadIdx.y*blockDim.x + threadIdx.x + i*numThreads;
		if (iloc < nhalos)
		{
			cosi[iloc]       = __ldg(&lens->anglecos             [shalos + iloc]);
			sinu[iloc]       = __ldg(&lens->anglesin             [shalos + iloc]);
			position_x[iloc] = __ldg(&lens->position_x           [shalos + iloc]);
			position_y[iloc] = __ldg(&lens->position_y           [shalos + iloc]);
			rci[iloc]        = __ldg(&lens->rcore                [shalos + iloc]);
			b0[iloc]         = __ldg(&lens->b0                   [shalos + iloc]);
			epsi[iloc]       = __ldg(&lens->ellipticity_potential[shalos + iloc]);
			rsqe[iloc]       = sqrt(epsi[iloc]);
		}
	}
	__syncthreads();
	//
	if ((row < nbgridcells) && (col < nbgridcells))
	{
		//
		for(int i = 0; i < nhalos; i++)
		{
			//int index  = iloc; 
#if 1
			type_t rc      = rci[i];
			type_t eps     = epsi[i];
			type_t onemeps = 1 - eps;
			type_t onepeps = 1 + eps;
			//
			type_t sqe = rsqe[i];
			type_t cx1 = onemeps/onepeps;
			//
			//
			//type_t zci_im = 1;
			type_t zci_im  = -0.5*(1. - eps*eps)/sqe;
			type_t inv_onepeps = 1./(onepeps*onepeps);
			type_t inv_onemeps = 1./(onemeps*onemeps);
#endif
			//
			{
				//KERNEL_8;
				type_t grad_x = grad_y = 0.;
				type_t true_coord_y = image_point_y - position_y[i];
				type_t true_coord_x = image_point_x - position_x[i];
				KERNEL_8_reg(0);
				grid_grad_x[iglob +  0] += grad_x;
				grid_grad_y[iglob +  0] += grad_y;
			}
			/*
			{
				//KERNEL_8;
				type_t grad_x = grad_y = 0.;
				type_t true_coord_y = image_point_y - position_y[i];
				type_t true_coord_x = image_point_x - position_x[i] + BLOCK_SIZE_X/2*dx;
				KERNEL_8_reg(0);
				grid_grad_x[iglob + BLOCK_SIZE_X/2] += grad_x;
				grid_grad_y[iglob + BLOCK_SIZE_X/2] += grad_y;
			}
			*/
		//
		}
	}
}
//
//
//
__global__
void
module_potentialDerivatives_totalGradient_8_SOA_GPU_SM4(type_t *grid_grad_x, type_t *grid_grad_y, const struct Potential_SOA
		lens, const struct grid_param *frame, int nbgridcells, int shalos, int nhalos/*, type_t* dtimer*/)
{
	//
	//asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
	// 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
	//
	type_t grad_x, grad_y;
	type_t clumpgrad_x, clumpgrad_y;
	type_t image_point_x, image_point_y;
	//
	//
	type_t* cosi       = &shared[0*nhalos];
	type_t* sinu       = &shared[1*nhalos];
	type_t* rc         = &shared[2*nhalos];
	type_t* b0         = &shared[3*nhalos];
	type_t* epsi       = &shared[4*nhalos];
	type_t* position_x = &shared[5*nhalos];
	type_t* position_y = &shared[6*nhalos];
	type_t* rsqe       = &shared[7*nhalos];

	//SHARED type_t sonepeps  [200];
	//SHARED type_t sonemeps  [200];
	//
	grad_x         = 0;
	grad_y         = 0;
	//
	int grid_size  =  nbgridcells/gridDim.y;
	int row        =  blockIdx.x*blockDim.x + threadIdx.x;
	int col        = (blockIdx.y*blockDim.y + threadIdx.y)/**grid_size*/;
	//
	//

#if 0
	//SHARED double sgrad_x [
	if (/*(threadIdx.x == 0) &&*/ (blockIdx.x == 0))
	{
		if (threadIdx.x == 0) printf("blockDim.x = %d, blockIdx.x = %d grimdDim.x = %d threadIdx.x = %d\n", blockDim.x, blockIdx.x, gridDim.x, threadIdx.x);
		if (threadIdx.x == 0) printf("blockDim.y = %d, blockIdx.y = %d grimdDim.y = %d threadIdx.y = %d\n", blockDim.y, blockIdx.y, gridDim.y, threadIdx.y);
		if (threadIdx.x == 0) printf("row = %d, col = %d, grid_size = %d\n", row, col, grid_size);
	}
	__syncthreads();
#endif
	//type_t* sgrad_x    = &shared[8*nhalos];
	type_t* sgrad_y    = &shared[8*nhalos + (grid_size + 1)*blockDim.x];
	//
	//
	//grid_grad_x[index] = 0.;
	//grid_grad_y[index] = 0.;
	//
	type_t dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
	type_t dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
	//
	//if (threadIdx.x == 0) printf("dx = %f, dy = %f\n", dx, dy);
	//
	image_point_x = frame->xmin + col*dx;
	image_point_y = frame->ymin + row*dy;
	return;
	//
	//int i = 0;
#if 0
	for (; i < nhalos; i = i + blockDim.x)	
	{
		int pos = threadIdx.x + i;
		/*if ((threadIdx.x == 0) && (blockIdx.x == 0))*/ printf("pos = %d\n"); 
		__syncthreads();
		//
		cosi[pos]       = __ldg(&lens->anglecos             [shalos + pos]);
		sinu[pos]       = __ldg(&lens->anglesin             [shalos + pos]);
		position_x[pos] = __ldg(&lens->position_x           [shalos + pos]);
		position_y[pos] = __ldg(&lens->position_y           [shalos + pos]);
		rc[pos]         = __ldg(&lens->rcore                [shalos + pos]);
		b0[pos]         = __ldg(&lens->b0                   [shalos + pos]);
		epsi[pos]       = __ldg(&lens->ellipticity_potential[shalos + pos]);
		rsqe[pos]       = sqrt(epsi[i]);
		//
	}
#endif
#if 0
	if (threadIdx.x == 0)
		for (; i < nhalos; i += 1)
		{
			cosi[i]       = __ldg(&lens->anglecos             [shalos + i]);
			sinu[i]       = __ldg(&lens->anglesin             [shalos + i]);
			position_x[i] = __ldg(&lens->position_x           [shalos + i]);
			position_y[i] = __ldg(&lens->position_y           [shalos + i]);
			rc[i]         = __ldg(&lens->rcore                [shalos + i]);
			b0[i]         = __ldg(&lens->b0                   [shalos + i]);
			epsi[i]       = __ldg(&lens->ellipticity_potential[shalos + i]);
			rsqe[i]       = sqrt(epsi[i]);
		}
#endif
	__syncthreads();
	//if ((row == col == 0)) printf("shared mem done...\n");
	//
	if (row < nbgridcells)
	{
		//for(int icol = 0; icol < grid_size; ++icol){
                //      if (col + icol < nbgridcells){
                //grad_x = grad_y = 0.;
		int  index  = row*nbgridcells + col;
                //
                for(int i = 0; i < nhalos; i++)
                {
                        int sindex  = threadIdx.x*grid_size; 
#if 0
                        type_t eps     = epsi[i];
                        type_t onemeps = 1 - eps;
                        type_t onepeps = 1 + eps;
                        //
                        type_t sqe = rsqe[i];
                        type_t cx1 = onemeps/onepeps;
                        //
                        //
                        //type_t x = true_coord_y*sinu[i];
                        //type_t y = true_coord_y*cosi[i];
                        type_t true_coord_y = image_point_y - position_y[i];
                        type_t true_coord_x = image_point_x - position_x[i] /*+ icol*dx*/;
                        //
                        complex zci;
                        zci.im  = -0.5*(1. - eps*eps)/sqe;
                        type_t inv_onepeps = 1./(onepeps*onepeps);
                        type_t inv_onemeps = 1./(onemeps*onemeps);
#endif
                        //
                        for(int icol = 0; icol < grid_size; ++icol)
                        {
                                if (col + icol < nbgridcells)
                                {
#if 0
                                        if ((row == 1) && (col == 1)) printf("%d %d: %f %f\n", row, col, true_coord_x, true_coord_y);

                                        true_coord_x = image_point_x - position_x[i] + icol*dx;
                                        //
                                        //x += true_coord_x*cosi[i];
                                        //y -= true_coord_x*sinu[i];
                                        type_t x = true_coord_x*cosi[i] + true_coord_y*sinu[i];
                                        type_t y = true_coord_y*cosi[i] - true_coord_x*sinu[i];
                                        //
                                        //if ((row == 1) && (col == 0)) printf("i = %d, eps = %f\n", i, eps);
                                        //
                                        //double eps     = epsi[i];
                                        //double onemeps = sonemeps[i];
                                        //double onepeps = sonepeps[i];
                                        //
                                        //double sqe  = sqrt(eps);
                                        //double rem2 = x*x/(onepeps*onepeps) + y*y/(onemeps*onemeps);
                                        type_t rem2 = x*x*inv_onepeps + y*y*inv_onemeps;
                                        //
                                        //
                                        //double znum_re, znum_im;
                                        //double zres_re, zres_im;
                                        //double zden_re, zden_im;
                                        //double  zis_re,  zis_im;
                                        type_t norm;
                                        //
                                        complex      zis;
                                        complex      znum;
                                        complex      zden;
                                        complex      zres;
                                        //
                                        //double cx1  = onemeps/onepeps;
                                        //
                                        znum.re = cx1*x;
                                        znum.im = 2.*sqe*sqrt(rc[i]*rc[i] + rem2) - y/cx1;
                                        //
                                        zden.re = x;
                                        zden.im = 2.*rc[i]*sqe - y;
                                        //
                                        norm    = (x*x + zden.im*zden.im);     // zis = znum/zden
                                        zis.re  = (znum.re*x + znum.im*zden.im)/norm;
                                        zis.im  = (znum.im*x - znum.re*zden.im)/norm;
                                        //
                                        norm    = zis.re;
                                        //
                                        zis.re  = log(sqrt(norm*norm + zis.im*zis.im));  // ln(zis) = ln(|zis|)+i.Arg(zis)
                                        zis.im  = atan2(zis.im, norm);
                                        //
                                        zres.re = - zci.im*zis.im;   // Re( zci*ln(zis) )
                                        zres.im =   zci.im*zis.re;   // Im( zci*ln(zis) )
                                        //
                                        grid_grad_x[index] += b0[i]*(zres.re*cosi[i] - zres.im*sinu[i]);
                                        grid_grad_y[index] += b0[i]*(zres.im*cosi[i] + zres.re*sinu[i]);
#endif
                                        //sgrad_x[sindex] += (float)  sindex;
                                        sgrad_y[sindex] += (float) -sindex;
                                        //sindex++;
                                        //
                                        //grid_grad_x[index] += grad_x;
                                        //grid_grad_y[index] += grad_y;
                                }
                                //
                        }
                }
                        __syncthreads();
	return;
                //
#if 0
		int sindex = threadIdx.x*grid_size;
		for(int icol = 0; icol < grid_size; ++icol)
		{
			if (col + icol < nbgridcells)
			{
                		grid_grad_x[index + col] = sgrad_x[sindex];
                		grid_grad_y[index + col] = sgrad_y[sindex];
				sindex++;
			}
		}
		__syncthreads();
#endif
        }
}



__global__
void
module_potentialDerivatives_totalGradient_8_SOA_GPU_v2(type_t *grid_grad_x, type_t *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int i, int nhalos)
{
	//asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
	// 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
	//
	struct point grad, image_point;
	grad.x = 0;
	grad.y = 0;
	//
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	//
	if ((row < nbgridcells) && (col < nbgridcells))
	{
		//
		int index = col*nbgridcells + row;
		//
		//grid_grad_x[index] = 0.;
		//grid_grad_y[index] = 0.;
		//
		type_t dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
		type_t dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
		//
#if 0
		/*SHARED*/ type_t img_pt[2];
		if ((row == 0) && (col == 0))
		{
			img_pt[0] = frame->xmin + col*dx;
			img_pt[1] = frame->ymin + row*dy;
		}
		__syncthreads();
#else
		image_point.x = frame->xmin + col*dx;
		image_point.y = frame->ymin + row*dy;
#endif
		//
		//
		//for(int i = shalos; i < shalos + nhalos; i++)
		//{
		//IACA_START;
		//
		struct point true_coord; //, result;
		//type_t       R, angular_deviation;
		complex      zis;
		//
		//result.x = result.y = 0.;
		//
#if 0
		true_coord.x = img_pt[0] - __ldg(&lens->position_x[i]);
		true_coord.y = img_pt[1] - __ldg(&lens->position_y[i]);
#else
		true_coord.x = image_point.x - __ldg(&lens->position_x[i]);
		true_coord.y = image_point.y - __ldg(&lens->position_y[i]);
#endif
		type_t cosi = __ldg(&lens->anglecos[i]);
		type_t sinu = __ldg(&lens->anglesin[i]);
		// positionning at the potential center
		// Change the origin of the coordinate system to the center of the clump
		type_t x = true_coord.x*cosi + true_coord.y*sinu;
		type_t y = true_coord.y*cosi - true_coord.x*sinu;
		//
		type_t eps = __ldg(&lens->ellipticity_potential[i]);
		//
		type_t sqe  = sqrt(eps);
		//
		type_t rem2 = x*x/((1. + eps)*(1. + eps)) + y*y/((1. - eps)*(1. - eps));
		//
		complex zci;
		complex znum, zden, zres;
		type_t norm;
		//
		zci.im  = -0.5*(1. - eps*eps)/sqe;
		//
		type_t rc  = __ldg(&lens->rcore[i]);
		type_t cx1  = (1. - eps)/(1. + eps);
		znum.re = cx1*x;
		znum.im = 2.*sqe*sqrt(rc*rc + rem2) - y/cx1;
		//
		zden.re = x;
		zden.im = 2.*rc*sqe - y;
		norm    = (zden.re*zden.re + zden.im*zden.im);     // zis = znum/zden
		//
		zis.re  = (znum.re*zden.re + znum.im*zden.im)/norm;
		zis.im  = (znum.im*zden.re - znum.re*zden.im)/norm;
		//
		type_t b0  = __ldg(&lens->b0[i]);
		//grad.x += b0*(zres.re*cosi - zres.im*sinu);
		//grad.y += b0*(zres.im*cosi + zres.re*sinu);
		//
		grid_grad_x[index] += grad.x;
		grid_grad_y[index] += grad.y;
	}
}

#endif

__device__ point module_potentialDerivatives_totalGradient_5_SOA_GPU(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos){
    //asm volatile("# module_potentialDerivatives_totalGradient_SIS_SOA_v2 begins");
    //printf("# module_potentialDerivatives_totalGradient_SIS_SOA_v2 begins\n");
    //
    struct point grad, result;
    grad.x = 0;
    grad.y = 0;
    for(int i = shalos; i < shalos + nhalos; i++)
    {
            //
		struct point true_coord;
		//
		true_coord.x = pImage->x - lens->position_x[i];
		true_coord.y = pImage->y - lens->position_y[i];
		//
		//true_coord_rotation = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[i]);
		type_t cose = lens->anglecos[i];
		type_t sine = lens->anglesin[i];
		//
		type_t x = true_coord.x*cose + true_coord.y*sine;
		type_t y = true_coord.y*cose - true_coord.x*sine;
		//
		type_t ell_pot = lens->ellipticity_potential[i];
		//
		type_t b0_inv_R = lens->b0[i]/sqrt(x*x*(1 - ell_pot) + y*y*(1 + ell_pot));
		//
		result.x = (1 - ell_pot)*x*b0_inv_R;
		result.y = (1 + ell_pot)*y*b0_inv_R;
		//
		grad.x += result.x*cose - result.y*sine;
		grad.y += result.y*cose + result.x*sine;
    }
    return grad;

}
__device__ point module_potentialDerivatives_totalGradient_8_SOA_GPU(const struct point *image_point, const struct Potential_SOA *lens, int shalos, int nhalos){
	struct point grad;
	grad.x = 0;
	grad.y = 0;
	//
	for(int i = shalos; i < shalos + nhalos; i++)
	{
        complex      zis;
        // positioning at the potential center
        // Change the origin of the coordinate system to the center of the clump
        //@@if ((row == Ny) && (col == Nx)) printf("image_x = %f, %f image_y = %f, %f\n",  image_point.x, frame->xmin, image_point.y,frame->ymin);

        type_t true_coord_x = image_point->x - __ldg(&lens->position_x[i]);
        type_t true_coord_y = image_point->y - __ldg(&lens->position_y[i]);
        //if ((row == Ny) && (col == Nx)) printf("x = %f y = %f\n",  true_coord_x, true_coord_y);
        //
        type_t cosi = __ldg(&lens->anglecos[i]);
        type_t sinu = __ldg(&lens->anglesin[i]);
        //
        type_t x = true_coord_x*cosi + true_coord_y*sinu;
        type_t y = true_coord_y*cosi - true_coord_x*sinu;
        //
        //if ((row == Ny) && (col == Nx)) printf("x = %f y = %f\n",  x, y);
        //
        type_t eps = __ldg(&lens->ellipticity_potential[i]);
        //
        type_t sqe  = sqrt(eps);
        //
        type_t rem2 = x*x/((1. + eps)*(1. + eps)) + y*y/((1. - eps)*(1. - eps));
        //
        complex zci;
        complex znum, zden, zres;
        type_t norm;
        //
        zci.re  = 0;
        zci.im  = -0.5*(1. - eps*eps)/sqe;
        //@@if ((col == Nx) && (row == Ny)) printf("%d %d, zis: %f %f\n", row, col, zci.re, zci.im);
        //
        type_t rc  = __ldg(&lens->rcore[i]);
        type_t cx1  = (1. - eps)/(1. + eps);

        znum.re = cx1*x;
        znum.im = 2.*sqe*sqrt(rc*rc + rem2) - y/cx1;
        //
        zden.re = x;
        zden.im = 2.*rc*sqe - y;
        norm    = (zden.re*zden.re + zden.im*zden.im);     // zis = znum/zden
        //@@if ((col == Nx) && (row == Ny)) printf("norm = %f\n", norm);
        //
        zis.re  = (znum.re*zden.re + znum.im*zden.im)/norm;
        zis.im  = (znum.im*zden.re - znum.re*zden.im)/norm;
        //
        //@@if ((col == Nx) && (row == Ny)) printf("%d %d, zis: %f %f\n", row, col, zis.re, zis.im);
        //
        norm    = zis.re;
        //
        zis.re  = log(sqrt(norm*norm + zis.im*zis.im));  // ln(zis) = ln(|zis|)+i.Arg(zis)
        zis.im  = atan2(zis.im, norm);
        //
        //@@if ((col == Nx) && (row == Ny)) printf("%d %d, zis: %f %f\n", row, col, zis.re, zis.im);
        //
        zres.re = zci.re*zis.re - zci.im*zis.im;   // Re( zci*ln(zis) )
        zres.im = zci.im*zis.re + zis.im*zci.re;   // Im( zci*ln(zis) )
        //
        //@@if ((col == Nx) && (row == Ny)) printf("%d %d, zres: %f %f\n", row, col, zres.re, zres.im);
        //
        type_t b0  = __ldg(&lens->b0[i]);

        grad.x += b0*(zres.re*cosi - zres.im*sinu);
        grad.y += b0*(zres.im*cosi + zres.re*sinu);
        //@@if ((col == Nx) && (row == Ny)) printf("grad: %f %f\n", grad.x, grad.y);
	}
	//IACA_END;
	//
	return(grad);

}
__device__ point module_potentialDerivatives_totalGradient_81_SOA_GPU(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos){
    //asm volatile("# module_potentialDerivatives_totalGradient_81_SOA begins");
    //std::cout << "# module_potentialDerivatives_totalGradient_81_SOA begins" << std::endl;
    // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
    //
    struct point grad;
    grad.x = 0;
    grad.y = 0;
    for(int i = shalos; i < shalos + nhalos; i++)
    {
		//IACA_START;
		//
		struct point true_coord; //, result;
		//type_t       R, angular_deviation;
		complex      zis;
		//
		//result.x = result.y = 0.;
		//
		true_coord.x = pImage->x - lens->position_x[i];
		true_coord.y = pImage->y - lens->position_y[i];
		/*positionning at the potential center*/
		// Change the origin of the coordinate system to the center of the clump
		type_t cose = lens->anglecos[i];
		type_t sine = lens->anglesin[i];
		type_t x = true_coord.x*cose + true_coord.y*sine;
		type_t y = true_coord.y*cose - true_coord.x*sine;
		//
		type_t eps  = lens->ellipticity_potential[i];
		type_t rc   = lens->rcore[i];
		type_t rcut = lens->rcut[i];
		type_t b0   = lens->b0[i];
		type_t t05  = b0*rcut/(rcut - rc);
		//
		type_t sqe  = sqrt(eps);
		//
		type_t cx1  = (1. - eps)/(1. + eps);
		type_t cxro = (1. + eps)*(1. + eps);
		type_t cyro = (1. - eps)*(1. - eps);
		//
		type_t rem2 = x*x/cxro + y*y/cyro;
		//
		complex zci, znum, zden, zres_rc, zres_rcut;
		type_t norm;
		//
		zci.re  = 0;
		zci.im  = -0.5*(1. - eps*eps)/sqe;
		//
		// step 1
		{
		KERNEL(rc, zres_rc)
		}
		// step 2
		{
		KERNEL(rcut, zres_rcut)
		}
		zis.re  = t05*(zres_rc.re - zres_rcut.re);
		zis.im  = t05*(zres_rc.im - zres_rcut.im);
		// rotation
		grad.x += (zis.re*cose - zis.im*sine);
		grad.y += (zis.im*cose + zis.re*sine);
            //
    }
    //
    return(grad);

}

#if 1
typedef struct point (*halo_func_GPU_t) (const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);

__constant__ halo_func_GPU_t halo_func_GPU[100] =
{
	0, 0, 0, 0, 0, module_potentialDerivatives_totalGradient_5_SOA_GPU, 0, 0, module_potentialDerivatives_totalGradient_8_SOA_GPU,  0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	   0,  module_potentialDerivatives_totalGradient_81_SOA_GPU, 0, 0, 0, 0, 0, 0, 0, 0,
	   0, 0, 0, 0, 0, 0, 0, 0, 0, 0
	   };
#endif

__global__
void
module_potentialDerivatives_totalGradient_SOA_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
        struct point grad, clumpgrad, image_point;
        //
        int col = blockIdx.x*blockDim.x + threadIdx.x;
        int row = blockIdx.y*blockDim.y + threadIdx.y;

        //
        if ((row + 0*jstart < nbgridcells_y) && (col + 0*istart < nbgridcells_x))
        {
                //
                //double dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
                //double dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
                //
                int index = row*nbgridcells_x + col;
                // Create temp grad variable to minimise writing to global memory grid_grad
                grad.x = 0.;
                grad.y = 0.;
                //
                image_point.x = frame->xmin + (col + istart)*dx;
                image_point.y = frame->ymin + (row + jstart)*dy;
                //
                int shalos = 0;
                while (shalos < nhalos)
                {
                        int lens_type = lens->type[shalos];
                        int count     = 1;
                        while (lens->type[shalos + count] == lens_type) count++;
                        //
                        //if(row == 0 && col == 0) printf("type = %d, count %d , shalos %d \n", lens_type,count,shalos );
                        //
                        clumpgrad = (*halo_func_GPU[lens_type])(&image_point, lens, shalos, count);
                        //
                        grad.x += clumpgrad.x;
                        grad.y += clumpgrad.y;
                        shalos += count;
                }
                // Write to global memory
                grid_grad_x[index] = grad.x;
                grid_grad_y[index] = grad.y;
                //if ((row == 0) && (col == 9))
                //printf("%f %f: %f %f\n",  image_point.x, image_point.y, grid_grad_x[index], grid_grad_y[index]);
        }
}

__global__
void
module_potentialDerivatives_totalGradient_SOA_GPU(type_t *grid_grad_x, type_t *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int nhalos)
{
    struct point grad, clumpgrad, image_point;
    //
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;

    //
    if ((row < nbgridcells) && (col < nbgridcells))
    {
        //
        type_t dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
        type_t dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
		//
		int index = row*nbgridcells + col;
		// Create temp grad variable to minimise writing to global memory grid_grad
		grad.x = 0;
		grad.y = 0;
		//
		image_point.x = frame->xmin + col*dx;
		image_point.y = frame->ymin + row*dy;
		//
		int shalos = 0;
		while (shalos < nhalos)
		{
			int lens_type = lens->type[shalos];
			int count     = 1;
			while (lens->type[shalos + count] == lens_type) count++;
			//
			//if(row == 0 && col == 0) printf("type = %d, count %d , shalos %d \n", lens_type,count,shalos );
			//
			clumpgrad = (*halo_func_GPU[lens_type])(&image_point, lens, shalos, count);
			//
			grad.x += clumpgrad.x;
			grad.y += clumpgrad.y;
			shalos += count;
		}
		// Write to global memory
		grid_grad_x[index] = grad.x;
		grid_grad_y[index] = grad.y;

	//if ((row == 0) && (col == 9))
	//printf("%f %f: %f %f\n",  image_point.x, image_point.y, grid_grad_x[index], grid_grad_y[index]);
    }

}
