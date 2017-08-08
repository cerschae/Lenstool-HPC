#include <fstream>
#include "grid_gradient_GPU.cuh"
#include "gradient_GPU.cuh"

#define BLOCK_SIZE_X 8
#define BLOCK_SIZE_Y 16

//#define ROT

#define _SHARED_MEM

#ifdef _SHARED_MEM
#define SHARED __shared__
#warning "shared memory"
extern __shared__ double shared[];
#else
#define SHARED 
#endif

#define Nx 1
#define Ny 0

extern "C" {
double myseconds();
}
void
module_potentialDerivatives_totalGradient_SOA_CPU_GPU(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_cpu, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos);

void
module_potentialDerivatives_totalGradient_SOA_CPU_GPU_v2(double *grid_grad_x, double *grid_grad_y, const struct grid_param *frame, const struct Potential_SOA *lens_cpu, const struct Potential_SOA *lens_gpu, int nbgridcells, int nhalos);

void calculate_cossin_values(double *theta_cos, double *theta_sin, double *angles, int nhalos ){
	for(int i = 0 ; i < nhalos; i++)
	{
		theta_cos[i] = cos(angles[i]);
		theta_sin[i] = sin(angles[i]);
	}
}



__global__
void
module_potentialDerivatives_totalGradient_8_SOA_GPU_cur(double *grid_grad_x, double *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int shalos, int nhalos)
{
        //asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        //
        struct point grad, clumpgrad, image_point;
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
                double dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
                double dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
                //
                image_point.x = frame->xmin + col*dx;
                image_point.y = frame->ymin + row*dy;
		//
                for(int i = shalos; i < shalos + nhalos; i++)
                {
                        struct point true_coord, true_coord_rot; //, result;
                        //double       R, angular_deviation;
                        complex      zis;
                        //
                        // positionning at the potential center
                        // Change the origin of the coordinate system to the center of the clump
                        //
			//@@if ((row == Ny) && (col == Nx)) printf("image_x = %f, %f image_y = %f, %f\n",  image_point.x, frame->xmin, image_point.y,frame->ymin);
			double true_coord_x = image_point.x - __ldg(&lens->position_x[i]);
                        double true_coord_y = image_point.y - __ldg(&lens->position_y[i]);
			//if ((row == Ny) && (col == Nx)) printf("x = %f y = %f\n",  true_coord_x, true_coord_y);	
			//
                        double cosi = __ldg(&lens->anglecos[i]);
                        double sinu = __ldg(&lens->anglesin[i]);
			//
                        double x = true_coord_x*cosi + true_coord_y*sinu;
                        double y = true_coord_y*cosi - true_coord_x*sinu;
			//
			//if ((row == Ny) && (col == Nx)) printf("x = %f y = %f\n",  x, y);	
                        //
                        double eps = __ldg(&lens->ellipticity_potential[i]);
                        //
                        double sqe  = sqrt(eps);
                        //
                        double rem2 = x*x/((1. + eps)*(1. + eps)) + y*y/((1. - eps)*(1. - eps));
                        //
                        complex zci;
                        complex znum, zden, zres;
                        double norm;
                        //
			zci.re  = 0;
                        zci.im  = -0.5*(1. - eps*eps)/sqe;
			//@@if ((col == Nx) && (row == Ny)) printf("%d %d, zis: %f %f\n", row, col, zci.re, zci.im);
                        //
                        double rc  = __ldg(&lens->rcore[i]);
                        double cx1  = (1. - eps)/(1. + eps);
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
                        double b0  = __ldg(&lens->b0[i]);
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
module_potentialDerivatives_totalGradient_8_SOA_GPU_SM2(double *grid_grad_x, double *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int shalos, int nhalos)
{
	//
        //asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        //
        double grad_x, grad_y;
	double clumpgrad_x, clumpgrad_y;
	double image_point_x, image_point_y;
	//
	SHARED double cosi	[200];
	SHARED double sinu	[200];
	SHARED double rc	[200];
	SHARED double b0	[200];
	SHARED double epsi	[200];
	SHARED double position_x[200];
	SHARED double position_y[200];
	SHARED double rsqe	[200];
	SHARED double sonepeps	[200];
	SHARED double sonemeps	[200];
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
	double dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
	double dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
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
			double true_coord_x = image_point_x - position_x[i];
			double true_coord_y = image_point_y - position_y[i];
			//
			double x = true_coord_x*cosi[i] + true_coord_y*sinu[i];
			double y = true_coord_y*cosi[i] - true_coord_x*sinu[i];
			//
			double eps     = epsi[i]; 
			//double onemeps = 1 - eps; 
			//double onepeps = 1 + eps; 
			//
			//double eps     = epsi[i]; 
			double onemeps = sonemeps[i]; 
			double onepeps = sonepeps[i];
			//
			//double sqe  = sqrt(eps);
			double sqe  = rsqe[i];
			double rem2 = x*x/(onepeps*onepeps) + y*y/(onemeps*onemeps);
			//
			complex      zis;
			//
			double znum_re, znum_im;
			double zres_re, zres_im;
			double norm;
			double zden_re, zden_im;
			double  zis_re,  zis_im;
			//
			double zci_im  = -0.5*(1. - eps*eps)/sqe;
			//
			double cx1  = onemeps/onepeps; 
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
module_potentialDerivatives_totalGradient_8_SOA_GPU_SM3(double *grid_grad_x, double *grid_grad_y, const struct Potential_SOA
 *lens, const struct grid_param *frame, int nbgridcells, int shalos, int nhalos)
{
        //
        //asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
        // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
        //
        double grad_x, grad_y;
        double clumpgrad_x, clumpgrad_y;
        double image_point_x, image_point_y;
        //
        SHARED double cosi      [200];
        SHARED double sinu      [200];
        SHARED double rci       [200];
        SHARED double b0        [200];
        SHARED double epsi      [200];
        SHARED double position_x[200];
        SHARED double position_y[200];
        SHARED double rsqe      [200];
	//SHARED double sgrad_x   [(BLOCK_SIZE_X + 1)*BLOCK_SIZE_Y];
	//SHARED double sgrad_y   [(BLOCK_SIZE_X + 1)*BLOCK_SIZE_Y]; 

        //SHARED double sonepeps  [200];
        //SHARED double sonemeps  [200];
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
	double dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
	double dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
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
			double rc      = rci[i];
			double eps     = epsi[i];
			double onemeps = 1 - eps;
			double onepeps = 1 + eps;
			//
			double sqe = rsqe[i];
			double cx1 = onemeps/onepeps;
			//
			//
			//double zci_im = 1;
			double zci_im  = -0.5*(1. - eps*eps)/sqe;
			double inv_onepeps = 1./(onepeps*onepeps);
			double inv_onemeps = 1./(onemeps*onemeps);
#endif
			//
			{
				//KERNEL_8;
				double grad_x = grad_y = 0.;
				double true_coord_y = image_point_y - position_y[i];
				double true_coord_x = image_point_x - position_x[i];
				KERNEL_8_reg(0);
				grid_grad_x[iglob +  0] += grad_x;
				grid_grad_y[iglob +  0] += grad_y;
			}
			/*
			{
				//KERNEL_8;
				double grad_x = grad_y = 0.;
				double true_coord_y = image_point_y - position_y[i];
				double true_coord_x = image_point_x - position_x[i] + BLOCK_SIZE_X/2*dx;
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
module_potentialDerivatives_totalGradient_8_SOA_GPU_SM4(double *grid_grad_x, double *grid_grad_y, const struct Potential_SOA
		lens, const struct grid_param *frame, int nbgridcells, int shalos, int nhalos/*, double* dtimer*/)
{
	//
	//asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
	// 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
	//
	double grad_x, grad_y;
	double clumpgrad_x, clumpgrad_y;
	double image_point_x, image_point_y;
	//
	//
	double* cosi       = &shared[0*nhalos];
	double* sinu       = &shared[1*nhalos];
	double* rc         = &shared[2*nhalos];
	double* b0         = &shared[3*nhalos];
	double* epsi       = &shared[4*nhalos];
	double* position_x = &shared[5*nhalos];
	double* position_y = &shared[6*nhalos];
	double* rsqe       = &shared[7*nhalos];

	//SHARED double sonepeps  [200];
	//SHARED double sonemeps  [200];
	//
	grad_x         = 0;
	grad_y         = 0;
	//
	int grid_size  =  nbgridcells/gridDim.y;
	int row        =  blockIdx.x*blockDim.x + threadIdx.x;
	int col        = (blockIdx.y*blockDim.y + threadIdx.y)/**grid_size*/;
	//
	//
#if 1
	//SHARED double sgrad_x [
	if (/*(threadIdx.x == 0) &&*/ (blockIdx.x == 0))
	{
		if (threadIdx.x == 0) printf("blockDim.x = %d, blockIdx.x = %d grimdDim.x = %d threadIdx.x = %d\n", blockDim.x, blockIdx.x, gridDim.x, threadIdx.x);
		if (threadIdx.x == 0) printf("blockDim.y = %d, blockIdx.y = %d grimdDim.y = %d threadIdx.y = %d\n", blockDim.y, blockIdx.y, gridDim.y, threadIdx.y);
		if (threadIdx.x == 0) printf("row = %d, col = %d, grid_size = %d\n", row, col, grid_size);
	}
	__syncthreads();
#endif
	double* sgrad_x    = &shared[8*nhalos];
	double* sgrad_y    = &shared[8*nhalos + (grid_size + 1)*blockDim.x];
	//
	//
	//grid_grad_x[index] = 0.;
	//grid_grad_y[index] = 0.;
	//
	double dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
	double dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
	//
	//if (threadIdx.x == 0) printf("dx = %f, dy = %f\n", dx, dy);
	//
	image_point_x = frame->xmin + col*dx;
	image_point_y = frame->ymin + row*dy;
	return;
	//
	int i = 0;
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
                        double eps     = epsi[i];
                        double onemeps = 1 - eps;
                        double onepeps = 1 + eps;
                        //
                        double sqe = rsqe[i];
                        double cx1 = onemeps/onepeps;
                        //
                        //
                        //double x = true_coord_y*sinu[i];
                        //double y = true_coord_y*cosi[i];
                        double true_coord_y = image_point_y - position_y[i];
                        double true_coord_x = image_point_x - position_x[i] /*+ icol*dx*/;
                        //
                        complex zci;
                        zci.im  = -0.5*(1. - eps*eps)/sqe;
                        double inv_onepeps = 1./(onepeps*onepeps);
                        double inv_onemeps = 1./(onemeps*onemeps);
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
                                        double x = true_coord_x*cosi[i] + true_coord_y*sinu[i];
                                        double y = true_coord_y*cosi[i] - true_coord_x*sinu[i];
                                        //
                                        //if ((row == 1) && (col == 0)) printf("i = %d, eps = %f\n", i, eps);
                                        //
                                        //double eps     = epsi[i];
                                        //double onemeps = sonemeps[i];
                                        //double onepeps = sonepeps[i];
                                        //
                                        //double sqe  = sqrt(eps);
                                        //double rem2 = x*x/(onepeps*onepeps) + y*y/(onemeps*onemeps);
                                        double rem2 = x*x*inv_onepeps + y*y*inv_onemeps;
                                        //
                                        //
                                        //double znum_re, znum_im;
                                        //double zres_re, zres_im;
                                        //double zden_re, zden_im;
                                        //double  zis_re,  zis_im;
                                        double norm;
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
module_potentialDerivatives_totalGradient_8_SOA_GPU_v2(double *grid_grad_x, double *grid_grad_y, const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int i, int nhalos)
{
	//asm volatile("# module_potentialDerivatives_totalGradient_SOA begins");
	// 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
	//
	struct point grad, clumpgrad, image_point;
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
		double dx = (frame->xmax - frame->xmin)/(nbgridcells-1);
		double dy = (frame->ymax - frame->ymin)/(nbgridcells-1);
		//
#if 0
		/*SHARED*/ double img_pt[2];
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
		struct point true_coord, true_coord_rot; //, result;
		//double       R, angular_deviation;
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
		double cosi = __ldg(&lens->anglecos[i]);
		double sinu = __ldg(&lens->anglesin[i]);
		// positionning at the potential center
		// Change the origin of the coordinate system to the center of the clump
		double x = true_coord.x*cosi + true_coord.y*sinu;
		double y = true_coord.y*cosi - true_coord.x*sinu;
		//
		double eps = __ldg(&lens->ellipticity_potential[i]);
		//
		double sqe  = sqrt(eps);
		//
		double rem2 = x*x/((1. + eps)*(1. + eps)) + y*y/((1. - eps)*(1. - eps));
		//
		complex zci;
		complex znum, zden, zres;
		double norm;
		//
		zci.im  = -0.5*(1. - eps*eps)/sqe;
		//
		double rc  = __ldg(&lens->rcore[i]);
		double cx1  = (1. - eps)/(1. + eps);
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
		double b0  = __ldg(&lens->b0[i]);
		grad.x += b0*(zres.re*cosi - zres.im*sinu);
		grad.y += b0*(zres.im*cosi + zres.re*sinu);
		//
		grid_grad_x[index] += grad.x;
		grid_grad_y[index] += grad.y;
	}
}
