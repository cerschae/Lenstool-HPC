#include <fstream>
#include "grid_gradient2_GPU.cuh"
#include "gradient2_GPU.cuh"
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

__device__ matrix mdci05_gpu(type_t x, type_t y, type_t eps, type_t rc, type_t b0)
{
	matrix res;
    type_t   ci, sqe, cx1, cxro, cyro, wrem;
    type_t  didyre, didyim, didxre;// didxim;
    type_t  cx1inv, den1, num2, den2;

    sqe = sqrt(eps);
    cx1 = (1. - eps) / (1. + eps);
    cx1inv = 1. / cx1;
    cxro = (1. + eps) * (1. + eps);     /* rem^2=x^2/(1+e^2) + y^2/(1-e^2) Eq 2.3.6*/
    cyro = (1. - eps) * (1. - eps);
    ci = 0.5 * (1. - eps * eps) / sqe;
    wrem = sqrt(rc * rc + x * x / cxro + y * y / cyro); /*wrem^2=w^2+rem^2 with w core radius*/
    den1 = 2.*sqe * wrem - y * cx1inv;
    den1 = cx1 * cx1 * x * x + den1 * den1;
    num2 = 2.*rc * sqe - y;
    den2 = x * x + num2 * num2;

    didxre = ci * ( cx1 * (2.*sqe * x * x / cxro / wrem - 2.*sqe * wrem + y * cx1inv) / den1 + num2 / den2 );
    didyre = ci * ( (2 * sqe * x * y * cx1 / cyro / wrem - x) / den1 + x / den2 );

    didyim = ci * ( (2 * sqe * wrem * cx1inv - y * cx1inv * cx1inv - 4 * eps * y / cyro +
                     2 * sqe * y * y / cyro / wrem * cx1inv) / den1 - num2 / den2 );
    //if(blockIdx.x*blockDim.x + threadIdx.x == 0 && blockIdx.y*blockDim.y + threadIdx.y == 0) printf("didxre %f didyre %f didyim %f ",didxre,didyre,didyim);
    res.a = b0 * didxre;
    res.b = res.d = b0 * didyre; //(didyre+didxim)/2.;
    res.c = b0 * didyim;

    return(res);


}


__device__ void printmat_gpu(matrix A){
	printf("A: %lf, B: %lf, C:%lf, D:%lf \n",A.a,A.b,A.c,A.d);
}


__device__ matrix module_potentialDerivatives_totalGradient2_81_SOA_GPU(const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos){
    //asm volatile("# module_potentialDerivatives_totalGradient_81_SOA begins");

    //std::cout << "# module_potentialDerivatives_totalGradient_81_SOA begins" << std::endl;
    // 6 DP loads, i.e. 48 Bytes: position_x, position_y, ellipticity_angle, ellipticity_potential, rcore, b0
    //
    type_t t05, vala,valb,valc,vald;
     struct matrix grad2, clump, clumpcore, clumpcut;
     grad2.a =  0;
     grad2.b =  0;
     grad2.c =  0;
     grad2.d =  0;
#if 1

     for(int i = shalos; i < shalos + nhalos; i++)
     {

		struct point true_coord, true_coord_rotation;
		//True coord
		true_coord.x = pImage->x - lens->position_x[i];
		true_coord.y = pImage->y - lens->position_y[i];
		//

		//std::cerr << "start" << std::endl;
		//if(blockIdx.x*blockDim.x + threadIdx.x == 0 && blockIdx.y*blockDim.y + threadIdx.y == 0) printf("start\n");
		//Rotation
		type_t cose = lens->anglecos[i];
		type_t sine = lens->anglesin[i];
		//
		type_t x = true_coord.x*cose + true_coord.y*sine;
		type_t y = true_coord.y*cose - true_coord.x*sine;
		// 81 comput
		t05 = lens->rcut[i] / (lens->rcut[i] - lens->rcore[i]);

		//if(blockIdx.x*blockDim.x + threadIdx.x == 0 && blockIdx.y*blockDim.y + threadIdx.y == 0)printf("t05 %f rcut: %f, rcore: %f, b0:%f \n",t05, lens->rcut[i],lens->rcore[i], lens->b0[i]);
		clumpcore = mdci05_gpu(x, y, lens->ellipticity_potential[i], lens->rcore[i], lens->b0[i]);
		////////////////////////////
#if 0
		//matrix res;
		type_t rc = lens->rcore[i];
		type_t eps =lens->ellipticity_potential[i];

	    type_t   ci, sqe, cx1, cxro, cyro, wrem;
	    type_t  didyre, didyim, didxre;// didxim;
	    type_t  cx1inv, den1, num2, den2;

	    sqe = sqrt( lens->ellipticity_potential[i]);
	    cx1 = (1. - eps) / (1. + eps);
	    cx1inv = 1. / cx1;
	    cxro = (1. + eps) * (1. + eps);     /* rem^2=x^2/(1+e^2) + y^2/(1-e^2) Eq 2.3.6*/
	    cyro = (1. - eps) * (1. - eps);
	    ci = 0.5 * (1. - eps * eps) / sqe;

	    wrem = sqrt(rc * rc + x * x / cxro + y * y / cyro); /*wrem^2=w^2+rem^2 with w core radius*/
	    den1 = 2.*sqe * wrem - y * cx1inv;
	    den1 = cx1 * cx1 * x * x + den1 * den1;
	    num2 = 2.*rc * sqe - y;
	    den2 = x * x + num2 * num2;

	    didxre = ci * ( cx1 * (2.*sqe * x * x / cxro / wrem - 2.*sqe * wrem + y * cx1inv) / den1 + num2 / den2 );
	    if(blockIdx.x*blockDim.x + threadIdx.x == 0 && blockIdx.y*blockDim.y + threadIdx.y == 0) printf("cyro %f ",cyro);
	    if(cyro == 0) printf("I %d fucked up!",threadIdx.x);
	    //didyre = ci * ( (2 * sqe * x * y * cx1 / wrem - x) + x  );
	    didyre = ci * ( (2 * sqe * x * y * cx1 / wrem - x) / den1 + x / den2 );
	    didyre = didyre /cyro;
	    //didyre = ci * ( (2 * sqe * x * y * cx1 / cyro / wrem - x) / den1 + x / den2 );

	    didyim = ci * ( (2 * sqe * wrem * cx1inv - y * cx1inv * cx1inv - 4 * eps * y / cyro +
	                     2 * sqe * y * y / cyro / wrem * cx1inv) / den1 - num2 / den2 );

	    //if(blockIdx.x*blockDim.x + threadIdx.x == 0 && blockIdx.y*blockDim.y + threadIdx.y == 0) printf("didxre %f didyre %f didyim %f ",didxre,didyre,didyim);

	    clumpcore.a = lens->b0[i] * didxre;
	    clumpcore.b = clumpcore.d = lens->b0[i] * didyre; //(didyre+didxim)/2.;
	    clumpcore.c = lens->b0[i] * didyim;
#endif
	    ////////////////////////////////////////////////
		clumpcut = mdci05_gpu(x, y, lens->ellipticity_potential[i], lens->rcut[i], lens->b0[i]);
		/////////////////////////////////////
#if 0
		rc = lens->rcut[i];
	    sqe = sqrt( lens->ellipticity_potential[i]);
	    cx1 = (1. - eps) / (1. + eps);
	    cx1inv = 1. / cx1;
	    cxro = (1. + eps) * (1. + eps);     /* rem^2=x^2/(1+e^2) + y^2/(1-e^2) Eq 2.3.6*/
	    cyro = (1. - eps) * (1. - eps);
	    ci = 0.5 * (1. - eps * eps) / sqe;
	    wrem = sqrt(rc * rc + x * x / cxro + y * y / cyro); /*wrem^2=w^2+rem^2 with w core radius*/
	    den1 = 2.*sqe * wrem - y * cx1inv;
	    den1 = cx1 * cx1 * x * x + den1 * den1;
	    num2 = 2.*rc * sqe - y;
	    den2 = x * x + num2 * num2;

	    didxre = ci * ( cx1 * (2.*sqe * x * x / cxro / wrem - 2.*sqe * wrem + y * cx1inv) / den1 + num2 / den2 );
	    didyre = ci * ( (2 * sqe * x * y * cx1 / cyro / wrem - x) / den1 + x / den2 );

	    didyim = ci * ( (2 * sqe * wrem * cx1inv - y * cx1inv * cx1inv - 4 * eps * y / cyro +
	                     2 * sqe * y * y / cyro / wrem * cx1inv) / den1 - num2 / den2 );
	    //if(blockIdx.x*blockDim.x + threadIdx.x == 0 && blockIdx.y*blockDim.y + threadIdx.y == 0) printf("didxre %f didyre %f didyim %f ",didxre,didyre,didyim);
	    clumpcut.a = lens->b0[i] * didxre;
	    clumpcut.b = clumpcut.d = lens->b0[i] * didyre; //(didyre+didxim)/2.;
	    clumpcut.c = lens->b0[i] * didyim;
#endif
		///////////
		//if(blockIdx.x*blockDim.x + threadIdx.x == 0 && blockIdx.y*blockDim.y + threadIdx.y == 0) printmat_gpu(clumpcore);
		//if(blockIdx.x*blockDim.x + threadIdx.x == 0 && blockIdx.y*blockDim.y + threadIdx.y == 0) printmat_gpu(clumpcut);
		//
#if 1
		clumpcore.a = t05 * (clumpcore.a - clumpcut.a);
		clumpcore.b = t05 * (clumpcore.b - clumpcut.b);
		clumpcore.c = t05 * (clumpcore.c - clumpcut.c);
		clumpcore.d = t05 * (clumpcore.d - clumpcut.d);
		//if(blockIdx.x*blockDim.x + threadIdx.x == 0 && blockIdx.y*blockDim.y + threadIdx.y == 0) printmat_gpu(clumpcore);

		//rotation matrix  1
		clumpcut.a = clumpcore.a * cose + clumpcore.b * sine;
		clumpcut.b = clumpcore.a * -sine + clumpcore.b * cose;
		clumpcut.c = clumpcore.d * -sine + clumpcore.c * cose;
		clumpcut.d = clumpcore.d * cose + clumpcore.c * sine;
		//if(blockIdx.x*blockDim.x + threadIdx.x == 0 && blockIdx.y*blockDim.y + threadIdx.y == 0) printmat_gpu(clumpcut);
		//rotation matrix  2
		clump.a = cose * clumpcut.a + sine * clumpcut.d;
		clump.b = cose * clumpcut.b + sine * clumpcut.c;
		clump.c = -sine * clumpcut.b + cose * clumpcut.c;
		clump.d = -sine * clumpcut.a + cose * clumpcut.d;
#endif
		//if(blockIdx.x*blockDim.x + threadIdx.x == 0 && blockIdx.y*blockDim.y + threadIdx.y == 0) printmat_gpu(clump);
		//vala += clump.a;
		//valb += clump.b;
		//valc += clump.c;
		//vald += clump.d;
		grad2.a += clump.a;
		grad2.b += clump.b;
		grad2.c += clump.c;
		grad2.d += clump.d;
		//if(blockIdx.x*blockDim.x + threadIdx.x == 0 && blockIdx.y*blockDim.y + threadIdx.y == 0) printmat_gpu(grad2);

     }
     //
     //grad2.a =  vala;
     //grad2.b =  valb;
     //grad2.c =  valc;
     //grad2.d =  vald;
#endif
     return(grad2);

}
#if 1
typedef struct matrix (*halo_func2_GPU_t) (const struct point *pImage, const struct Potential_SOA *lens, int shalos, int nhalos);

__constant__ halo_func2_GPU_t halo_func2_GPU[100] =
{
	0, 0, 0, 0, 0, 0, 0, 0, 0,  0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	   0,  module_potentialDerivatives_totalGradient2_81_SOA_GPU, 0, 0, 0, 0, 0, 0, 0, 0,
	   0, 0, 0, 0, 0, 0, 0, 0, 0, 0
	   };
#endif
__global__
void
module_potentialDerivatives_totalGradient2_SOA_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct Potential_SOA *lens, const struct grid_param *frame, int nhalos, type_t dx, type_t dy, int nbgridcells_x, int nbgridcells_y, int istart, int jstart)
{
    struct point image_point;
	struct matrix grad, clumpgrad;
        //
        int col = blockIdx.x*blockDim.x + threadIdx.x;
        int row = blockIdx.y*blockDim.y + threadIdx.y;
        //if(row == 0 && col == 0) printf("Start GPU \n");

        //
        if ((row + jstart < nbgridcells_y) && (col + istart < nbgridcells_x))
        {
                int index = row*nbgridcells_x + col;
                // Create temp grad variable to minimise writing to global memory grid_grad
        		grad.a = 0;
        		grad.b = 0;
        		grad.c = 0;
        		grad.d = 0;
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
                        //if(row == 0 && col == 0) printf("point = %p \n", halo_func2_GPU[lens_type] );
                        //
                        //clumpgrad = (*halo_func2_GPU[lens_type])(&image_point, lens, shalos, count);
                        clumpgrad = module_potentialDerivatives_totalGradient2_81_SOA_GPU(&image_point, lens, shalos, count);

                        //
            			grad.a += clumpgrad.a;
            			grad.b += clumpgrad.b;
            			grad.c += clumpgrad.c;
            			grad.d += clumpgrad.d;
            			shalos += count;
            		}
            		// Write to global memory
            		grid_grad2_a[index] = grad.a;
            		grid_grad2_b[index] = grad.b;
            		grid_grad2_c[index] = grad.c;
            		grid_grad2_d[index] = grad.d;

        }
}

__global__
void
module_potentialDerivatives_totalGradient2_SOA_GPU(type_t *grid_grad2_a, type_t *grid_grad2_b, type_t *grid_grad2_c, type_t *grid_grad2_d, const struct Potential_SOA *lens, const struct grid_param *frame, int nbgridcells, int nhalos)
{
    struct point image_point;
	struct matrix grad, clumpgrad;
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
		grad.a = 0;
		grad.b = 0;
		grad.c = 0;
		grad.d = 0;
		//
		image_point.x = frame->xmin + col*dx;
		image_point.y = frame->ymin + row*dy;
		//
		int shalos = 0;
		while (shalos < nhalos)
		{
			int lens_type = lens->type[shalos];
			int count     = 1;
			while (lens->type[shalos + count] == lens_type and shalos + count < nhalos) count++;
			//
			clumpgrad = (*halo_func2_GPU[lens_type])(&image_point, lens, shalos, count);
			//
			//clumpgrad = module_potentialDerivatives_totalGradient2_81_SOA_GPU(&image_point, lens, shalos, count);
			grad.a += clumpgrad.a;
			grad.b += clumpgrad.b;
			grad.c += clumpgrad.c;
			grad.d += clumpgrad.d;
			shalos += count;
		}
		// Write to global memory
		grid_grad2_a[index] = grad.a;
		grid_grad2_b[index] = grad.b;
		grid_grad2_c[index] = grad.c;
		grid_grad2_d[index] = grad.d;
	//if ((row == 0) && (col == 9))
	//printf("%f %f: %f %f\n",  image_point.x, image_point.y, grid_grad_x[index], grid_grad_y[index]);
    }

}
