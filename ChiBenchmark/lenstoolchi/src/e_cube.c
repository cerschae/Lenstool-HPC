#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"
#include<string.h>


void o_cube(double ***cube, double xmin, double xmax, double ymin, double ymax,
		double llmin, double llmax, int nx, int ny, int nz, 
	       	double scale, struct galaxie *source);

/****************************************************************/
/*      nom:        e_cube               */
/*      auteur:     Vincent Binet        */
/*      date:       26/06/2012           */
/*      place:      Lyon                 */
/****************************************************************
 * Create a FITS file containing the computed cube of the arcs
 * and arclets in the source plane
 */
void    e_cube(int np, int nz, char *iname, char *sname, struct galaxie *source)
{
	const extern struct g_frame   F;
     	const extern struct g_mode    M;
	const extern struct g_observ  O;
	//const extern    struct  g_large L;

        extern struct vfield vf;
	double xmin, ymin, xmax, ymax, llmin, llmax;
	double dx, dy;
	double f, scale;
	int nx, ny;
        int j,k,kk;
	
	double ***cube=NULL;
	double **ima_temp;
	const char termformat[] = ".fits";

	const extern struct g_cosmo C;
	double scalekpc;

	double xc,yc;

	xmin = F.xmin;
        xmax = F.xmax;
	ymin = F.ymin;
	ymax = F.ymax;
	llmin = F.lmin;
	llmax = F.lmax;

	dx = xmax - xmin;
	dy = ymax - ymin;	
	nx = np;    // default: assign np to nx and adjust ny
	scale = dx / (nx-1);
	ny = dy / scale + 1;  // warning: trunc(ny) behavior 
	//f = dx / scale;
	//nx = (int) f;  // BUG in C
	//f = dy / scale;
	//ny = (int) f;
	
	if (O.setbin)
	{
		if ( O.bin * nx <= NTMAX && O.bin * ny <= NTMAX )
		{
			nx = O.bin * nx;
			ny = O.bin * ny;
			scale = scale / ((double) O.bin);
		}
		else
		{
			fprintf(stderr, "ERROR: reaching maximal size of array. Increase NTMAX in dimension.h\n");	
			exit(-1);										
		}
	}

	dx = (nx - 1) * scale;
	dy = (ny - 1) * scale;
	xmax = xmin + dx;
	ymax = ymin + dy;

	NPRINTF(stderr, "\tcube (%d,%d,%d) s=%.3lf [%.3lf,%.3lf] [%.3lf,%.3lf] [%.3lf,%.3lf]\n",
			nx, ny, nz, scale, xmin, xmax, ymin, ymax,llmin,llmax);
	
	ima_temp = (double **) alloc_square_double(ny, nx);
	cube = (double ***)alloc_cubic_double(ny,nx,nz);

	scalekpc = 1.;//d0/C.h*distcosmo1(source[0].z);
	
	xc=scalekpc*(xmin+xmax)/2.;
	vf.C.x=xc;
	yc=scalekpc*(ymin+ymax)/2.;
	vf.C.y=yc;

	
	o_cube(cube ,xmin, xmax, ymin, ymax, llmin, llmax, nx, ny, nz, scale, source);
 
    if (O.setseeing || O.bruit)
    {
        for (kk=0  ; kk < nz ; kk ++)
        {
            for (j = 0; j < ny; j++)
                for (k = 0; k < nx; k++)
                     ima_temp[j][k]=cube[j][k][kk];

            if (O.setseeing)
                 d_seeing(ima_temp, nx, ny, scale);

            if (O.bruit)
                 d_bruiter(ima_temp, nx, ny);

            for (j = 0; j < ny; j++)
                for (k = 0; k < nx; k++)
                     cube[j][k][kk]=ima_temp[j][k];
        }
    }

        wrf_cube_fits(iname,cube,nx,ny,nz,xmin,xmax,ymin,ymax,llmin,llmax);
	free_square_double(ima_temp,ny);
	free_cubic_double(cube, ny, nx);
}


/* Function to simulate an image from a source
 * Called by o_chi() and e_pixel()
 * cube is filled in a field defined by champ section with xmin, ymin, llmin, xmax, ymax, llmax
 *
 */ 
void o_cube(double ***cube, double xmin, double xmax, double ymin,
	       	double ymax,double llmin, double llmax, int nx, int ny, int nz, 
	        double scale, struct galaxie *source)
{

/* variables pour le champ de vitesse  */
    
	const extern struct g_cosmo C;
	double  cosi, sini, cost, sint;
	double  scalekpc, doppler, normalisation; 
	double  sigma2, cosinus, sinus, dist, psxkpc, psykpc;
        double lstep;
    	int j, k, kk;
	double **zzz;

    	extern struct vfield vf;

        lstep=(llmax-llmin)/nz;

	cosi=cos(vf.i);
        sini=sin(vf.i);
        cost=cos(vf.theta);
        sint=sin(vf.theta);
    
        zzz = (double **) alloc_square_double(ny, nx);
	o_pixel(zzz, nx, ny, scale, xmin, xmax, ymin, ymax, source);

        scalekpc = d0/C.h*distcosmo1(source[0].z);

//#pragma omp parallel for schedule(static)
    for (j = 0; j < ny; j++)
    {
	struct point pi;
        pi.y = ymin + j * scale;

        for (k = 0; k < nx; k++)
        {
            struct point ps;
            pi.x = xmin + k * scale;

            e_dpl(&pi, source[0].dr, &ps);
	    psykpc=ps.y*scalekpc;
	    psxkpc=ps.x*scalekpc;
	    dist=sqrt(pow((psxkpc-vf.C.x)*cost+(psykpc-vf.C.y)*sint,2.)+
			    pow(((psykpc-vf.C.y)*cost-(psxkpc-vf.C.x)*sint)/cosi,2.));

	    cosinus=((psxkpc-vf.C.x)*cost+(psykpc-vf.C.y)*sint)/dist;
	    sinus=((psykpc-vf.C.y)*cost-(psxkpc-vf.C.x)*sint)/(cosi*dist);
	    if (vf.profile ==1)
	    {
		    doppler=1.+2.*vf.vt*atan(2.*dist/vf.rt)*cosinus*sini/(M_PI*vol);
	    }
	    else if (vf.profile ==2)
	    {
		    doppler=1.+2.*vf.vt*atan(2.*dist/vf.rt)*cosinus*sini/(M_PI*vol);
	    }
	    else if (vf.profile ==3)
	    {
		    doppler=1.+2.*vf.vt*atan(2.*dist/vf.rt)*cosinus*sini/(M_PI*vol);
	    }
	    else if (vf.profile==4)
	    {
		    doppler=1.+2.*vf.vt*atan(2.*dist/vf.rt)*cosinus*sini/(M_PI*vol); 
	    }
	    sigma2=4.;//pow(sini*l0/vol,2.)*(pow(cosinus*4.*vf.vt/(M_PI*vf.rt*
	    		//    (1.+pow(2.*dist/vf.rt,2.))),2.)+
	    	 //pow(2.*vf.vt*sinus*atan(2*dist/vf.rt)/(M_PI*dist),2.))/12.;
            sigma2=pow(vf.lcent*vf.sigma*(1+source[0].z)/vol,2.0);
	    normalisation=zzz[j][k]/sqrt(2.*M_PI*sigma2);

	    for ( kk = 0 ; kk < nz ; kk++ )
	    {
		    cube[j][k][kk]=normalisation*exp(-pow(llmin+(double)kk*lstep - vf.lcent*doppler,2.)/(2.*sigma2));
	    }
	}
    }

    free_square_double(zzz, ny);
}
