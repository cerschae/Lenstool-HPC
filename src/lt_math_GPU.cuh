/**
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch), Gilles Fourestey (gilles.fourestey@epfl.ch)
* @date   July 2017
* @version 0,1
*
*/
#ifndef LT_MATH_GPU_CUH_
#define LT_MATH_GPU_CUH_

//#include "cudafunctions.cuh"
#include <fstream>
#include <structure_hpc.hpp>
//#include <cuda.h>



/* COMPLEX FUNCTIONS ---------------------------------------------------*/

/*
* make complex number
* Global variables used :
* - none
*/

__device__ complex cpx_GPU(double re, double im)
{
    complex a;
    a.re = re;
    a.im = im;
    return(a);
}

/*
* square complex --------------------------------------------------
*/

__device__ complex sqcpx_GPU(complex c)
{
    complex res;

    res.re = c.re * c.re - c.im * c.im;
    res.im = 2.*c.re * c.im;

    return(res);
}

/*
* complex addition  --------------------------------------------------
* Global variables used :
* - none
*/
__device__ complex acpx_GPU(complex c1, complex c2)
{
    complex res;

    res.re = c1.re + c2.re;
    res.im = c1.im + c2.im;

    return(res);
}

/*
* complex substraction --------------------------------------------------
* Global variables used :
* - none
*/
__device__ complex scpx_GPU(complex c1, complex c2)
{
    complex res;

    res.re = c1.re - c2.re;
    res.im = c1.im - c2.im;

    return(res);
}

/*
* double substraction to complex ----------------------------------------
*/
__device__ complex scpxflt_GPU(complex c1, double f2)
{
    complex res;

    res.re = c1.re - f2;
    res.im = c1.im;

    return(res);
}

/*
* double addition to complex ----------------------------------------------
*/
__device__ complex acpxflt_GPU(complex c1, double f2)
{
    complex res;

    res.re = c1.re + f2;
    res.im = c1.im;

    return(res);
}

/*
* complex product----------------------------------------------
* Global variables used :
* - none
*/
__device__ complex pcpx_GPU(complex c1, complex c2)
{
    complex res;

    res.re = c1.re * c2.re - c1.im * c2.im;
    res.im = c1.im * c2.re + c2.im * c1.re;

    return(res);
}

/*
* complex times double ----------------------------------------------
* Global variables used :
* - none
*/
__device__ complex pcpxflt_GPU(complex c, double f)
{
    complex res;

    res.re = c.re * f;
    res.im = c.im * f;

    return(res);
}

/*
*  complex divided by double ------------------------------------------
* Global variables used :
* - none
*/
__device__ complex dcpxflt_GPU(complex c, double f)
{
    complex res;

    if (f != 0)
    {
        res.re = c.re / f;
        res.im = c.im / f;
        return(res);
    }
    else
    {
        //fprintf(stderr, "dcpxflt: Division by Zero!\n");
        return(c);
    }
}
/* norme complexe au carre------------------------------------------
 * Global variables used :
 * - none
 */
__device__ double ncpx_GPU(complex c)
{
    double res;

    res = c.re * c.re + c.im * c.im;

    return(res);
}

/* complex inverse ------------------------------------------
 * Global variables used :
 * - none
 */
__device__ complex icpx_GPU(complex c)
{
    complex res;
    double norm;

    norm = ncpx_GPU(c);

    if (norm != 0.)
    {
        res.re = c.re / norm;
        res.im = -c.im / norm;
    }
    else
    {
        //fprintf(stderr, "ERROR:(icpx) division by zero\n");
        res.re = 0.;
        res.im = 0.;
    };

    return(res);
}

/* division complexe ------------------------------------------
 * Return a complex number with
 * re = Real(c1*Conj(c2))/Norm(c2)
 * im = Img(c1*Conj(c2))/Norm(c2)
 *
 * If Norma(c2)=0 then return a complex = 0
 *
 * Global variables used :
 * - none
*/
__device__ complex dcpx_GPU(complex c1, complex c2)
{
    complex res;
    double norm;

    norm = ncpx_GPU(c2);

    if (norm != 0.)
    {
        res.re = (c1.re * c2.re + c1.im * c2.im) / norm;
        res.im = (c1.im * c2.re - c1.re * c2.im) / norm;
    }
    else
    {
        //fprintf(stderr, "ERROR:(dcpx) division by zero\n");
        res.re = 0.;
        res.im = 0.;
    };

    return(res);
}

/*
* square root complex --------------------------------------------------
* Global variables used :
* - none
*/

__device__ complex sqrtcpx_GPU(complex c)
{
    complex res;
    double   arg, nc;

    nc = sqrt(ncpx_GPU(c));
    arg = atan2(c.im, c.re);

    res.re = sqrt(nc) * cos(arg / 2.);
    res.im = sqrt(nc) * sin(arg / 2.);

    return(res);
}

/* Global variables used :
 * - none
 */
__device__ double  sgn_darg_GPU(complex z1, complex z2)
{
    double   arg1, arg2;

    arg1 = atan2(z1.im, z1.re);
    arg2 = atan2(z2.im, z2.re);

    /* fprintf (stderr,"%.3lf %.3lf %.3lf %.3lf %.2lf %.2lf %.2lf\n",
            z1.re,z1.im,z2.re,z2.im,arg1,arg2,fabs(arg1-arg2));
    */
    if (fabs(arg1 - arg2) < M_PI / 2.)
        return(1.);
    else
        return(-1.);
}

/* exponentiel complexe ------------------------------------------
*/
__device__ complex ecpx_GPU(complex c)
{
    complex res;
    double expo;

    expo = exp(c.re);
    res.re = expo * cos(c.im);
    res.im = expo * sin(c.im);

    return(res);
}

/* cosh complexe ------------------------------------------
*/
__device__ complex coshcpx_GPU(complex c)
{
    complex res;
    double  ex, emx;

    ex = exp(c.re);
    emx = exp(-c.re);
    res.re = 0.5 * (ex + emx) * cos(c.im);
    res.im = 0.5 * (ex - emx) * sin(c.im);
    return(res);
}

/* sinh complexe ------------------------------------------
*/
__device__ complex sinhcpx_GPU(complex c)
{
    complex res;
    double  ex, emx;

    ex = exp(c.re);
    emx = exp(-c.re);
    res.re = 0.5 * (ex - emx) * cos(c.im);
    res.im = 0.5 * (ex + emx) * sin(c.im);
    return(res);
}


/* cos complexe */
__device__ complex coscpx_GPU(complex c)
{
    complex res;
    double  ey, emy;

    ey = exp(c.im);
    emy = exp(-c.im);
    res.re = 0.5 * (ey + emy) * cos(c.re);
    res.im = -0.5 * (ey - emy) * sin(c.re);
    return(res);
}

/* sin complexe */
__device__ complex sincpx_GPU(complex c)
{
    complex res;
    double  ey, emy;

    ey = exp(c.im);
    emy = exp(-c.im);
    res.re = -0.5 * (ey - emy) * cos(c.re);
    res.re = 0.5 * (ey + emy) * sin(c.re);
    // TODO: Define res.im
    return(res);
}

/* tan complexe */
__device__ complex tancpx_GPU(complex c)
{
    return(dcpx_GPU(coscpx_GPU(c), sincpx_GPU(c)));
}

/* log complexe ------------------------------------------
*/
__device__ complex lncpx_GPU(complex c)
{
    complex res;

    res.re = log(sqrt(c.re * c.re + c.im * c.im));
    res.im = atan2(c.im, c.re);

    return(res);
}

/* END OF COMPLEX FUNCTIONS ----------------------------------------------*/


#endif
