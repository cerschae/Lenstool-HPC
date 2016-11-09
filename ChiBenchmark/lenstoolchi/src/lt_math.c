#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include"dimension.h"
#include "structure.h"
#include "fonction.h"
#include "lt.h"

/* Return the value of the main mode of a pdf
 * Use the Freedman & Diaconis rule to find the best bin size */
double mode(int n, double *array)
{
    int i;
    double *copy, *histo;
    double binsize;
    int bin, nbin, biNMAX;
    double xmax;


    // make a copy of array
    copy = (double *)calloc(n, sizeof(double));
    for ( i = 0; i < n; i++ )
        copy[i] = array[i];

    arraySort(copy, n );
    binsize = copy[(int)(n*0.75)] - copy[(int)(n*0.25)];
    binsize *= 2 * pow(n, -0.33333);
    nbin = (int) round( (copy[n-1] - copy[0]) / binsize );

    histo = (double *)calloc(nbin, sizeof(double));

    xmax = copy[0] + binsize;
    bin = 0;
    for ( i = 0; i < n; i++ )
    {
        histo[bin]++;
        if ( copy[i] > xmax && bin < nbin - 1 )
        {
            histo[bin]--;
            bin++;
            histo[bin]++;
            xmax += binsize;
        }
    }

    // Find the index of the largest bin
    xmax = 0;
    biNMAX = 0;
    for ( i = 0; i < nbin; i++ )
        if ( histo[i] > xmax )
        {
            xmax = histo[i];
            biNMAX = i;
        }

    xmax =  copy[0] + binsize / 2. + biNMAX * binsize;

    free(histo);
    free(copy);


    return(xmax);
}

/* Return the mean value of an array of double */
double mean(int n, double *array)
{
    double sum;
    int i;

    sum = 0;
    // Compute the mean value
    for ( i = 0; i < n; i++)
    {
        sum += array[i];
    }

    sum /= n;
    return( sum );
}

// Return the median value of an array
double median(int n, double *array)
{
    return( median_WIRTH(n, array) );
}

/* Return the bias corrected standard deviation of an array
 *   stddev = sqrt( 1/N-1 * Sum (x - mu)^2 )
 */
double stddev(int n, double *array)
{
    double  stddev, mean;
    int     i;

    stddev = 0;
    mean = 0;
    for ( i = 0; i < n; i++ )
    {
        mean += array[i];
        stddev += array[i] * array[i];
    }
    stddev /= n - 1;
    stddev -= mean * mean / n / (n - 1);
    return( sqrt(stddev) );
}

/* Return the minimum of a list
 */
double min(int n, double *array)
{
    int i;
    double min = DBL_MAX;

    for ( i = 0; i < n; i++ )
        if ( array[i] < min ) min = array[i];

    return min;
}

/* Return the maximum of a list
 */
double max(int n, double *array)
{
    int i;
    double max = DBL_MIN;

    for ( i = 0; i < n; i++ )
        if ( array[i] > max ) max = array[i];

    return max;
}

/* return the sign of a double as a double -1/0./1. */

double  sgn(double x)
{
    if (x < 0)
        return(-1.);
    else if (x > 0)
        return(1.);
    else
        return(0.);
}

/* COMPLEX FUNCTIONS ---------------------------------------------------*/

/*
* make complex number
* Global variables used :
* - none
*/

complex cpx(double re, double im)
{
    complex a;
    a.re = re;
    a.im = im;
    return(a);
}

/*
* square complex --------------------------------------------------
*/

complex sqcpx(complex c)
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
complex acpx(complex c1, complex c2)
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
complex scpx(complex c1, complex c2)
{
    complex res;

    res.re = c1.re - c2.re;
    res.im = c1.im - c2.im;

    return(res);
}

/*
* double substraction to complex ----------------------------------------
*/
complex scpxflt(complex c1, double f2)
{
    complex res;

    res.re = c1.re - f2;
    res.im = c1.im;

    return(res);
}

/*
* double addition to complex ----------------------------------------------
*/
complex acpxflt(complex c1, double f2)
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
complex pcpx(complex c1, complex c2)
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
complex pcpxflt(complex c, double f)
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
complex dcpxflt(complex c, double f)
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
        fprintf(stderr, "dcpxflt: Division by Zero!\n");
        return(c);
    }
}
/* norme complexe au carre------------------------------------------
 * Global variables used :
 * - none
 */
double ncpx(complex c)
{
    double res;

    res = c.re * c.re + c.im * c.im;

    return(res);
}

/* complex inverse ------------------------------------------
 * Global variables used :
 * - none
 */
complex icpx(complex c)
{
    complex res;
    double norm;

    norm = ncpx(c);

    if (norm != 0.)
    {
        res.re = c.re / norm;
        res.im = -c.im / norm;
    }
    else
    {
        fprintf(stderr, "ERROR:(icpx) division by zero\n");
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
complex dcpx(complex c1, complex c2)
{
    complex res;
    double norm;

    norm = ncpx(c2);

    if (norm != 0.)
    {
        res.re = (c1.re * c2.re + c1.im * c2.im) / norm;
        res.im = (c1.im * c2.re - c1.re * c2.im) / norm;
    }
    else
    {
        fprintf(stderr, "ERROR:(dcpx) division by zero\n");
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

complex sqrtcpx(complex c)
{
    complex res;
    double   arg, nc;

    nc = sqrt(ncpx(c));
    arg = atan2(c.im, c.re);

    res.re = sqrt(nc) * cos(arg / 2.);
    res.im = sqrt(nc) * sin(arg / 2.);

    return(res);
}

/* Global variables used :
 * - none
 */
double  sgn_darg(complex z1, complex z2)
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
complex ecpx(complex c)
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
complex coshcpx(complex c)
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
complex sinhcpx(complex c)
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
complex coscpx(complex c)
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
complex sincpx(complex c)
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
complex tancpx(complex c)
{
    return(dcpx(coscpx(c), sincpx(c)));
}

/* log complexe ------------------------------------------
*/
complex lncpx(complex c)
{
    complex res;

    res.re = log(sqrt(c.re * c.re + c.im * c.im));
    res.im = atan2(c.im, c.re);

    return(res);
}

/* END OF COMPLEX FUNCTIONS ----------------------------------------------*/

