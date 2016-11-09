#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      programme   d_profil.c          */
/*      auteur      Henri Bonnet            */
/*      place       OMP             */
/*      date        12.11.1991          */
/*      version     1               */
/****************************************************************
 *  d_profil(x,y,k)
 *      double x,y,k;
 *
 *  evalue, en (x,y), la fonction:
 *      f(x,y)=Io/(1+alpha**2*(e**2*(x-C.x)**2+(y-C.y)**2))
 *  k est l'indice de l'objet qu'on integre
 *
 *  les parametres sont definis dans "para"
 *      alpha=1/(e*l)
 */

double d_profil(double x, double y, const struct galaxie *gal)
{
    double xx, yy, xxx, yyy;
    double res;
    const extern  struct g_observ O;
    const extern  struct g_large L;

    res = 0;

    if (gal->c == 's')
    {
        xx = x - gal->C.x;
        yy = y - gal->C.y;
        xxx = xx * xx + yy * yy;
        res = 100 * exp((-1.) * xxx / O.r0st);
    }
    else if (L.vitesse == 0 || gal->type == 0 )
    {
        xx = ((x - gal->C.x) * sin(gal->E.theta) - (y - gal->C.y) * cos(gal->E.theta)) / gal->E.b;
        yy = ((y - gal->C.y) * sin(gal->E.theta) + (x - gal->C.x) * cos(gal->E.theta)) / gal->E.a;
        res = 1. + (xx * xx + yy * yy); 
        if (gal->mag > 0)
            res = pow(10., (26. - gal->mag) / 2.5) / res;
        else
            res = gal->I0 / res;
    }
    else if (L.vitesse == 1 || gal->type == 1 )
    {
        xx = ((x - gal->C.x) * sin(gal->E.theta) - (y - gal->C.y) * cos(gal->E.theta)) / gal->E.b;
        yy = ((y - gal->C.y) * sin(gal->E.theta) + (x - gal->C.x) * cos(gal->E.theta)) / gal->E.a;
        if (xx*xx + yy*yy < 3.)
        {
            res = gal->I0 * sqrt(fabs(yy));
            if (yy < 0.)
                res = -res;
        }
    }
	 /*  Exponential Disk Profile  */
	 else if (L.vitesse == 2 || gal->type == 2 )
	 {
	     xx = ((x - gal->C.x) * sin(gal->E.theta) - (y - gal->C.y) * cos(gal->E.theta)) / gal->E.b;
        yy = ((y - gal->C.y) * sin(gal->E.theta) + (x - gal->C.x) * cos(gal->E.theta)) / gal->E.a;
        res =exp(sqrt(xx * xx + yy * yy)); 
        if (gal->mag > 0)
            res = pow(10., (26. - gal->mag) / 2.5) / res;
        else
            res = gal->I0 / res;
	 }
	 /*  2D Gaussian  */
	 else if (L.vitesse == 3 || gal->type == 3 ) 
	 {
	     xx = ((x - gal->C.x) * sin(gal->E.theta) - (y - gal->C.y) * cos(gal->E.theta)) / gal->E.b;
        yy = ((y - gal->C.y) * sin(gal->E.theta) + (x - gal->C.x) * cos(gal->E.theta)) / gal->E.a;
        res = exp(0.5*(xx * xx + yy * yy)); 
        if (gal->mag > 0)
        {
            //res = pow(10., (26. - gal->mag) / 2.5) / res;
     	    //mag AB Oke 1983 (B. Clement 2011)
     	    //magnitude given in catalogue corresponds to the peak of the gaussian in source plan
     	    //TODO: should be the integral
     	    //flux in image is uJy
            res = pow(10.,(-0.4*(gal->mag+48.57)))/1E-23/1E-6 / res;
        }
        else
            res = gal->I0 / res;
	 }
	 /*  Sersic  */
	 else if (L.vitesse == 4 || gal->type == 4 ) 
	 {
	     xx = ((x - gal->C.x) * sin(gal->E.theta) - (y - gal->C.y) * cos(gal->E.theta)) / gal->E.b;
        yy = ((y - gal->C.y) * sin(gal->E.theta) + (x - gal->C.x) * cos(gal->E.theta)) / gal->E.a;
        res = exp(pow(xx * xx + yy * yy, 0.5/gal->var1)); 
        if (gal->mag > 0)
            res = pow(10., (26. - gal->mag) / 2.5) / res; 
        else
            res = gal->I0 / res;
	 }
     // Uniform disk at gal->I0 within radius R
     else if (L.vitesse == 5 || gal->type == 5)
     {
	     xx = ((x - gal->C.x) * sin(gal->E.theta) - (y - gal->C.y) * cos(gal->E.theta)) / gal->E.b;
        yy = ((y - gal->C.y) * sin(gal->E.theta) + (x - gal->C.x) * cos(gal->E.theta)) / gal->E.a;
        if ( xx + yy < gal->E.a * gal->E.b )
            res = pow(10., (26. - gal->mag) / 2.5);
        else
            res = 0.;
     }
     else
     {
         fprintf(stderr, "ERROR: source %s brightness profil type %d unknown\n", gal->n, gal->type);
         exit(1);
     }


	 return(res);

}
