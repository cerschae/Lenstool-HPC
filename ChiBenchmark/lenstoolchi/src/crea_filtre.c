#include<stdio.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/* Create a square seeing filter for convolution in d_seeing() function.
 * filtre[y][x] = ln(2)/(PI*seeing^2) * 2^{ -(x*x+y*y)/seeing^2 }  (!!c'est normal??)
 */
void crea_filtre(double seeing, double scale, double **filtre, int n)
{
//  char name[30];
    register int i, j;
    double x, y, s, sepix, nd, norm;

    printf("GAUSSIAN\n");
    sepix = seeing / scale;   // seeing in pixel
    s = (sepix * sepix / (4*log(2.)));
    norm = PI * s;  // PI*R^2 in pixels^2

    nd = ((double) n) / 2.;

    for (i = 0; i < n; i++)
    {
        y = i - nd;
        for (j = 0; j < n; j++)
        {
            x = j - nd;
            filtre[i][j] = exp(-1.*(x * x + y * y) / s) / norm;
        }
    }
//  strcpy(name,"filtre");

    /*
    dim=2;
    size[0]=n;
    size[1]=n;
    strcpy(type,"double");
    strcpy(nature,"real");
    strcpy(mode,"bin");
    sprintf(comments,"filtre seeing");
    wr_ipxs(name,filtre,dim,size,type,mode,nature,comments);
    wrf_fits(name,filtre,n,n,0.,1.,0.,1.);
    */

}

/* Create a square seeing filter for convolution in d_seeing() function.
 * filtre[y][x] = ln(2)/(PI*seeing^2) * 2^{ -(x*x+y*y)/seeing^2 }  (!!c'est normal??)
 */
void crea_filtre_e(double seeing_a, double seeing_b, double seeing_angle, double scale, double **filtre, int n)
{
//  char name[30];
    register int i, j;
    double x, y, s, s1, s2, sepix1, sepix2, nd, norm, ctheta, stheta, s2theta;

    sepix1 = (seeing_a / scale)/(2.0*sqrt(2.0*log(2.0)));   // seeing in pixel
    sepix2 = (seeing_b / scale)/(2.0*sqrt(2.0*log(2.0)));   // seeing in pixel
    s = sepix1 * sepix2;
    s1=2*sepix1*sepix1;
    s2=2*sepix2*sepix2;
    stheta=pow(sin((PI/180.0)*seeing_angle),2.0);
    ctheta=pow(cos((PI/180.0)*seeing_angle),2.0);
    s2theta=sin(2*(PI/180.0)*seeing_angle);

    norm = 2*PI * s;  // 2*PI*sigma1*sigma2 in pixels^2

    nd = ((double) n) / 2.;

    for (i = 0; i < n; i++)
    {
        y = i - nd;
        for (j = 0; j < n; j++)
        {
            x = j - nd;
            filtre[i][j] = exp(-1.*x*x*(ctheta/s1+stheta/s2)-1.*x*y*s2theta*(-1/(s2)+1/(s1))-1.*y*y*(stheta/s1+ctheta/s2)) / norm;
        }
    }
//  strcpy(name,"filtre");

    /*
    dim=2;
    size[0]=n;
    size[1]=n;
    strcpy(type,"double");
    strcpy(nature,"real");
    strcpy(mode,"bin");
    sprintf(comments,"filtre seeing");
    wr_ipxs(name,filtre,dim,size,type,mode,nature,comments);
    wrf_fits(name,filtre,n,n,0.,1.,0.,1.);
    */

}

