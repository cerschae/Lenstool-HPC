#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        e_grad2             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

/*
* SIEMD KK
* Global variables used :
* - none
*/
void mdcsiemd(double x, double y, double eps, double b0, struct matrix *res)
{
    double  cx, cy;
    double  theta;
    struct  matrix Q;

    theta = atan2(y, x);
    cx = x / (1 + eps);
    cy = y / (1 - eps);
    Q.a = Q.b = Q.d = 0.;
    Q.c = b0 / sqrt(cx * cx + cy * cy);
    *res = rotmatrix(&Q, theta);
}

/*
* derivates of I0.5 KK
* Parameters :
* - (x,y) is the computation position of the potential
* - eps is the ellepticity (a-b)/(a+b)
* - rc is the core radius
* - b0 asymptotic Einstein radius E0. (6pi*vdisp^2/c^2)
*
* Return a the 4 second derivatives of the PIEMD potential
*/
void mdci05(double x, double y, double eps, double rc, double b0, struct matrix *res)
{
    double   ci, sqe, cx1, cxro, cyro, wrem;
    double  didyre, didyim, didxre;// didxim;
    double  cx1inv, den1, num2, den2;
//  complex znum,zden,zdidx,zdidy;
//  struct  matrix  res;

    sqe = sqrt(eps);
    cx1 = (1. - eps) / (1. + eps);
    cx1inv = 1. / cx1;
    cxro = (1. + eps) * (1. + eps);     /* rem^2=x^2/(1+e^2) + y^2/(1-e^2) Eq 2.3.6*/
    cyro = (1. - eps) * (1. - eps);
    ci = 0.5 * (1. - eps * eps) / sqe;
    wrem = sqrt(rc * rc + x * x / cxro + y * y / cyro); /*wrem^2=w^2+rem^2 with w core radius*/
    /*
        zden=cpx(x,(2.*rc*sqe-y)); //denominator
        znum=cpx(cx1*x,(2.*sqe*wrem-y/cx1)); // numerator

        zdidx=acpx(dcpx(cpx(2.*ci*sqe*x/cxro/wrem,-cx1*ci),znum),
                   dcpx(cpx(0.,ci),zden)); // dI/dx with I in Eq 4.1.2
        zdidy=acpx(dcpx(cpx(-ci/cx1+2.*ci*sqe*y/cyro/wrem,0.),znum),
                   dcpx(cpx(ci,0.),zden)); // dI/dy with I in Eq 4.1.2
        //in Eq 4.1.2 I=A*ln(u/v)) ==> dI/dx=A*(u'/u-1/v) because v'=1
        res.a=b0*zdidx.re;
        res.b=res.d=b0*(zdidy.re+zdidx.im)/2.;
        res.c=b0*zdidy.im;
    */
    den1 = 2.*sqe * wrem - y * cx1inv;
    den1 = cx1 * cx1 * x * x + den1 * den1;
    num2 = 2.*rc * sqe - y;
    den2 = x * x + num2 * num2;

    didxre = ci * ( cx1 * (2.*sqe * x * x / cxro / wrem - 2.*sqe * wrem + y * cx1inv) / den1 + num2 / den2 );

//  didxim = ci * ( (2*sqe*x*y*cx1inv/cxro/wrem - cx1*cx1*x - 4*eps*x/cxro)/den1 + x/den2 );
    didyre = ci * ( (2 * sqe * x * y * cx1 / cyro / wrem - x) / den1 + x / den2 );

    didyim = ci * ( (2 * sqe * wrem * cx1inv - y * cx1inv * cx1inv - 4 * eps * y / cyro +
                     2 * sqe * y * y / cyro / wrem * cx1inv) / den1 - num2 / den2 );

    res->a = b0 * didxre;
    res->b = res->d = b0 * didyre; //(didyre+didxim)/2.;
    res->c = b0 * didyim;

//  return(res);


}

/*
* derivates of I1.5 KK
*/

struct matrix mdci15(double x, double y, double eps, double rc, double b0)
{
    double   sqe, cx1, cxro, cyro, wrem2, wrem;
    complex zaa, zbb, zcc, zdd, zee, zff, zdidx, zdidy;
    struct  matrix  res;

    sqe = sqrt(eps);
    cx1 = (1. - eps) / (1. + eps);
    cxro = (1. + eps)*(1. + eps);
    cyro = (1. - eps)*(1. - eps);
    wrem2 = rc*rc + x*x / cxro + y*y / cyro;
    wrem = sqrt(wrem2);

    zaa = cpx(x, -y + 2*rc*sqe);
    zbb = cpx(cx1*x, -y / cx1 + 2*wrem*sqe);
    zcc = pcpxflt(icpx(sqcpx(zaa)), eps*eps - 1.);
    zdd = pcpxflt(icpx(sqcpx(zbb)), rc / wrem2 / wrem);
    zee = acpxflt(pcpxflt( cpx(cx1*x, -y / cx1 + 4.*wrem*sqe), cx1*x), cyro*wrem2);
    zff = acpx(pcpxflt( cpx(cx1*x, -y / cx1 + 4.*wrem*sqe), y / cx1), cpx(0., -cxro*wrem2));
    zdidx = acpx(zcc, pcpx(zdd, zee));
    zdidy = acpx(pcpx(zcc, cpx(0., -1.)), pcpx(zdd, zff));

    res.a = b0*rc*zdidx.re;
    res.b = res.d = b0*rc*(zdidy.re + zdidx.im) / 2.;
    res.c = b0*rc*zdidy.im;
    return(res);
}

/*
* derivates of I1.0 KK
*/

struct matrix mdci10(double x, double y, double eps, double rc, double b0)
{
    double   t, cx1, cxro, cyro, wrem2, wrem; //,sqe
    complex zz, zeps, zbar, ztot, zcc, zdd, zdidx, zdidy;
    struct  matrix  res;

    cx1 = (1. - eps) / (1. + eps);
    cxro = (1. + eps)*(1. + eps);
    cyro = (1. - eps)*(1. - eps);
    wrem2 = rc*rc + x*x / cxro + y*y / cyro;
    wrem = sqrt(wrem2);

    zbar = cpx(x, -y);
    zeps = cpx(cx1*x, -y / cx1);
    zz = acpxflt(sqcpx(zbar), 4*eps*rc*rc);
    ztot = pcpx(zbar, ci10(x, y, eps, rc, b0));
    zcc = acpxflt(pcpxflt(zeps, x / cxro / wrem2),
                  2.*eps / (1. + eps));
    zdd = acpx(pcpxflt(zeps, y / cyro / wrem2),
               cpx(0., 2.*eps / (1. - eps)) );
    t = b0*(1. - eps*eps);
    zdidx = dcpx(scpx(pcpxflt(zcc, t), ztot), zz);
    zdidy = dcpx(acpx(pcpxflt(zdd, t), pcpx(ztot, cpx(0., 1.))), zz);

    res.a = zdidx.re;
    res.b = res.d = (zdidy.re + zdidx.im) / 2.;
    res.c = zdidy.im;
    return(res);
}

