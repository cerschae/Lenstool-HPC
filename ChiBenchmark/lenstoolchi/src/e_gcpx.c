#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        e_grad              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

/*
* SIEMD KK
* Global variables used :
* - none
*/

complex csiemd(double x, double y, double eps, double b0)
{
    double  sqe, ce0, ce1, ce2;
    double  r;
    complex zres;

    r = sqrt(x * x + y * y);
    sqe = sqrt(eps);
    ce0 = b0 * (1 - eps * eps) / 2. / sqe;
    ce1 = 2.*sqe / (1. + eps);
    ce2 = 2.*sqe / (1. - eps);
    zres.re = ce0 * asin(ce1 * x / r);
    zres.im = ce0 * asinh(ce2 * y / r);

    return(zres);
}

/*
* I*w,v=0.5 Kassiola & Kovner, 1993 PIEMD, paragraph 4.1
* Same as ci05() but with circular clumps. Much faster!!
*
*
* Global variables used :
* - none
*/

/*
* I*w,v=0.5 Kassiola & Kovner, 1993 PIEMD, paragraph 4.1
*
* Global variables used :
* - none
*/

complex ci05(double x, double y, double eps, double rc)
{
    double  sqe, cx1, cxro, cyro, rem2;
    complex zci, znum, zden, zis, zres;
    double norm;

    sqe = sqrt(eps);
    cx1 = (1. - eps) / (1. + eps);
    cxro = (1. + eps) * (1. + eps);
    cyro = (1. - eps) * (1. - eps);
    rem2 = x * x / cxro + y * y / cyro;
    /*zci=cpx(0.,-0.5*(1.-eps*eps)/sqe);
    znum=cpx(cx1*x,(2.*sqe*sqrt(rc*rc+rem2)-y/cx1));
    zden=cpx(x,(2.*rc*sqe-y));
    zis=pcpx(zci,lncpx(dcpx(znum,zden)));
    zres=pcpxflt(zis,b0);*/

    // --> optimized code
    zci.re = 0;
    zci.im = -0.5 * (1. - eps * eps) / sqe;
    znum.re = cx1 * x;
    znum.im = 2.*sqe * sqrt(rc * rc + rem2) - y / cx1;
    zden.re = x;
    zden.im = 2.*rc * sqe - y;
    norm = zden.re * zden.re + zden.im * zden.im;     // zis = znum/zden
    zis.re = (znum.re * zden.re + znum.im * zden.im) / norm;
    zis.im = (znum.im * zden.re - znum.re * zden.im) / norm;
    norm = zis.re;
    zis.re = log(sqrt(norm * norm + zis.im * zis.im));  // ln(zis) = ln(|zis|)+i.Arg(zis)
    zis.im = atan2(zis.im, norm);
//  norm = zis.re;
    zres.re = zci.re * zis.re - zci.im * zis.im;   // Re( zci*ln(zis) )
    zres.im = zci.im * zis.re + zis.im * zci.re;   // Im( zci*ln(zis) )
    //zres.re = zis.re*b0;
    //zres.im = zis.im*b0;

    return(zres);
}


/*
* I1.5 KK
* Global variables used :
* - none
*/

complex ci15(double x, double y, double eps, double rc, double b0)
{
    double  sqe, cx1, cxro, cyro, wrem;
    complex zaa, zbb, ztot, zres;

    cx1 = (1. - eps) / (1. + eps);
    cxro = (1. + eps) * (1. + eps);
    cyro = (1. - eps) * (1. - eps);
    sqe = sqrt(eps);
    wrem = sqrt(rc * rc + x * x / cxro + y * y / cyro);
    zaa = cpx(x, -y + 2 * rc * sqe);
    zbb = cpx(cx1 * x, -y / cx1 + 2 * wrem * sqe);
    ztot = scpx(icpx(zaa), pcpxflt(icpx(zbb), rc / wrem));

    zres = pcpxflt(ztot, b0 * rc * (1. - eps * eps));
    return(zres);
}

/*
* I1.0 KK -we actually return I_1/rc
* Global variables used :
* - none
*/

complex ci10(double x, double y, double eps, double rc, double b0)
{
    double  cx1, cxro, cyro, wrem;
    complex sqzs, sqz, zbar, zaa, zbb, ztot, zres;

    cx1 = (1. - eps) / (1. + eps);
    cxro = (1. + eps) * (1. + eps);
    cyro = (1. - eps) * (1. - eps);
    wrem = sqrt(rc * rc + x * x / cxro + y * y / cyro);
    zbar = cpx(x, -y);
    sqz = sqrtcpx(acpxflt(pcpx(zbar, zbar), (4.*eps * rc * rc)));
    sqzs = pcpxflt(sqz, sgn_darg(sqz, zbar));
    zaa = acpx(zbar, sqzs);
    zbb = acpx(cpx(cx1 * x, -y / cx1), sqzs);
    ztot = dcpx(lncpx(pcpxflt(dcpx(zaa, zbb), wrem / rc)), sqzs);

    zres = pcpxflt(ztot, b0 * (1. - eps * eps));
    return(zres);
}
