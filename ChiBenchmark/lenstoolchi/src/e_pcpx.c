#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/*
file e_pcpx.c

JP Kneib
25 Avril 97
OMP, Toulouse

content: function to compute the projevcted gravitational potentiel
using the complex PIEMD formulas (cf Hjorth & Kneib 97)
*/

/*
* SIEMD KK
*/

double  psiemd(double x, double y, double eps, double b0)
{
    double  sqe, ci, r, cx, cy, cxe, cye;

    sqe = 2.*sqrt(eps);
    ci = (1. - eps * eps) / sqe;
    r = sqrt(x * x + y * y);
    cx = x / r;
    cy = y / r;
    cxe = sqe / (1. + eps) * cx;
    cye = sqe / (1. - eps) * cy;

    return( b0*ci*( cx*asin(cxe) + cy*asinh(cye) ) );
}


/*
* I0.5 KK
*/

double  pi05(double x, double y, double eps, double rc, double b0)
{
    double  sqe, ci, cxro, cyro, rem2, e1, e2, z;
    complex eta, zeta, b1, b2, a1, a2, c1, c2, ckk;


    sqe = sqrt(eps);
    ci = .5 * (1. - eps * eps) / sqe;
    cxro = (1. + eps) * (1. + eps);
    cyro = (1. - eps) * (1. - eps);
    rem2 = x * x / cxro + y * y / cyro;
    e1 = 2.*sqe / (1 - eps);
    e2 = 2.*sqe / (1 + eps);
    z = sqrt(x * x + y * y);
    eta = cpx(-.5 * asinh(e1 * y / z), .5 * asin(e2 * x / z));
    zeta = cpx( 0.5 * log( (sqrt(rem2) + sqrt(rc * rc + rem2)) / rc), 0. );
    b1 = coshcpx(acpx(eta, zeta));
    b2 = coshcpx(scpx(eta, zeta));
    a1 = lncpx( dcpx(sqcpx(coshcpx(eta)), pcpx(b1, b2)) );
    a2 = lncpx(dcpx(b1, b2));
    c1 = pcpx(sinhcpx(pcpxflt(eta, 2.)), a1);
    c2 = pcpx(sinhcpx(pcpxflt(zeta, 2.)), a2);
    ckk = acpx(c1, c2);

    return( b0*ci*rc / sqrt(rem2)*(ckk.im*x - ckk.re*y) );
}


/*
* I1.5 KK
*/

double  pi15(double x, double y, double eps, double rc, double b0)
{
    double  sqe, ci, cxro, cyro, rem2, e1, e2, z;
    complex eta, zeta, b1, b2, a1, a2, c1, c2, ckk;
    complex s1, s2, t1, t2, u1, u2, c3;


    sqe = sqrt(eps);
    ci = .5 * (1. - eps * eps) / sqe;
    cxro = (1. + eps) * (1. + eps);
    cyro = (1. - eps) * (1. - eps);
    rem2 = x * x / cxro + y * y / cyro;
    e1 = 2.*sqe / (1 - eps);
    e2 = 2.*sqe / (1 + eps);
    z = sqrt(x * x + y * y);
    eta = cpx(-.5 * asinh(e1 * y / z), .5 * asin(e2 * x / z));
    zeta = cpx( 0.5 * log( (sqrt(rem2) + sqrt(rc * rc + rem2)) / rc), 0. );
    b1 = coshcpx(acpx(eta, zeta));
    b2 = coshcpx(scpx(eta, zeta));
    a1 = lncpx( dcpx(sqcpx(coshcpx(eta)), pcpx(b1, b2)) );
    a2 = lncpx(dcpx(b1, b2));
    s1 = sinhcpx(pcpxflt(eta, 2.));
    s2 = sinhcpx(pcpxflt(zeta, 2.));
    t1 = tancpx(acpx(eta, zeta));
    t2 = tancpx(scpx(eta, zeta));
    u1 = pcpx(t1, scpx(s2, s1));
    u2 = pcpx(t2, acpx(s2, s1));
    c1 = pcpx(s1, a1);
    c2 = pcpxflt(pcpx(s2, a2), 2.);
    c3 = pcpxflt(acpx(u1, u2), 0.5);
    ckk = acpx(acpx(c1, c2), c3);

    return( -b0*ci*rc / sqrt(rem2)*(ckk.im*x - ckk.re*y) );
}

/*
* I1.0 KK -we actually return I_1/rc
*/

double  pi10(double x, double y, double eps, double rc, double b0)
{
//double    cx1,cxro,cyro,wrem,sgn_darg();
//complex   sqzs,sqz,zbar,zaa,zbb,ztot,zres;

    return(0.);
}


/*
END
*/

