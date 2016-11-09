#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_set_lmax          */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void o_set_lmax(int i, int ipx, double x)
{
    extern struct pot lmax[];
    extern struct g_cosmo clmax;
    extern struct galaxie source[NFMAX];
    extern struct galaxie smax[NFMAX];
    extern struct vfield vfmax;

    switch (ipx)
    {
        case(CX):
            lmax[i].C.x = x;
            break;
        case(CY):
            lmax[i].C.y = x;
            break;
        case(EPOT):
            lmax[i].epot = x;
            break;
        case(EMASS):
            lmax[i].emass = x;
            update_epot(i, &lmax[i].epot);
            break;
        case(THETA):
            lmax[i].theta = x;
            break;
		case(PHI):
			lmax[i].phi = x;
			break;
        case(RC):
            lmax[i].rc = x;
            break;
        case(B0):
            lmax[i].b0 = x;
            break;
        case(ALPHA):
            lmax[i].alpha = x;
            break;
        case(BETA):
            lmax[i].beta = x;
            break;
        case(RCUT):
            lmax[i].rcut = x;
            break;
        case(MASSE):
            lmax[i].masse = x;
            break;
        case(ZLENS):
            lmax[i].z = x;
            break;
        case(RCSLOPE):
            lmax[i].rcslope = x;
            break;
        case(PMASS):
            lmax[i].pmass = x;
            break;
        case(OMEGAM):
            clmax.omegaM = x;
            break;
        case(OMEGAX):
            clmax.omegaX = x;
            break;
        case(WX):
            clmax.wX = x;
            break;
        case(WA):
            clmax.wa = x;
            break;
        case(SCX):
            smax[i].C.x = x;
            break;
        case(SCY):
            smax[i].C.y = x;
            break;
        case(SA):
            smax[i].E.a = x;
            break;
        case(SB):
            smax[i].E.b = x;
            break;
        case(STHETA):
            smax[i].E.theta = x;
            break;
        case(SFLUX):
            smax[i].mag = x;
            break;
        case(VFCX):
            vfmax.C.x = x;
            break;
        case(VFCY):
            vfmax.C.y = x;
            break;
        case(VFVT):
            vfmax.vt = x;
            break;
        case(VFRT):
            vfmax.rt = x;
            break;
        case(VFI):
            vfmax.i = x;
            break;
        case(VFTHETA):
            vfmax.theta = x;
            break;
        case(VFLCENT):
            vfmax.lcent = x;
            break;
        case(VFSIGMA):
            vfmax.sigma = x;
            break;
        default:
            break;
    }
}

double o_get_lmax(int i, int ipx)
{
    extern struct pot lmax[];
    extern struct g_cosmo clmax;
    extern struct vfield vfmax;
    double x;

    switch (ipx)
    {
        case(CX):
            x = lmax[i].C.x;
            break;
        case(CY):
            x = lmax[i].C.y;
            break;
        case(EPOT):
            x = lmax[i].epot;
            break;
        case(EMASS):
            //x=lmax[i].epot;
            x = lmax[i].emass;
            break;
        case(THETA):
            x = lmax[i].theta;
            break;
		case(PHI):
			x = lmax[i].phi;
			break;
        case(RC):
            x = lmax[i].rc;
            break;
        case(B0):
            x = lmax[i].b0;
            break;
        case(ALPHA):
            x = lmax[i].alpha;
            break;
        case(BETA):
            x = lmax[i].beta;
            break;
        case(RCUT):
            x = lmax[i].rcut;
            break;
        case(MASSE):
            x = lmax[i].masse;
            break;
        case(ZLENS):
            x = lmax[i].z;
            break;
        case(RCSLOPE):
            x = lmax[i].rcslope;
            break;
        case(PMASS):
            x = lmax[i].pmass;
            break;
        case(OMEGAM):
            x = clmax.omegaM;
            break;
        case(OMEGAX):
            x = clmax.omegaX;
            break;
        case(WX):
            x = clmax.wX;
            break;
        case(WA):
            x = clmax.wa;
            break;
        case(VFCX):
            x = vfmax.C.x;
            break;
        case(VFCY):
            x = vfmax.C.y;
            break;
        case(VFVT):
            x = vfmax.vt;
            break;
        case(VFRT):
            x = vfmax.rt;
            break;
        case(VFI):
            x = vfmax.i;
            break;
        case(VFTHETA):
            x = vfmax.theta;
            break;
        case(VFLCENT):
            x = vfmax.lcent;
            break;
        case(VFSIGMA):
            x = vfmax.sigma;
            break;
        default:
            break;
    }

    return x;
}


