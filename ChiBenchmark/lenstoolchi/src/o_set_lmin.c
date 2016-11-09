#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_set_lmin          */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void  o_set_lmin(int i, int ipx, double x)
{
    extern struct pot lmin[];
    extern struct g_cosmo clmin;
    extern struct galaxie source[NFMAX];
    extern struct galaxie smin[NFMAX];
    extern struct vfield vfmin;

    switch (ipx)
    {
        case(CX):
            lmin[i].C.x = x;
            break;
        case(CY):
            lmin[i].C.y = x;
            break;
        case(EPOT):
            lmin[i].epot = x;
            break;
        case(EMASS):
            //lmin[i].epot=x;
            lmin[i].emass = x;
            update_epot(i, &lmin[i].epot);
            break;
        case(THETA):
            lmin[i].theta = x;
            break;
		case(PHI):
			lmin[i].phi = x;
			break;
        case(RC):
            lmin[i].rc = x;
            break;
        case(B0):
            lmin[i].b0 = x;
            break;
        case(ALPHA):
            lmin[i].alpha = x;
            break;
        case(BETA):
            lmin[i].beta = x;
            break;
        case(RCUT):
            lmin[i].rcut = x;
            break;
        case(MASSE):
            lmin[i].masse = x;
            break;
        case(ZLENS):
            lmin[i].z = x;
            break;
        case(RCSLOPE):
            lmin[i].rcslope = x;
            break;
        case(PMASS):
            lmin[i].pmass = x;
            break;
        case(OMEGAM):
            clmin.omegaM = x;
            break;
        case(OMEGAX):
            clmin.omegaX = x;
            break;
        case(WX):
            clmin.wX = x;
            break;
        case(WA):
            clmin.wa = x;
            break;
        case(SCX):
            smin[i].C.x = x;
            break;
        case(SCY):
            smin[i].C.y = x;
            break;
        case(SA):
            smin[i].E.a = x;
            break;
        case(SB):
            smin[i].E.b = x;
            break;
        case(STHETA):
            smin[i].E.theta = x;
            break;
        case(SFLUX):
            smin[i].mag = x;
            break;
        case(VFCX):
            vfmin.C.x = x;
            break;
        case(VFCY):
            vfmin.C.y = x;
            break;
        case(VFVT):
            vfmin.vt = x;
            break;
        case(VFRT):
            vfmin.rt = x;
            break;
        case(VFI):
            vfmin.i = x;
            break;
        case(VFTHETA):
            vfmin.theta = x;
            break;
        case(VFLCENT):
            vfmin.lcent = x;
            break;
        case(VFSIGMA):
            vfmin.sigma = x;
            break;
        default:
            break;
    }
}

double o_get_lmin(int i, int ipx)
{
    extern struct pot lmin[];
    extern struct g_cosmo clmin;
    extern struct vfield vfmin;
    double x;

    switch (ipx)
    {
        case(CX):
            x = lmin[i].C.x;
            break;
        case(CY):
            x = lmin[i].C.y;
            break;
        case(EPOT):
            x = lmin[i].epot;
            break;
        case(EMASS):
//          x=lmin[i].epot;
            x = lmin[i].emass;
            break;
        case(THETA):
            x = lmin[i].theta;
            break;
		case(PHI):
			x = lmin[i].phi;
			break;
        case(RC):
            x = lmin[i].rc;
            break;
        case(B0):
            x = lmin[i].b0;
            break;
        case(ALPHA):
            x = lmin[i].alpha;
            break;
        case(BETA):
            x = lmin[i].beta;
            break;
        case(RCUT):
            x = lmin[i].rcut;
            break;
        case(MASSE):
            x = lmin[i].masse;
            break;
        case(ZLENS):
            x = lmin[i].z;
            break;
        case(RCSLOPE):
            x = lmin[i].rcslope;
            break;
        case(PMASS):
            x = lmin[i].pmass;
            break;
        case(OMEGAM):
            x = clmin.omegaM;
            break;
        case(OMEGAX):
            x = clmin.omegaX;
            break;
        case(WX):
            x = clmin.wX;
            break;
        case(WA):
            x = clmin.wa;
            break;
        case(VFCX):
            x = vfmin.C.x;
            break;
        case(VFCY):
            x = vfmin.C.y;
            break;
        case(VFVT):
            x = vfmin.vt;
            break;
        case(VFRT):
            x = vfmin.rt;
            break;
        case(VFI):
            x = vfmin.i;
            break;
        case(VFTHETA):
            x = vfmin.theta;
            break;
        case(VFLCENT):
            x = vfmin.lcent;
            break;
        case(VFSIGMA):
            x = vfmin.sigma;
            break;
        default:
            break;
    }

    return x;
}

