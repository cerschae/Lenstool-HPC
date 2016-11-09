#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_set_lens          */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void  o_set_lens(int i, int ipx, double x)
{
    extern struct pot lens[NLMAX];
    extern struct g_cosmo C;

    switch (ipx)
    {
        case(CX):
            lens[i].C.x = x;
            break;
        case(CY):
            lens[i].C.y = x;
            break;
        case(EPOT):
            lens[i].epot = x;
            break;
        case(EMASS):
            lens[i].emass = x;
            if( lens[i].type != 121 )
                update_epot(i, &lens[i].epot);
//          lens[i].epot=x;
//          update_emass(i);
            break;
        case(THETA):
            lens[i].theta = x;
            break;
		case(PHI):
			lens[i].phi = x;
			break;
        case(RC):
            lens[i].rc = x;
            break;
        case(B0):
            lens[i].b0 = x; //6.*pia_c2*x*x;
            break;
        case(ALPHA):
            lens[i].alpha = x;
            break;
        case(BETA):
            lens[i].beta = x;
            break;
        case(RCUT):
            lens[i].rcut = x;
            updatecut(i);
            break;
        case(RCSLOPE):
            lens[i].rcslope = x;
            break;
        case(PMASS):
            lens[i].pmass = x;
            break;
        case(ZLENS):
            lens[i].z = x;
            break;
        case(MASSE):
            lens[i].masse = x;
            break;
        case(OMEGAM):
            C.omegaM = x;
            break;
        case(OMEGAX):
            C.omegaX = x;
            break;
        case(WX):
            C.wX = x;
            break;
        case(WA):
            C.wa = x;
            break;
        default:
            break;
    };
}

void update_epot_ptr(struct pot *ilens, double *epot);
void updatecut_ptr(struct pot *ilens);
void  o_set_lens_ptr(struct pot *ilens, int ipx, double x)
{
    extern struct g_cosmo C;

    switch (ipx)
    {
        case(CX):
            ilens->C.x = x;
            break;
        case(CY):
            ilens->C.y = x;
            break;
        case(EPOT):
            ilens->epot = x;
            break;
        case(EMASS):
            ilens->emass = x;
            if( ilens->type != 121 )
                update_epot_ptr(ilens, &ilens->epot);
//          ilens->epot=x;
//          update_emass(i);
            break;
        case(THETA):
            ilens->theta = x;
            break;
		case(PHI):
			x = ilens->theta;
			break;
        case(RC):
            ilens->rc = x;
            break;
        case(B0):
            ilens->b0 = x; //6.*pia_c2*x*x;
            break;
        case(ALPHA):
            ilens->alpha = x;
            break;
        case(BETA):
            ilens->beta = x;
            break;
        case(RCUT):
            ilens->rcut = x;
            updatecut_ptr(ilens);
            break;
        case(RCSLOPE):
            ilens->rcslope = x;
            break;
        case(PMASS):
            ilens->pmass = x;
            break;
        case(ZLENS):
            ilens->z = x;
            break;
        case(MASSE):
            ilens->masse = x;
            break;
        case(OMEGAM):
            C.omegaM = x;
            break;
        case(OMEGAX):
            C.omegaX = x;
            break;
        case(WX):
            C.wX = x;
            break;
        case(WA):
            C.wa = x;
            break;
        default:
            break;
    }
}


double  o_get_lens(int i, int ipx)
{
    extern struct pot lens[NLMAX];
    extern struct g_cosmo C;
    double x;

    switch (ipx)
    {
        case(CX):
            x = lens[i].C.x;
            break;
        case(CY):
            x = lens[i].C.y;
            break;
        case(EPOT):
            x = lens[i].epot;
            break;
        case(EMASS):
            //x=lens[i].epot;
            x = lens[i].emass;
            break;
        case(THETA):
            x = lens[i].theta;
            break;
        case(PHI):
            x = lens[i].phi;
            break;
        case(RC):
            x = lens[i].rc;
            break;
        case(B0):
            x = lens[i].b0;
            break;
        case(ALPHA):
            x = lens[i].alpha;
            break;
        case(BETA):
            x = lens[i].beta;
            break;
        case(RCUT):
            x = lens[i].rcut;
            break;
        case(RCSLOPE):             
            x = lens[i].rcslope;
            break;
        case(PMASS):
            x = lens[i].pmass;
            break;
        case(ZLENS):
            x = lens[i].z;
            break;
        case(MASSE):
            x = lens[i].masse;
            break;
        case(OMEGAM):
            x = C.omegaM;
            break;
        case(OMEGAX):
            x = C.omegaX;
            break;
        case(WX):
            x = C.wX;
            break;
        case(WA):
            x = C.wa;
            break;
        default:
            break;
    }

    return x;
}
