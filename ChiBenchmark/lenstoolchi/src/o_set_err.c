#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_set_err           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void  o_set_err(int i, int ipx, double x)
{
    extern struct pot  prec[];

    switch (ipx)
    {
        case(CX):
            prec[i].C.x = x;
            break;
        case(CY):
            prec[i].C.y = x;
            break;
        case(EPOT):
            prec[i].epot = x;
            break;
        case(EMASS):
//          prec[i].epot=x;
            prec[i].emass = x;
            update_epot(i, &prec[i].epot);
            break;
        case(THETA):
            prec[i].theta = x;
            break;
        case(PHI):
            prec[i].phi = x;
            break;
        case(RC):
            prec[i].rc = x;
            break;
        case(B0):
            prec[i].b0 = x;
            break;
        case(ALPHA):
            prec[i].alpha = x;
            break;
        case(BETA):
            prec[i].beta = x;
            break;
        case(RCUT):
            prec[i].rcut = x;
            break;
        case(MASSE):
            prec[i].masse = x;
            break;
        case(ZLENS):
            prec[i].z = x;
            break;
        case(RCSLOPE):
            prec[i].rcslope = x;
            break;
        case(PMASS):
            prec[i].pmass = x;
            break;
        default:
            break;
    }
}


double  o_get_err(int i, int ipx)
{
    extern struct pot  prec[];
    double x;

    switch (ipx)
    {
        case(CX):
            x = prec[i].C.x;
            break;
        case(CY):
            x = prec[i].C.y;
            break;
        case(EPOT):
            x = prec[i].epot;
            break;
        case(EMASS):
            x = prec[i].emass;
            break;
        case(THETA):
            x = prec[i].theta;
            break;
        case(PHI):
            x = prec[i].phi;
            break;
        case(RC):
            x = prec[i].rc;
            break;
        case(B0):
            x = prec[i].b0;
            break;
        case(ALPHA):
            x = prec[i].alpha;
            break;
        case(BETA):
            x = prec[i].beta;
            break;
        case(RCUT):
            x = prec[i].rcut;
            break;
        case(MASSE):
            x = prec[i].masse;
            break;
        case(ZLENS):
            x = prec[i].z;
            break;
        case(RCSLOPE):
            x = prec[i].rcslope;
            break;
        case(PMASS):
            x = prec[i].pmass;
            break;
        default:
            break;
    }

    return x;
}
