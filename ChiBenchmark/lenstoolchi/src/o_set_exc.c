#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_set_exc           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/
/* Check each parameter for lens i and decide if the optimization
 * has to be continued */

void o_set_exc( int i,
                double excu[NLMAX][NPAMAX], double excd[NLMAX][NPAMAX],
                int block[NLMAX][NPAMAX])
{
    extern struct pot      lens[], lmin[], lmax[], prec[];
    extern struct ipot     ip;
    register int j;

    for (j = 0; j < ip.pmax; j++)
    {
        excu[i][j] = .3;
        excd[i][j] = .3;
    }

    if ((block[i][CX] == 0) || (lmax[i].C.x - lmin[i].C.x <= prec[i].C.x))
    {
        block[i][CX] = 0;
        lmax[i].C.x = lmin[i].C.x = lens[i].C.x;
        prec[i].C.x = 1.;
        excd[i][CX] = 0.;
        excu[i][CX] = 0.;
    }
    if ((block[i][CY] == 0) || (lmax[i].C.y - lmin[i].C.y <= prec[i].C.y))
    {
        block[i][CY] = 0;
        lmax[i].C.y = lmin[i].C.y = lens[i].C.y;
        prec[i].C.y = 1.;
        excu[i][CY] = 0.;
        excd[i][CY] = 0.;
    }
    if ((block[i][EPOT] == 0) || (lmax[i].epot - lmin[i].epot <= prec[i].epot))
    {
        block[i][EPOT] = 0;
        lmax[i].epot = lmin[i].epot = lens[i].epot;
        prec[i].epot = 1.;
        excu[i][EPOT] = 0.;
        excd[i][EPOT] = 0.;
    }
    if ((block[i][EMASS] == 0) || (lmax[i].emass - lmin[i].emass <= prec[i].emass))
    {
        block[i][EMASS] = 0;
        lmax[i].emass = lmin[i].emass = lens[i].emass;
        prec[i].emass = 1.;
        excu[i][EMASS] = 0.;
        excd[i][EMASS] = 0.;
    }
    if ((block[i][THETA] == 0) || (lmax[i].theta - lmin[i].theta <= prec[i].theta))
    {
        block[i][THETA] = 0;
        lmax[i].theta = lmin[i].theta = lens[i].theta;
        prec[i].theta = 1.;
        excu[i][THETA] = 0.;
        excd[i][THETA] = 0.;
    }
    if ((block[i][RC] == 0) || (lmax[i].rc - lmin[i].rc <= prec[i].rc))
    {
        block[i][RC] = 0;
        lmax[i].rc = lmin[i].rc = lens[i].rc;
        prec[i].rc = 1.;
        excu[i][RC] = 0.;
        excd[i][RC] = 0.;
    }
    if ((block[i][B0] == 0) || (lmax[i].b0 - lmin[i].b0 <= prec[i].b0))
    {
        block[i][B0] = 0;
        lmax[i].b0 = lmin[i].b0 = lens[i].b0;
        prec[i].b0 = 1.;
        excu[i][B0] = 0.;
    }
    if ((block[i][ALPHA] == 0) || (lmax[i].alpha - lmin[i].alpha <= prec[i].alpha))
    {
        block[i][ALPHA] = 0;
        lmax[i].alpha = lmin[i].alpha = lens[i].alpha;
        prec[i].alpha = 1.;
        excu[i][ALPHA] = 0.;
        excd[i][ALPHA] = 0.;
    }
    if ((block[i][BETA] == 0) || (lmax[i].beta - lmin[i].beta <= prec[i].beta))
    {
        block[i][BETA] = 0;
        lmax[i].beta = lmin[i].beta = lens[i].beta;
        prec[i].beta = 1.;
        excu[i][BETA] = 0.;
    }
    if ((block[i][RCUT] == 0) || (lmax[i].rcut - lmin[i].rcut <= prec[i].rcut))
    {
        block[i][RCUT] = 0;
        lmax[i].rcut = lmin[i].rcut = lens[i].rcut;
        prec[i].rcut = 1.;
        excu[i][RCUT] = 0.;
    }
    if ((block[i][MASSE] == 0) || (lmax[i].masse - lmin[i].masse <= prec[i].masse))
    {
        block[i][MASSE] = 0;
        lmax[i].masse = lmin[i].masse = lens[i].masse;
        prec[i].masse = 1.;
        excu[i][MASSE] = 0.;
    }
    if ((block[i][PMASS] == 0) || (lmax[i].pmass - lmin[i].pmass <= prec[i].pmass))
    {
        block[i][PMASS] = 0;
        lmax[i].pmass = lmin[i].pmass = lens[i].pmass;
        prec[i].pmass = 1.;
        excu[i][PMASS] = 0.;
    }
	if ((block[i][RCSLOPE] == 0) || (lmax[i].rcslope - lmin[i].rcslope <= prec[i].rcslope))
    {
        block[i][RCSLOPE] = 0;
        lmax[i].rcslope = lmin[i].rcslope = lens[i].rcslope;
        prec[i].rcslope = 1.;
        excu[i][RCSLOPE] = 0.;
    }	
}
