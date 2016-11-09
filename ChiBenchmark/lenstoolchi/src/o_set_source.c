#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_set_source          */
/*      auteur:     Eric Jullo         */
/*      date:       3/10/2011            */
/*      place:      Marseille            */
/****************************************************************/

void  o_set_source(struct galaxie *source, int ipx, double val)
{

    switch (ipx)
    {
        case(SCX):
            source->C.x = val;
            break;
        case(SCY):
            source->C.y = val;
            break;
        case(SA):
            source->E.a = val;
            break;
        case(SB):
            source->E.b = val;
            break;
        case(SEPS):
            source->eps = val;
            break;
        case(STHETA):
            source->E.theta = val;
            break;
        case(SINDEX):
            source->var1 = val;
            break;
        case(SFLUX):
            source->mag = val;
            break;
    }
}

double  o_get_source(struct galaxie *source, int ipx)
{
    double val = 0;

    switch (ipx)
    {
        case(SCX):
            val = source->C.x;
            break;
        case(SCY):
            val = source->C.y;
            break;
        case(SA):
            val = source->E.a;
            break;
        case(SB):
            val = source->E.b;
            break;
        case(SEPS):
            val = source->eps;
            break;            
        case(STHETA):
            val = source->E.theta;
            break;
        case(SINDEX):
            val = source->var1;
            break;
        case(SFLUX):
            val = source->mag;
            break;
    }

    return val;
}
