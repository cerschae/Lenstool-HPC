#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        ctimp               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

double iter(double phi0, double coeur, double ct)
{
    double  newct;

    newct = sqrt(phi0 / 2.*log(1. + ct * ct / coeur / coeur));
    if (fabs(newct - ct) < 1e-4)
        return(newct);
    else
        ct = iter(phi0, coeur, newct);

    return 0.; //just to avoid the warning
}

