#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        update_emass            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       8/2005              */
/*      place:      Caltech             */
/****************************************************************/

void    update_emass(int i)
{
    extern struct pot lens[];


    switch (lens[i].type)
    {
        case(12): // NFW
            lens[i].emass = 3.*lens[i].epot;
            break;
        case(13): // Sersic
            lens[i].emass = 3.*lens[i].epot;
            break;
        default:
            break;
    }

}
