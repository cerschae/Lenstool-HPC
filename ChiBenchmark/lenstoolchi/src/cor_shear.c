#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        cor_shear           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    cor_shear()
{
    extern  struct  g_mode      M;
//  extern  struct  g_source    S;

    struct  galaxie     shear[NAMAX], coshear[NAMAX];
    register long int    i;
    long int nsh = 0;

    NPRINTF(stderr, "COR: shear\r");

    if (M.icorshear == 1)
        f_shape(&nsh, shear, M.corshfile, 1);
    else
        exit(-1);

    pro_arclet(nsh, shear);
    o_mag(nsh, shear);

    for (i = 0; i < nsh; i++)
    {
        coshear[i] = shear[i];
        coshear[i].E.theta = shear[i].E.theta - shear[i].thp;
        coshear[i].E.a = fabs(shear[i].tau * shear[i].tp * sin(2.*coshear[i].E.theta));
        coshear[i].E.b = coshear[i].E.a / 2.;
    };


    ecrire_r(0, nsh, coshear, "coshear.dat", 1);

}
