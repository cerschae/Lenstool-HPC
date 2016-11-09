#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        fz_dlsds            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/*                              */
/*  inversion de l'equation f(z)=z_dlsds    with the zero() function
 *  - z_dlsds is defined in o_print_res() function.
 ****************************************************************/

double  fz_dlsds(double z)
{
    const extern  double  z_dlsds;
    const extern  struct  pot lens[];

    return( z_dlsds - dratio(lens[0].z, z) );
}
