#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_run               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void  o_run()
{
    extern  struct  g_mode  M;
    extern  double  chip, chis, chil, chia; //,chix,chiy;
    extern  int     optim_z;
    extern  int     nwarn;

    double  chi0;
    int     stop = 0;   /* Trigger the end of the iterative optimization process */
    int     i = 0;      /* Number of iterations */
    int     ntot = 0;

    /* Initialize variables for the optimisation */
    chi0 = o_prep();
    NPRINTF(stderr, "INFO: initial chi %lf\n", chi0);

    /* boucle d'optim, jusqu'a atteindre la prec desire sur les parametres */
    optim_z = 0;

    do
    {
        /* o_step() return the min chi0 obtained from all the parameters*/
        chi0 = o_step(chi0);
        if (chi0 < 0.)
        {
            chi0 = -chi0;
            stop = 1;
        };

        i++;

        ntot += nwarn;  /* No use */

        NPRINTF(stderr, "%d/%d(%d) %.3lf p:%.3lf s:%.3lf l:%.3lf a:%.3lf\n",
                i, M.itmax, stop, chi0, chip, chis, chil, chia);
        NPRINTF(stderr, "chia=%.3lf  chip=%.3lf\n", chia, chip);
        NPRINTF(stderr, "n WARNING= %d n_total=%d\n\n", nwarn, ntot);
    }
    while ((chi0 > M.minchi0) && (stop != 1) && (i < M.itmax));

    // Free the structures
}
