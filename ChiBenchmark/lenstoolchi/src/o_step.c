#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

double  x1min, x1max, y1min, y1max;
double  x2min, x2max, y2min, y2max;
int imapmin, imapmax;
int jmapmin, jmapmax;
int ipmin, ipmax; /*parameter index that have minimise or maximise chi2*/
int ilsmin, ilsmax;
int izmin, izmax;

/****************************************************************/
/*      nom:        o_step              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            *
 ***************************************************************
 * Return the min chi0 obtained from all the parameters
 */
double  o_step(double chi0)
{
//  extern  struct  g_mode          M;
    extern  struct  pot lens[];
    extern  struct  g_grille G;

    double   y0;
    int stop;   /*number of parameters which variation ranges */
    /* are in the error bars*/

    izmin = izmax = -1; /* arclet parameter related to z_m_limit in image keyword*/
    ipmin = ipmax = -1; /* lens parameter related to limit keyword */
    ilsmin = ilsmax = -1;   /* index of the lens */
    imapmin = imapmax = -1; /* map parameter if opt limit parameter is negative */
    jmapmin = jmapmax = -1;

    y0 = chi0;
    y1min = chi0;
    y1max = chi0;

    /* o_big_slope return the number of parametres that are not in the error
     * bars and have to be modified*/
    if (lens[0].type != 10)
    {
        stop = o_big_slope();

        if (stop != 0)
        {
            /* Check if chi2 min/max has been modified by a potential,
             * a zm_limit or a map parameter*/
            if ((ipmin != -1) || (izmin != -1) || (imapmin != -1))
            {
                chi0 = o_min_slope(y0); /*chi0 is not modified but... it's
                                        safer to keep its value in y0*/
                return(chi0);
            }
            /* If a parameter has modified the max chi2 limit*/
            else if ((ipmax != -1) || (izmax != -1) || (imapmax != -1))
            {
                chi0 = o_min_loc(y0);
                return(chi0);
            }
            /* Warning if no chi2 modification and stop>0*/
            else
            {
                fprintf(stderr, "WARNING: o_step not correct\n");
                return(chi0);
            }
        }
        /* No chi2 modification */
        else
            return(-chi0);
    }
    /* lens[0].type is 10*/
    else
    {
        stop = o_slope_sp(&chi0);
        if (stop == 0)
            G.exc *= .5;
        if (G.exc < G.excmin)
            return(-chi0);
        else
            return(chi0);
    };
} /*end*/
