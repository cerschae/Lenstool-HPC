#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        o_print             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            *
 *
 * Called by the o_run*() functions to write the potential parameters
 * in the OUT file.
 *
 * Parameters :
 * - OUT : FILE pointer to the output file
 * - chi0 : the chi2 value obtained for the current potential parameters
 */
void    o_print(FILE *OUT, double chi0)
{
    extern  struct  g_grille    G;
    extern  struct  g_cosmo     C;
    extern  struct  pot lens[];
    register int    i;

//    if (I.forme == 2)
//        fprintf(OUT, "chi2 S: %.2lf I: %.2lf\n", chi0, o_chi_pos());
//    else
        fprintf(OUT, "chi2 S:%lf\n", chi0);

    for (i = 0; i < G.no_lens; i++)
        if (lens[i].type != 10)
        {
            if (lens[i].type == 1 || lens[i].type == -1 || lens[i].type == -2)
            {
                lens[i].sigma = sqrt(lens[i].b0 / 4. / pia_c2);
                fprintf(OUT, "%.3lf %.3lf  ", lens[i].C.x, lens[i].C.y);
                fprintf(OUT, "%.4lf %.3lf  ", lens[i].emass, RTD*lens[i].theta);
                fprintf(OUT, "%.2lf  ", lens[i].sigma);
            }
            else if (lens[i].type == 7)
            {
                lens[i].masse = lens[i].b0 / (4.*RTA * GM_c2) * (D0 / C.h * distcosmo1(lens[i].z));
                fprintf(OUT, "%.3lf %.3lf  ", lens[i].C.x, lens[i].C.y);
                fprintf(OUT, "%.3lf ", lens[i].masse);
            }
            else if (lens[i].type == 9)
            {
                lens[i].pmass = lens[i].b0 * cH2piG * C.h / distcosmo1(lens[i].z);
                fprintf(OUT, "%.3lf ", lens[i].pmass);
            }
            else if (lens[i].type == 5)
            {
                lens[i].sigma = vol * sqrt(lens[i].b0 / 18. / RTA);
                fprintf(OUT, "%.3lf %.3lf  ", lens[i].C.x, lens[i].C.y);
                fprintf(OUT, "%.4lf %.3lf  ", lens[i].emass, RTD*lens[i].theta);
                fprintf(OUT, "%.3lf %.2lf  ", lens[i].rc, lens[i].sigma);
            }
            else if (lens[i].type == 8 || lens[i].type == 83 || lens[i].type == 85)
            {
                lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                fprintf(OUT, "%.3lf %.3lf ", lens[i].C.x, lens[i].C.y);
                fprintf(OUT, "%.4lf %.3lf ", lens[i].emass, RTD*lens[i].theta);
                fprintf(OUT, "%.4lf %.2lf ", lens[i].rc, lens[i].sigma);
            }
            else if (lens[i].type == 81 || lens[i].type == 82 || lens[i].type == 86)
            {
                lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                fprintf(OUT, "%.3lf %.3lf ", lens[i].C.x, lens[i].C.y);
                fprintf(OUT, "%.4lf %.3lf ", lens[i].emass, RTD*lens[i].theta);
                fprintf(OUT, "%.4lf %.2lf ", lens[i].rc, lens[i].sigma);
                fprintf(OUT, "%.3lf ", lens[i].rcut);
            }
            else if (lens[i].type == 12)
            {
                lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                fprintf(OUT, "%.3lf %.3lf ", lens[i].C.x, lens[i].C.y);
                fprintf(OUT, "%.4lf %.3lf ", lens[i].emass, RTD*lens[i].theta);
                fprintf(OUT, "%.4lf %.2lf ", lens[i].rc, lens[i].sigma);
                fprintf(OUT, "%.3lf ", lens[i].alpha);
            }
            else if ( lens[i].type == 84 || lens[i].type == 87 || lens[i].type == 88 )
            {
                lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                fprintf(OUT, "%.3lf %.3lf ", lens[i].C.x, lens[i].C.y);
                fprintf(OUT, "%.4lf %.3lf ", lens[i].emass, RTD*lens[i].theta);
                fprintf(OUT, "%.3lf %.2lf ", lens[i].rc, lens[i].sigma);
                fprintf(OUT, "%.4lf %.3lf", lens[i].rcut, lens[i].alpha);
            }
            else if ( lens[i].type == 89)
            {
                lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                fprintf(OUT, "%.3lf %.3lf ", lens[i].C.x, lens[i].C.y);
                fprintf(OUT, "%.4lf %.3lf ", lens[i].emass, RTD*lens[i].theta);
                fprintf(OUT, "%.4lf %.2lf ", lens[i].rc, lens[i].sigma);
                fprintf(OUT, "%.3lf %.3lf", lens[i].beta, lens[i].alpha);
            }
            else
            {
                lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);

                fprintf(OUT, "%.3lf %.3lf ", lens[i].C.x, lens[i].C.y);
                fprintf(OUT, "%.4lf %.3lf ", lens[i].emass, RTD*lens[i].theta);
                fprintf(OUT, "%.4lf %.2lf ", lens[i].rc, lens[i].sigma);
                if ((lens[i].type == 3) || (lens[i].type == 6))
                    fprintf(OUT, "%.4lf %.4lf\t", lens[i].alpha, lens[i].beta);

                if (lens[i].type > 20)
                {
                    updatecut(i);   //update lens parameters for type 41
                    fprintf(OUT, "%.3lf ", lens[i].rcut);
                }

            }
            fprintf(OUT, "\n");
        }; // end of loop over the lenses with type != 10
    fprintf(OUT, "\n");
}
