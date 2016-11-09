#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "fonction.h"
#include "constant.h"
#include "dimension.h"
#include "structure.h"
#include<lt.h>

/****************************************************************/
/*                 Program  : grille                */
/*                 Version  : 1 mai 1992            */
/*                 Location : Obs. Toulouse         */
/*                 Auteur   : jean-paul             */
/****************************************************************/

void r_potfile(FILE *IN, FILE *OUT, struct g_pot *pot)
{
    extern  struct  g_mode          M;
    char    second[20], third[120];

    strcpy(pot->potfile, "");

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(120)

        if (!strcmp(second, "filein"))
        {
            sscanf(third, "%d%s", &pot->ftype, pot->potfile);
            fprintf(OUT, "\t%s\t%d %s\n", second, pot->ftype, pot->potfile);
        }
        else if (!strcmp(second, "type"))
        {
            sscanf(third, "%d", &pot->type);
            fprintf(OUT, "\t%s\t%d\n", second, pot->type);
        }
        else if (!strcmp(second, "zlens") || !strcmp(second, "z_lens") )
        {
            sscanf(third, "%lf", &pot->zlens);
            fprintf(OUT, "\t%s\t%lf\n", second, pot->zlens);
        }
        else if (!strcmp(second, "mag0") || !strcmp(second, "r200"))
        {
            sscanf(third, "%lf", &pot->mag0);
            fprintf(OUT, "\t%s\t%lf\n", second, pot->mag0);
        }
        else if (!strcmp(second, "select"))
        {
            sscanf(third, "%d", &pot->select);
            fprintf(OUT, "\t%s\t%d\n", second, pot->select);
        }
        else if (!strcmp(second, "core"))
        {
            sscanf(third, "%lf", &pot->core);
            fprintf(OUT, "\t%s\t%lf\n", second, pot->core);
        }
        else if (!strcmp(second, "corekpc"))
        {
            sscanf(third, "%lf", &pot->corekpc);
            fprintf(OUT, "\t%s\t%lf\n", second, pot->corekpc);
        }
        else if (!strcmp(second, "cut"))
        {
            sscanf(third, "%d%lf%lf", &pot->ircut, &pot->cut1, &pot->cut2);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, pot->ircut, pot->cut1, pot->cut2);
        }
        else if (!strcmp(second, "cutkpc"))
        {
            sscanf(third, "%d%lf%lf", &pot->ircut, &pot->cutkpc1, &pot->cutkpc2);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, pot->ircut, pot->cutkpc1, pot->cutkpc2);
        }
        else if (!strcmp(second, "slope") || !strcmp(second, "m200slope"))
        {
            sscanf(third, "%d%lf%lf", &pot->islope, &pot->slope1, &pot->slope2);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, pot->islope, pot->slope1, pot->slope2);
        }
        else if (!strcmp(second, "sigma"))
        {
            sscanf(third, "%d%lf%lf", &pot->isigma, &pot->sigma1, &pot->sigma2);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, pot->isigma, pot->sigma1, pot->sigma2);
        }
        else if (!strcmp(second, "vdslope") || !strcmp(second, "c200slope"))
        {
            sscanf(third, "%d%lf%lf", &pot->ivdslope, &pot->vdslope1, &pot->vdslope2);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, pot->ivdslope,
                    pot->vdslope1, pot->vdslope2);
        }
        else if (!strncmp(second, "vdscat", 6))
        {
            sscanf(third, "%d%lf%lf", &pot->ivdscat, &pot->vdscat1, &pot->vdscat2);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, pot->ivdscat,
                    pot->vdscat1, pot->vdscat2);
        }
        else if (!strncmp(second, "rcutscat", 8))
        {
            sscanf(third, "%d%lf%lf", &pot->ircutscat, &pot->rcutscat1, &pot->rcutscat2);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, pot->ircutscat,
                    pot->rcutscat1, pot->rcutscat2);
        }
        else if (!strcmp(second, "a") || !strcmp(second, "m200"))
        {
            sscanf(third, "%d%lf%lf", &pot->ia, &pot->a1, &pot->a2);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, pot->ia,pot->a1,pot->a2);
        }
        else if (!strcmp(second, "b") || !strcmp(second, "c200"))
        {
            sscanf(third, "%d%lf%lf", &pot->ib, &pot->b1, &pot->b2);
            fprintf(OUT, "\t%s\t%d %lf %lf\n", second, pot->ib,pot->b1,pot->b2);
        }


        // Read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);

    // In the case, potfile is just used to define mag0
    if ( !strcmp(pot->potfile, "") )
        return;

}
