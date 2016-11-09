#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include "fonction.h"
#include "constant.h"
#include "dimension.h"
#include "structure.h"

/****************************************************************/
/*      nom:        g_shearf            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/*   Modified :                         */
/*      EJ (02/09/05)--Print absolute coordinates       */
/****************************************************************/

void    g_shearf(int ishearf, double z, char *file, int nshearf)
{
    const extern struct g_mode   M;
    const extern struct g_frame  F;
    const extern struct pot      lens[];

    struct point     pi;
    struct pointell  **shear;
    struct pointell  *s;

    register int i, j;
    double e, q; 
    double dl0s, dos, dlsds;
    double Cx, Cy;
    FILE *OUT;

    // Create the shear array
    shear = (struct pointell **)(malloc((unsigned) nshearf * sizeof(struct pointell *)));
    for ( i = 0; i < nshearf; i++ )
        shear[i] = (struct pointell *)(malloc((unsigned) nshearf * sizeof(struct pointell)));


    if ( ishearf == 1 )
    {
        NPRINTF(stderr, "COMP: shear_field (eps) in the Image Plane for z_s=%.3lf=>%s\n",
                z, file);
    }
    else if ( ishearf == 2 )
    {
        NPRINTF(stderr, "COMP: shear_field (ori) in the Image Plane for z_s=%.3lf=>%s\n",
                z, file);
    }
    else if ( ishearf == 3 )
    {
        NPRINTF(stderr, "COMP: shear_field (tau) in the Image Plane for z_s=%.3lf=>%s\n",
                z, file);
    }
    else if ( ishearf == 4 )
    {
        NPRINTF(stderr, "COMP: shear_field (ellipse) in the Image Plane for z_s=%.3lf=>%s\n",
                z, file);
    }

    dl0s = distcosmo2( lens[0].z, z);
    dos = distcosmo1( z );
    dlsds = dl0s / dos;


    for (j = 0; j < nshearf; j++)
    {
        pi.y = j * (F.ymax - F.ymin) / (nshearf - 1) + F.ymin;
        for (i = 0; i < nshearf; i++)
        {
            pi.x = i * (F.xmax - F.xmin) / (nshearf - 1) + F.xmin;
            shear[i][j].E = e_unmag(&pi, dl0s, dos, z);
            shear[i][j].C = pi;
        }
    }

    // Write file on disk
    OUT = fopen(file, "w");

    if ( M.iref == 3 )
        fprintf(OUT, "#REFERENCE 3 %lf %lf\n", M.ref_ra, M.ref_dec);
    else
        fprintf(OUT, "#REFERENCE 0\n");

    for (j = 0; j < nshearf; j++)
        for (i = 0; i < nshearf; i++)
        {
            s = &shear[i][j];
            q = s->E.a / s->E.b;
            e = (q * q - 1.) / (q * q + 1.);
            Cy = s->C.y;
            Cx = s->C.x;

            if ( ishearf == 1 )
            {
                if (e > 0.)
                    fprintf(OUT, "%d %.6lf %.6lf %.4lf %.4lf %7.2lf %.1lf %.1lf\n", i,
                            Cx, Cy, e, 0., RTD*s->E.theta, 0., 0.);
                else
                    fprintf(OUT, "%d %.6lf %.6lf %.4lf %.4lf %7.2lf %.1lf %.1lf\n", i,
                            Cx, Cy, 0., -e, RTD*s->E.theta, 0., 0.);
            }
            else if ( ishearf == 2 )
            {
                if (e > 0.)
                    fprintf(OUT, "%d %.6lf %.6lf %.4lf %.4lf %7.2lf %.1lf %.1lf\n", i,
                            Cx, Cy, 1., 0., RTD*s->E.theta, 0., 0.);
                else
                    fprintf(OUT, "%d %.6lf %.6lf %.4lf %.4lf %7.2lf %.1lf %.1lf\n", i,
                            Cx, Cy, 0., 1., RTD*s->E.theta, 0., 0.);
            }
            else if ( ishearf == 3 )
            {
                e = (q * q - 1.) / 2 / q;
                if (e > 0.)
                    fprintf(OUT, "%d %.6lf %.6lf %.4lf %.4lf %7.2lf %.1lf %.1lf\n", i,
                            Cx, Cy, e, 0., RTD*s->E.theta, 0., 0.);
                else
                    fprintf(OUT, "%d %.6lf %.6lf %.4lf %.4lf %7.2lf %.1lf %.1lf\n", i,
                            Cx, Cy, 0., -e, RTD*s->E.theta, 0., 0.);
            }
            else if ( ishearf == 4 )
            {
                // +90 because the shear angle is radial and we want ellipses
                // tangeantially oriented
                fprintf(OUT, "%d %.6lf %.6lf %.5lf %.5lf %7.2lf %.1lf %.1lf\n", i,
                        Cx, Cy, s->E.a, s->E.b, RTD*s->E.theta + 90, 0., 0.);
            }
        };

    fclose(OUT);
}
