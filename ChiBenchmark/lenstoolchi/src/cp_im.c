#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        comp_im             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    cp_im(  double **im,
                int nx, int ny,
                double xmin, double xmax, double ymin, double ymax,
                struct galaxie *source,
                int nbs )
{
    register    int i, j, k;
    double  scalex, scaley, dx, dy;
    double   f;
    struct  point   ps, pi;

    dx = xmax - xmin;
    dy = ymax - ymin;
    scalex = dx / (nx - 1);
    scaley = dy / (ny - 1);

    for (j = 0; j < nx; j++)
    {
        pi.x = xmin + ((double)j) * scalex;
        for (k = 0; k < ny; k++)
        {
            pi.y = ymin + ((double)k) * scaley;
            e_dpl(&pi, source[0].dr, &ps);
            for (i = 0; i < nbs; i++)
            {
                f = g_profil(ps.x, ps.y, source[i]);
                im[k][j] += f;
            }
        };
    };

}
