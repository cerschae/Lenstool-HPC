#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<constant.h>
#include<fonction.h>
#include<lt.h>

/*
*
* s_sourcebox - lens-tool
* JP Kneib
* IoA Cambridge
* oct 1994
*
* define the source box to have as center
* the barycentre of the center_images file
*
* Global variables used :
* - M
* - in e_dpl() : G, lens, lens_table
*/

void    s_sourcebox(struct g_pixel *ps, const char *centerfile, double dlsds)
{
    const extern struct   g_mode M;
    FILE    *IN;
    char    line[128];
    char    n[10];
    int     np;
    double dmax = 0.;
    struct point A, B, C;
    double crval1, crval2;

    /*
    open the center images file
    */

    IN = fopen(centerfile, "r");
    if ( IN != NULL )
    {
        /*
        * read the center image file; compute the sources and the barycenter
        */
        np = 0;
        while (fscanf(IN, "%s%lf%lf", n, &A.x, &A.y) != -1)
        {
            flire(IN, line);    /* read the end of the line */
            np++;
            if (np == 1)
                e_dpl(&A, dlsds, &B);
            else
            {
                e_dpl(&A, dlsds, &C);
                B = wcenter(B, ((double) np), C, 1.);
                dmax = Max(dmax, dist(B, C));
            }
        }

        //printf("(B.x, B.y)  = (%lf, %lf)  dmax = %lf\n", B.x, B.y, dmax);
        /*
        * check the size of the box
        */

        if (dmax / 2. > ps->pixelx*ps->nx || dmax / 2. > ps->pixely*ps->ny)
        {
            NPRINTF(stderr,
                    "\nWARNING: The size of the Source Box is probably too small !!!\n");
        }

        /*
        * compute the actual dimension of the box
        * Hypothesis : positions in source plane are relative to [B.x, B.y]
        */

        ps->xmin = B.x - ps->pixelx * (ps->nx - 1) / 2.;
        ps->ymin = B.y - ps->pixely * (ps->ny - 1) / 2.;
        ps->xmax = ps->xmin + ps->pixelx * (ps->nx - 1);
        ps->ymax = ps->ymin + ps->pixely * (ps->ny - 1);

        /*NPRINTF(stderr, "INFO: Source box relative bounds (%.3lf:%.3lf %.3lf:%.3lf)\n",
                ps->xmin, ps->xmax, ps->ymin, ps->ymax); */

        /*
         * Create a mapping from wcs coordinates to pixel... new ra/dec in source plane : [B.x,B.y]
         */
        strcpy(n, "RA---TAN");
        strcpy(line, "DEC--TAN");
        crval1 = -B.x / 3600. / cos(M.ref_dec * DTR) + M.ref_ra;
        crval2 = B.y / 3600. + M.ref_dec;
        ps->wcsinfo = wcskinit(ps->nx, ps->ny, n, line, ps->nx / 2., ps->ny / 2., crval1, crval2,
                               NULL, ps->pixelx / 3600., ps->pixely / 3600., 0., 2000, 0.);

    }

    /*
    * file not found, quit.
    */
    else
    {
        fprintf(stderr, "ERROR: file %s not found\n", centerfile);
        exit(-1);
    }
    fclose(IN);

}
