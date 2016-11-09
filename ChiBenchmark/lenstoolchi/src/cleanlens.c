#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"

/****************************************************************/
/*      nom:        cleanlens           */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/03/94            */
/*      place:      ESO LaSilla         */
/* 
 * Not finished.
 * From the moment, sample the image plane in (x,y), comptue the 
 * source position, and get the counter image pixels in pi0->pi4
 ****************************************************************/

void    cleanlens(double zimage)
{
    const extern struct   g_mode M;
    extern  struct  g_pixel imFrame; //PSF,
    const extern  struct  pot lens[];
    extern struct point gsource_global[NGGMAX][NGGMAX];
    //extern    struct  chaine  *Tsol;

    struct bitriplet Tsol[NIMAX];
    register    int i, j;
    int nimage;
    double  err2, dlsds;
    double  **image; //,**psf;
    struct  point   **ps, **pi0, **pi1, **pi2, **pi3, **pi4;
    struct  point   imap[NIMAX];


    NPRINTF(stderr, "DO: cleanlens\n");

    /* calcul du dlsds */
    dlsds = dratio(lens[0].z, zimage);

    /* reading images */
    NPRINTF(stderr, "READ: image frame\n");
    image = (double **) readimage(&imFrame);

    /* NPRINTF(stderr,"READ: psf frame\n");
       psf=(double **) readimage(&PSF);
    */

    /* Calcul du tableau des points sources */

    grid();
    e_unlensgrid(gsource_global, dlsds);

    NPRINTF(stderr, "CREATE: source point grid (%d,%d)\n", imFrame.nx, imFrame.ny);

    ps = (struct point **)al_sq_point(imFrame.nx, imFrame.ny);
    pi0 = (struct point **)al_sq_point(imFrame.nx, imFrame.ny);
    pi1 = (struct point **)al_sq_point(imFrame.nx, imFrame.ny);
    pi2 = (struct point **)al_sq_point(imFrame.nx, imFrame.ny);
    pi3 = (struct point **)al_sq_point(imFrame.nx, imFrame.ny);
    pi4 = (struct point **)al_sq_point(imFrame.nx, imFrame.ny);


    err2 = (imFrame.xmax - imFrame.xmin) /
           (imFrame.nx - 1) * (imFrame.ymax - imFrame.ymin) /
           (imFrame.ny - 1) / 9.;

    for (i = 0; i < imFrame.nx; i++)
    {
        NPRINTF(stderr, "%d\r", i);

        for (j = 0; j < imFrame.ny; j++)
        {
            pi0[i][j].x = i * (imFrame.xmax - imFrame.xmin) / (imFrame.nx - 1) + imFrame.xmin;
            pi0[i][j].y = j * (imFrame.ymax - imFrame.ymin) / (imFrame.ny - 1) + imFrame.ymin;

            e_dpl(&pi0[i][j], dlsds, &ps[i][j]);

            //Tsol=NULL;
	    nimage = inverse((const struct point (*)[NGGMAX]) gsource_global, &ps[i][j], Tsol);
            if (nimage > 1)
            {
                nimage = e_test_P(Tsol, nimage, &ps[i][j], imap, dlsds, err2);
                pi0[i][j] = imap[0];
                pi1[i][j] = imap[1];
                pi2[i][j] = imap[2];
                pi3[i][j] = imap[3];
                pi4[i][j] = imap[4];
            };
        };
    };


    fr_sq_point(ps, imFrame.nx);
    fr_sq_point(pi0, imFrame.nx);
    fr_sq_point(pi1, imFrame.nx);
    fr_sq_point(pi2, imFrame.nx);
    fr_sq_point(pi3, imFrame.nx);
    fr_sq_point(pi4, imFrame.nx);
    free_square_double(image, imFrame.ny);

}
