#include<stdio.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"

/****************************************************************/
/*      nom:        g_prop              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    zonemult()
{
    extern struct   g_mode M;
    extern  struct  g_cline     CL;
    extern  struct  g_frame     F;
    extern  struct  pot lens[];
    extern struct biline radial[], tangent[];
    extern int nrline, ntline;//,flagr,flagt;

    register int    i, j, k;
    int iline, nline[NLMAX], ncount, flag;
    double  dlsds;
    struct  point   ps, pi, line[NIMAX][NPOINT];
    //int   size[4];
    double  xmin, ymin, xmax, ymax;
    int **nimage;

    /* Definition de la fenetre de calcul */

    if (CL.dmax != 0.)
    {
        xmin = -CL.dmax;
        xmax = CL.dmax;
        ymin = -CL.dmax;
        ymax = CL.dmax;
    }
    else
    {
        xmin = F.xmin;
        xmax = F.xmax;
        ymin = F.ymin;
        ymax = F.ymax;
    }


    /*    verification du nombre de plan en z pour le calcul des zones */

    if (CL.nplan != 1)
    {
        NPRINTF(stderr, "WARNING: Image zone not computed. Too many critical lines\n");
        return;
    }

    if ( strcmp(CL.algorithm, "SNAKE") )
    {
        NPRINTF(stderr, "WARNING: image zone not computed. Critical line algorithm must be SNAKE\n");
        return;
    }

    /* calcul des zones images pour zs=CL.cz[0] */
    NPRINTF(stderr, "COMP: multiple images area in the Image plane for sources at z=%.3lf\n", CL.cz[0]);

    dlsds = dratio(lens[0].z, CL.cz[0]);

    iline = 0;
    if (ntline > 2)
    {
        ncount = 0;
        flag = tangent[0].i;
        line[iline][ncount++] = tangent[0].S;
        for (k = 1; k < ntline; k++)
        {
            if (tangent[k].i == flag)
                line[iline][ncount++] = tangent[k].S;
            else
            {
                nline[iline++] = ncount;
                ncount = 0;
                flag = tangent[k].i;
                line[iline][ncount++] = tangent[k].S;
            }
        }

        if (iline != 2)
            nline[iline++] = ncount;
    }

    if (nrline > 2)
    {
        ncount = 0;
        flag = radial[0].i;
        line[iline][ncount++] = radial[0].S;
        for (k = 1; k < nrline; k++)
        {
            if (radial[k].i == flag)
                line[iline][ncount++] = radial[k].S;
            else
            {
                nline[iline++] = ncount;
                ncount = 0;
                flag = radial[k].i;
                line[iline][ncount++] = radial[k].S;
            }
        }
        nline[iline++] = ncount;
    }

#ifdef DEBUG
    NPRINTF(stderr, "DBG: Before alloc in zone_mult\n");
    NPRINTF(stderr, "DBG: CL.npzone %d\n", CL.npzone);
#endif

    nimage = (int **) alloc_square_int(CL.npzone, CL.npzone);
#ifdef DEBUG
    NPRINTF(stderr, "DBG: After alloc in zone_mult\n");
#endif

    for (j = 0; j < CL.npzone; j++)
    {
        pi.y = j * (ymax - ymin) / (CL.npzone - 1) + ymin;

        for (i = 0; i < CL.npzone; i++)
        {
            pi.x = i * (xmax - xmin) / (CL.npzone - 1) + xmin;
            e_dpl(&pi, dlsds, &ps);

            nimage[j][i] = 1;

            for (k = 0; k < iline; k++)
                nimage[j][i] += 2 * inconvexe(ps, nline[k], line[k]);

        };
    };

    if (M.iref >0) wri_fits_abs(CL.zonefile, nimage, CL.npzone, CL.npzone, F.xmin, F.xmax, F.ymin, F.ymax,M.ref_ra,M.ref_dec);
    else wri_fits(CL.zonefile, nimage, CL.npzone, CL.npzone, F.xmin, F.xmax, F.ymin, F.ymax);

    /*
    size[0]=CL.npzone;
    size[1]=CL.npzone;
    wr_ipx(CL.zonefile,nimage,2,size,"int","bin","real",
        "zone_mult_of_image_plane",xmin,xmax,ymin,ymax);
    */

    free_square_int(nimage, CL.npzone);

}
