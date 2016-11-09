#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*      nom:        e-giant             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Contour the giant arc in the source plane and for each point try to
 * get a familly of arclets.
 *
 * Parameters :
 * - ng : current number of giant arcs
 * - giants : source of the giant arc
 * - image : array of arclets for the source giants
 * */
void    e_giant(int *ng, struct galaxie giants, struct galaxie *image )
{
    const extern struct g_mode    M;
    const extern struct g_large   L;
    const extern struct g_frame   F;
    const extern struct g_source  S;
    //const extern    struct  chaine  *Tsol;
    const extern struct g_observ  O;
    const extern struct point gsource_global[NGGMAX][NGGMAX];
   
    struct bitriplet Tsol[NIMAX];
    int    i, j, k;
    int  nimage, nig;
    int  nx, ny, ny_sav;
    double b, scale;
    double aa, bb, theta, cot, sit, coz, siz, J, K, NR, NC;
    double xmin, ymin, xmax, ymax, delta, dx, dy;
    double **zz;
    double f;
    double r;
    char   file[50];
    struct point P, ps, pi;
    FILE   *OUT;


    nig = (*ng);

    NPRINTF(stderr, "COMP: large distorsion %s (%.3lf,%.3lf) (%.3lf,%.3lf,%.3lf)\n",
            giants.n, giants.C.x, giants.C.y, giants.E.a, giants.E.b, giants.E.theta);

    /*L.ncourbe is the number of ellipses inside the giant arc*/
    if (L.ncourbe != 0)
    {
        NC = ((double)L.npt) / 2. / PI; /*number of points per radian to contour the arc*/
        NR = L.ncourbe;
        siz = sin(giants.E.theta);
        coz = cos(giants.E.theta);

        /*For each ellipse inside the giant arc*/
        for (j = 0; j < L.ncourbe; j++)
        {
            J = ((double)(j + 1)) / NR; /*scale a and b from the inner to the outer ellipse*/
            aa = giants.E.a * J;
            bb = giants.E.b * J;

            for (k = 0; k < L.npt; k++)
            {
                K = k;
                theta = K / NC;
                cot = cos(theta);
                sit = sin(theta);
                P.x = giants.C.x + aa * cot * coz - bb * sit * siz;
                P.y = giants.C.y + aa * cot * siz + bb * sit * coz;

                //Tsol=NULL;
                nimage = inverse(gsource_global, &P, Tsol); /*return the number of arclets for the source P*/
                if (nimage > 0)
                {
                    e_testg(nig, 0, Tsol, nimage, &P, giants.dr); /*test if its possible to precise the position of the arclets*/
                    nig++;
                };
            };
        };
    }
    /*Representation as a gaussian profile for the source*/
    else if (L.profil != 0)
    {
        OUT = fopen("giants.dat", "w");
        i = S.rand;
        siz = sin(giants.E.theta);
        coz = cos(giants.E.theta);
        aa = giants.E.a;
        bb = giants.E.b;
        for (j = 0; j < L.pt; j++)
        {
            do
                r = .5 * sqrt(-2.*log(d_random(&i)));
            while (r > 1.);

            theta = 2.*PI * d_random(&i);
            cot = cos(theta);
            sit = sin(theta);
            P.x = giants.C.x + r * (aa * cot * coz - bb * sit * siz);
            P.y = giants.C.y + r * (aa * cot * siz + bb * sit * coz);
            fprintf(OUT, "%.3lf %.3lf \n", P.x, P.y);
            //Tsol=NULL;
            nimage = inverse(gsource_global, &P, Tsol);
            if (nimage > 0)
            {
                e_testg(nig, 0, Tsol, nimage, &P, giants.dr); /*Set n of the arclet 0 of the familly nig to "1"*/
                nig++;
            };
        };
    }
    /*Draw the giant arc in a serie of fits file*/
    else if (L.iso != 0)
    {
        i = 0;
        sprintf(file, "%s%s.fits", L.iname, giants.n);

        xmin = F.xmin;
        xmax = F.xmax;
        ymin = F.ymin;
        ymax = F.ymax;

        if (L.iso == 1)
        {
            while (!strcmp(image[i].n, ""))
            {
                P = image[i].C;
                if (i == 0)
                {
                    xmax = xmin = P.x;
                    ymax = ymin = P.y;
                    i++;
                }
                else
                {
                    xmin = Min(P.x, xmin);
                    xmax = Max(P.x, xmax);
                    ymin = Min(P.y, ymin);
                    ymax = Max(P.y, ymax);
                    i++;
                };

//              Tsol=(Tsol->F);
            }
        }

        if (xmin == xmax)
        {
            xmin = xmin - giants.E.a * 3.;
            xmax = xmax + giants.E.a * 3.;
        }
        if (ymin == ymax)
        {
            ymin = ymin - giants.E.a * 3.;
            ymax = ymax + giants.E.a * 3.;
        }

        dx = xmax - xmin;
        dy = ymax - ymin;
        xmin = xmin - L.zonex * dx;
        xmax = xmax + L.zonex * dx;
        ymin = ymin - L.zoney * dy;
        ymax = ymax + L.zoney * dy;
        dx = xmax - xmin;
        dy = ymax - ymin;

        // delta : larger image side
        delta = Max(dx, dy);

        // set scale (the pixel size) so that the larger
        // image side (delta) is smaller than L.nmaxiso
        if ( delta < L.scale*(L.nmaxiso - 1) )
            scale = L.scale;
        else
            scale = delta / (L.nmaxiso - 1);

        nx = 1 + ((int)(dx / scale));
        ny = 1 + ((int)(dy / scale));

        // extend the image size if binning...
        if (O.setbin)
        {
            b = O.bin;
            nx = (int) b * nx;
            ny = (int) b * ny;
            delta = scale * (b - 1) / 2 / b;
            xmin -= delta;
            xmax += delta;
            ymin -= delta;
            ymax += delta;
            scale = scale / ((double)b); // ... and the pixel size
        }

        zz = (double **)alloc_square_double(ny, nx);
        for (j = 0; j < ny; j++)
        {
            pi.y = ymin + ((double)j + 1) * scale;
            for (k = 0; k < nx; k++)
            {
//              NPRINTF(stderr, "%d %d\n",j,k);
                pi.x = xmin + ((double)k + 1) * scale;
                e_dpl(&pi, giants.dr, &ps);
                giants.c = 'g';
                giants.I0 = 50.;
                f = d_profil(ps.x, ps.y, &giants);
                /*
                            if ((f>sqrt(O.SKY/O.gain)) || (f>giants.I0/100.))
                                            zz[k][j]=d_integrer(giants,pi.x,pi.y,scale/2.,f,0);
                            else
                */
                zz[j][k] = f;
            }
        }

        if (O.setseeing)
            d_seeing(zz, nx, ny, scale);

        ny_sav = ny;
        if (O.setbin)
        {
            d_binning(zz, &nx, &ny, O.bin);
            xmin += delta;
            xmax -= delta;
            ymin += delta;
            ymax -= delta;
        }

        if (O.bruit)
            d_bruiter(zz, nx, ny);

        if ( M.iref > 0 )
            wrf_fits_abs(file, zz, nx, ny, xmin, xmax, ymin, ymax, M.ref_ra, M.ref_dec);
        else
            wrf_fits(file, zz, nx, ny, xmin, xmax, ymin, ymax);

        free_square_double(zz, ny_sav);
    } /*end of the iso representation*/

    *ng = nig;
}
