#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>


/****************************************************************/
/*      nom:        g_shear             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    g_shear(int ishear, int np, double z, char *file)
{
    const extern  struct  g_mode          M;
    const extern  struct  g_frame     F;
    const extern  struct  pot             lens[];

    register int    i, j, ii, jj;
    double  q, conv;
    double  dl0s, dos, dlsds; 
    struct  point   ps, pi;
    struct  ellipse ampli;
    double  **shear;
    int **nsh;
    struct matrix  grad2;
    double  dx, dy, xmin, xmax, ymin, ymax;

    if (ishear == 1)
    {
        NPRINTF(stderr, "COMP:shear_map (gamma) in the Image Plane for z_s=%.3lf =>%s\n", z, file);
    }
    else if (ishear == 2)
    {
        NPRINTF(stderr, "COMP:shear_map (eps) in the Image Plane for z_s=%.3lf=>%s\n",
                z, file);
    }
    else if (ishear == 3)
    {
        NPRINTF(stderr, "COMP:shear_map (gamma1) in the Image Plane for z_s=%.3lf=>%s\n",
                z, file);
    }
    else if (ishear == 4)
    {
        NPRINTF(stderr, "COMP:shear_map (gamma2) in the Image Plane for z_s=%.3lf=>%s\n",
                z, file);
    }
    else if (ishear == -1)
    {
        NPRINTF(stderr, "COMP:shear_map (gamma) in the Source Plane for z_s=%.3lf =>%s\n",
                z, file);
    }
    else if (ishear == -2)
    {
        NPRINTF(stderr, "COMP:shear_map (eps) in the Source Plane for z_s=%.3lf =>%s\n",
                z, file);
    }

    dl0s = distcosmo2( lens[0].z, z );
    dos = distcosmo1( z );
    dlsds = dl0s / dos;

    shear = (double **) alloc_square_double(np, np);
    nsh = (int **) alloc_square_int(np, np);

    // Compute the shear at the positions in the image plane
    if (ishear > 0)
    {
        xmin = F.xmin;
        xmax = F.xmax;
        ymin = F.ymin;
        ymax = F.ymax;

        // Compute a shear map
        if ( ishear == 1 )
        {
            for (j = 0; j < np; j++)
            {
                pi.y = j * (ymax - ymin) / (np - 1) + ymin;
                for (i = 0; i < np; i++)
                {
                    pi.x = i * (xmax - xmin) / (np - 1) + xmin;
                    ampli = e_unmag(&pi, dl0s, dos, z);
                    shear[j][i] = (ampli.a - ampli.b) / 2.;
                }
            }
        }
        // Compute an ellipticity map
        else if ( ishear == 2 )
        {
            for (j = 0; j < np; j++)
            {
                pi.y = j * (ymax - ymin) / (np - 1) + ymin;
                for (i = 0; i < np; i++)
                {
                    pi.x = i * (xmax - xmin) / (np - 1) + xmin;
                    ampli = e_unmag(&pi, dl0s, dos, z);
                    q = ampli.a / ampli.b;
                    shear[j][i] = fabs((q * q - 1.) / (q * q + 1.));
                }
            }
        }
        // Compute the gamma1 shear component
        else if ( ishear == 3 )
        {
            //conv = vol * vol / 4.*distcosmo1(lens[0].z) * dlsds;
            //conv = dlsds / 2.;
	    conv = 1./ dos / 2.;
            for (j = 0; j < np; j++)
            {
                pi.y = j * (ymax - ymin) / (np - 1) + ymin;
                for (i = 0; i < np; i++)
                {
                    pi.x = i * (xmax - xmin) / (np - 1) + xmin;
                    grad2 = e_grad2(&pi, dl0s, z);
                    shear[j][i] = conv * (grad2.a - grad2.c);
                }
            }
        }
        // Compute the gamma2 shear component
        else if ( ishear == 4 )
        {
            //conv = vol * vol / 4.*distcosmo1(lens[0].z) * dlsds;
            //conv = dlsds;
	    conv = 1. / dos;
            for (j = 0; j < np; j++)
            {
                pi.y = j * (ymax - ymin) / (np - 1) + ymin;
                for (i = 0; i < np; i++)
                {
                    pi.x = i * (xmax - xmin) / (np - 1) + xmin;
                    grad2 = e_grad2(&pi, dl0s, z);
                    shear[j][i] = conv * grad2.b;
                }
            }
        }
    }
    else
        // Compute the shear at the positions in the source plane
        // (what will be felt by each pixel of the source plane)
    {
        dx = (F.xmax - F.xmin) / 6.;
        dy = (F.ymax - F.ymin) / 6.;
        xmin = F.xmin + dx; /*define a smaller window in the Source plane*/
        xmax = F.xmax - dx;
        ymin = F.ymin + dy;
        ymax = F.ymax - dy;
        dx = (xmax - xmin) / (np - 1);
        dy = (ymax - ymin) / (np - 1);

        for (j = 0; j < (np*1.8); j++)
        {
            pi.y = j * (F.ymax - F.ymin) / (np * 1.8 - 1) + F.ymin;
            for (i = 0; i < (np*1.8); i++)
            {
                pi.x = i * (F.xmax - F.xmin) / (np * 1.8 - 1) + F.xmin;
                ampli = e_unmag(&pi, dl0s, dos, z);
                e_dpl(&pi, dlsds, &ps);
                ii = (int) (0.5 + (ps.x - xmin) / dx);
                jj = (int) (0.5 + (ps.y - ymin) / dy);

                if ((ii >= 0) && (ii < np) && (jj >= 0) && (jj < np))
                {
                    if ( ishear == -1)
                        shear[jj][ii] += (ampli.a - ampli.b) / 2.;
                    else if ( ishear == -2)
                    {
                        q = ampli.a / ampli.b;
                        shear[jj][ii] += fabs((q * q - 1.) / (q * q + 1.));
                    }
                    nsh[jj][ii]++;
                }
            }
        }

        for (ii = 0; ii < np; ii++)
            for (jj = 0; jj < np; jj++)
                if (nsh[jj][ii] > 0)
                    shear[jj][ii] /= nsh[jj][ii];
    }

    if (M.iref > 0)
        wrf_fits_abs(file, shear, np, np, xmin, xmax, ymin, ymax, M.ref_ra, M.ref_dec);
    else
        wrf_fits(file, shear, np, np, xmin, xmax, ymin, ymax);

    free_square_double(shear, np);
    free_square_int(nsh, np);
}


