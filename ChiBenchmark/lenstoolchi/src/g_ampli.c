#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*      nom:        g_ampli             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 *
 * Global variables used :
 * - F, M, lens, imFrame
 */
void g_ampli_m3_boucle(double  **ampli, int **namp, 
		       const double sxmin, const double symin, 
		       const double dx,    const double dy,
		       const double dlsds, const double dl0s, 
		       const double dos, const double z, const int np,
		       double  **imsens);

void    g_ampli(int iamp, int np, double z, char *file)
{
    const extern  struct  g_frame     F;
    const extern  struct  g_mode      M;
//  const extern  struct  g_cosmo     C;
    extern  struct  g_pixel     imFrame;
    const extern  struct  pot lens[];
    extern struct point gsource_global[NGGMAX][NGGMAX];
   
    register int    i, j, k, ii, jj;
    double  dl0s,dos,dlsds;
    struct  point   pi, ps;
    struct  point   pImage[NIMAX]; // list of image positions for a given source position
    struct  ellipse amp;
    struct  matrix  MA;
    double  kappa, ga1, ga2, gam, gp;
    double  **ampli;
    double  xw, yw;     // WCS coordinates in degrees for -3 mode
    int **namp;
    int ni;
    double  dx, dy, sxmin, sxmax, symin, symax;

    // initialise variables
    amp.a = amp.b = 0.;

    if (iamp == 1)
    {
        NPRINTF(stderr, "COMP: Amp in the Image Plane for z_s=%.3lf =>%s\n", z, file);
    }
    else if (iamp == 2)
    {
        NPRINTF(stderr, "COMP: abs(Amp) in the Image Plane for z_s=%.3lf =>%s\n",
                z, file);
    }
    else if (iamp == 3)
    {
        NPRINTF(stderr, "COMP:-2.5log(abs(Amp)) in the Image Plane for z_s=%.3lf =>%s\n",
                z, file);
    }
    else if (iamp == 5)
    {
        NPRINTF(stderr, "COMP:kappa in the Image Plane (not normalized by dlsds) for z_s=%.3lf =>%s\n",
                z, file);
    }
    else if (iamp == 6)
    {
        NPRINTF(stderr, "COMP:gamma in the Image Plane (not normalized by dlsds) for z_s=%.3lf =>%s\n",
                z, file);
    }
    else if (iamp == -1)
    {
        NPRINTF(stderr, "COMP:-2.5log((abs(Amp)) in the Source Plane for z_s=%.3lf =>%s\n", z, file);
    }
    else if (iamp == -3)
    {
        NPRINTF(stderr, "COMP: Correct -2.5log((abs(Amp)) in the Source Plane for z_s=%.3f =>%s using the detection map %s\n", z, file, imFrame.pixfile);
    }

    else
    {
        NPRINTF(stderr, "COMP: 1/Amp in the Source Plane %d\n", iamp);
    }

    dl0s = distcosmo2(lens[0].z, z);
    dos = distcosmo1(z);
    dlsds = dl0s / dos;

    // warning message
    if ( iamp == 5 || iamp == 6 )
    {
        extern struct g_grille G;
        double oldz = lens[0].z;
        long int i = 0;
        while ( oldz == lens[i].z && i < G.nlens ) i++;
        
        if ( i < G.nlens )
            fprintf(stderr, "WARN: case ampli %d not valid for lenses at different redshifts\n", iamp);
    }

    ampli = (double **) alloc_square_double(np, np);
    namp = (int **) alloc_square_int(np, np);

    /* Make sure we have empty arrays */
    for (j = 0; j < np; j++)
        for (i = 0; i < np; i++)
        {
            ampli[i][j] = 0.;
            namp[i][j] = 0;
        }

    if (iamp > 0)
    {
        for (j = 0; j < np; j++)
        {
            pi.y = j * (F.ymax - F.ymin) / (np - 1) + F.ymin;
            for (i = 0; i < np; i++)
            {
                pi.x = i * (F.xmax - F.xmin) / (np - 1) + F.xmin;
                amp = e_unmag(&pi, dl0s, dos, z);
                /*amplification*/
                if (iamp == 1)
                    ampli[j][i] = 1. / (amp.a * amp.b);
                /*absolute value of amplification*/
                else if (iamp == 2)
                    ampli[j][i] = 1. / fabs(amp.a * amp.b);
                /*amplification in magnitudes*/
                else if (iamp == 3)
                    ampli[j][i] = -2.5 * log10(fabs(amp.a * amp.b));
                /**/
                else if (iamp == 4)
                {
                    MA = e_grad2(&pi, dl0s, z);
                    MA.a /= dos;
                    MA.b /= dos;
                    MA.c /= dos;

                    kappa = (MA.a + MA.c) / 2.;
                    ga1 = (MA.a - MA.c) / 2.;
                    ga2 = MA.b;
                    gam = sqrt(ga1 * ga1 + ga2 * ga2); /*gamma*/
                    gp = gam / (1 - kappa);
                    ampli[j][i] = (1 - kappa) * (1 + gp * gp) / (1 - gp * gp);
                }
                else if (iamp == 5 || iamp == 6)
                {
                    MA = e_grad2(&pi, dl0s, z);
                    /*
                                            MA.a*=dlsds;
                                            MA.b*=dlsds;
                                            MA.c*=dlsds;
                                            MA.d*=dlsds;
                    */
                    MA.a /= dl0s;
                    MA.b /= dl0s;
                    MA.c /= dl0s;

                    kappa = (MA.a + MA.c) / 2.;
                    ga1 = (MA.a - MA.c) / 2.;
                    ga2 = MA.b;
                    gam = sqrt(ga1 * ga1 + ga2 * ga2);
                    if (iamp == 5)
                        ampli[j][i] = kappa;
                    else if (iamp == 6)
                        ampli[j][i] = gam;
                }
                /*amplification^-1*/
                else
                    ampli[j][i] = (amp.a * amp.b);
            };
        };

        if (M.iref > 0)
        {
            wrf_fits_abs(file, ampli, np, np, F.xmin, F.xmax, F.ymin, F.ymax, M.ref_ra, M.ref_dec);
        }
        else
        {
            wrf_fits(file, ampli, np, np, F.xmin, F.xmax, F.ymin, F.ymax);
        }

    }
    /*amplification in the source plane*/
    else if (iamp < 0)
    {
        /*define a smaller window in the Source plane*/
        dx = (F.xmax - F.xmin) / 6.;
        dy = (F.ymax - F.ymin) / 6.;
        sxmin = F.xmin + dx;
        sxmax = F.xmax - dx;
        symin = F.ymin + dy;
        symax = F.ymax - dy;
        dx = (sxmax - sxmin) / (np - 1);
        dy = (symax - symin) / (np - 1);

        if (iamp == -1)
        {
            for (j = 0; j < (np*1.8); j++)
            {
                pi.y = j * (F.ymax - F.ymin) / (np * 1.8 - 1) + F.ymin;
                for (i = 0; i < (np*1.8); i++)
                {
                    pi.x = i * (F.xmax - F.xmin) / (np * 1.8 - 1) + F.xmin;
                    amp = e_unmag(&pi, dl0s, dos, z);
                    e_dpl(&pi, dlsds, &ps);
                    ii = (int) (0.5 + (ps.x - sxmin) / dx);
                    jj = (int) (0.5 + (ps.y - symin) / dy);

                    if ((ii >= 0) && (ii < np) && (jj >= 0) && (jj < np))
                    {
                        ampli[jj][ii] += 1. / (fabs(amp.a * amp.b));
                        namp[jj][ii]++;
                    }
                }
            }

            for (ii = 0; ii < np; ii++)
                for (jj = 0; jj < np; jj++)
                    if (namp[jj][ii] > 0)
                    {
                        ampli[jj][ii] /= namp[jj][ii];
                        ampli[jj][ii] = 2.5 * log10(ampli[jj][ii]);
                    }
        }
        /*amplification total of all arclets of a familly in the source plane
         * not finished*/
        else if (iamp == -2)
        {
            e_unlensgrid(gsource_global, dlsds);

            for (j = 0; j < np; j++)
            {
                ps.y = j * dy + symin;
                for (i = 0; i < np; i++)
                {
                    ps.x = i * dx + sxmin;
                    ni = e_lens_P(ps, pImage, dlsds);
                    for (k = 0; k < ni; k++)
                    {
                        amp = e_unmag(&pImage[k], dl0s, dos, z);
                        ampli[j][i] += 1. / (fabs(amp.a * amp.b));
                    }
                }
            }

            for (ii = 0; ii < np; ii++)
                for (jj = 0; jj < np; jj++)
                    if (namp[jj][ii] > 0)
                        ampli[jj][ii] = 2.5 * log10(ampli[jj][ii]);
        }
        else if ( iamp == -3 )
        {
	   extern struct point gsource_global[NGGMAX][NGGMAX];   
	   extern  struct  g_pixel     imFrame;
	   const extern  struct  g_frame     F;
	   
	   double  **imsens;   // sensitivity map for the -3 mode in image plane
	   
	   //read image
	   imsens = (double **) readimage(&imFrame);
	   e_unlensgrid(gsource_global, dlsds);
	   dx = (F.xmax - F.xmin) / 6.;
	   dy = (F.ymax - F.ymin) / 6.;
	   sxmin = F.xmin + dx;    /*define a smaller window in the Source plane*/
	   sxmax = F.xmax - dx;
	   symin = F.ymin + dy;
	   symax = F.ymax - dy;
	   dx = (sxmax - sxmin) / (np - 1);
	   dy = (symax - symin) / (np - 1);
   
	   fprintf(stderr, "SIZE X:%d(pix) Y:%d(pix) DX:%.2lf(arcsec/pix) DY:%.2lf(arcsec/pix)\n",
		   np, np, dx, dy);
	   
	   g_ampli_m3_boucle(ampli, namp, 
			     sxmin, symin, 
			     dx,   dy,
			     dlsds,  dl0s, 
			     dos, z, np, imsens);
	   
	   free_square_double(imsens, imFrame.ny);
	   fprintf(stderr, "> done \n");
        }

        // write the fits amplification map
        if (M.iref > 0)
        {
            wrf_fits_abs(file, ampli, np, np, sxmin, sxmax, symin, symax, M.ref_ra, M.ref_dec);
        }
        else
        {
            wrf_fits(file, ampli, np, np, sxmin, sxmax, symin, symax);
        }
    } // end if iamp<0
    else
    {
        NPRINTF(stderr, "WARNING: Command ampli not recognise\n");
    }

    free_square_double(ampli, np);
    free_square_int(namp, np);
}


void g_ampli_m3_boucle(double  **ampli, int **namp, 
		       const double sxmin, const double symin, 
		       const double dx,    const double dy,
		       const double dlsds, const double dl0s, 
		       const double dos, const double z, const int np,
		       double  **imsens)
{
   extern  struct  g_pixel     imFrame;
   const extern  struct  g_mode      M;
   
   //We have the following non-constant parameters:
   //ampli, namp   --- main result --- shared
   //imFrame       --- wcs2pix problem (solved by critical), the rest used as constant
   //imsens        --- used as constant
   
   //we call the following function:
   //e_lens_P  --- thread save (don't have non constant extern and static in all calls inside)
   //e_unmag   --- thread save (don't have non constant extern and static in all calls inside)
   //wcs2pix  --- unfortunately change imFrame.wcsinfo, probably use it as temporaly storage
   //             MUST be in critical block
   //
   
   int j;      
   #pragma omp parallel for schedule(dynamic,1)
   for (j = 0; j < np; j++)
     {
	struct  point   ps;
	ps.y = j * dy + symin;
	
	int i;
	for (i = 0; i < np; i++)
	  {
	     ps.x = i * dx + sxmin;
	     
	     struct  point   pImage[NIMAX]; // list of image positions for a given source position
	     int ni = e_lens_P(ps, pImage, dlsds);
	     
	     int k;
	     for (k = 0; k < ni; k++)
	       {
		  struct  ellipse amp = e_unmag(&pImage[k], dl0s, dos, z);
		  /* compute the corresponding pixel in the image plane map*/
		  
		  int imii, imij; // coordinates in pixels of the FITS file for the sensitivity map
		  if ( imFrame.wcsinfo != NULL )
		    {
		       double yw = pImage[k].y / 3600. + M.ref_dec;
		       double xw = -pImage[k].x / 3600. / cos(M.ref_dec * DTR) + M.ref_ra;
		       
		       double  xpix, ypix; // Image coordinates in pixels for -3 mode
		       int offscl;     // offset for wcs2pix
		       
		       //in principle imFrame.wcsinfo should be constant parameter, but
		       //actually it is not like this... for some reason wcs2pix change 
		       //imFrame.wcsinfo.... there is possiblity that he use some fields like 
		       //temporaly variables... 
		       // This MUST be critical!
#pragma omp critical
			 {
			    wcs2pix(imFrame.wcsinfo, xw, yw, &xpix, &ypix, &offscl);
			 }
		       
		       imii = (int) ypix;
		       imij = (int) xpix;
		    }
		  else
		    {
		       imii = (int) ((pImage[k].y - imFrame.ymin) / imFrame.pixely);
		       imij = (int) ((pImage[k].x - imFrame.xmin) / imFrame.pixelx);
		    }
		       
		       if ((imii >= 0) && (imij >= 0) && (imii < imFrame.ny) && (imij < imFrame.nx))
			 {
			    // sensitivity of 1 pixel in source plane
                            double ssens = imsens[imii][imij] * (fabs(amp.a * amp.b));
                            // look for the minimum ssens
                            if (ssens > 0)
			      {
				 if ( ampli[j][i] == 0)
				   ampli[j][i] = ssens;
				 else
				   {
				      if ( ampli[j][i] > ssens) ampli[j][i] = ssens;
				   }
				 namp[j][i]++;
			      }
			 }
	       }
	  } //end of i loop
	//amp modified inside k loop!!! we cannot print it outside i loop!
//	NPRINTF(stderr, "%d/%d %.2lf %.2lf\r",
//		j, np, (*pampli)[j][i], fabs(amp.a*amp.b) );
     } // end of j loop 
}
