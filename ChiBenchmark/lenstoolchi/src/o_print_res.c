#include<stdio.h>
#include<string.h>
#include<time.h>
#include<math.h>
#include<float.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*      nom:        o_print_res         */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/*                              */
/*   Modified :                         */
/*      EJ (01/09/05)--reference in best.par, index of the objects
 *
 ****************************************************************/

static void writePotentiel(FILE *best, long int i, int flag);
static void writeLimit(FILE *best, long int i);

void  o_print_res(double chi0, double evidence)
{
    extern struct g_mode    M;
    extern struct g_grille  G;
    extern struct g_image   I;
    extern struct g_frame   F;
    extern struct g_cline   CL;
    extern struct g_source  S;
    extern struct g_large   L;
    extern struct g_observ  O;
    extern struct g_cosmo   C;
    extern struct vfield    vf;
	extern struct	g_dyn	Dy;       //I added this  TV
    extern struct g_pot     P[NPOTFILE];
    extern struct pot       lens[NLMAX], lmin[], lmax[], prec[];
    extern struct pot       lmin_s[NLMAX], lmax_s[], prec_s[];
    extern struct galaxie   smin[NFMAX], smax[NFMAX];
    extern struct z_lim     zlim[];
    extern struct z_lim     zalim;
    extern struct z_lim     zlim_s[];
    extern struct cline     cl[];
    extern struct galaxie   multi[NFMAX][NIMAX];
    extern struct galaxie   source[NFMAX];
    extern struct galaxie   arclet[NAMAX];
    extern struct g_cosmo   clmin, clmax;
    extern struct vfield    vfmin, vfmax;
    extern int cblock[NPAMAX];
    extern int vfblock[NPAMAX];

    extern int block[][NPAMAX];
    extern int block_s[][NPAMAX];
    extern int sblock[NFMAX][NPAMAX];
    extern int nwarn;
    extern struct sigposStr sigposAs;

    extern double z_dlsds;
    extern double chip, chix, chiy, chis, chil, chi_vel,chi_mass;    //I added chi_vel and chi_mass TV
    extern double **map_p, **map_axx, **map_ayy;
    extern struct g_pixel ps, imFrame;
    int     i, j;
    FILE    *best, *besto;
    char    limages[ZMBOUND][IDSIZE];
    time_t  rawtime;

    NPRINTF(stderr, "\n*********************************************************\n");
    NPRINTF(stderr, "Chi2:%.3lf\t p:%.3lf s:%.3lf l:%.3lf\n", chi0, chip, chis, chil);
    if ( M.inverse == 3 )
        NPRINTF( stderr, "log(Evidence) : %.3lf\n", evidence);

    if (lens[0].type != 10)
    {
        for (i = 0; i < G.no_lens; i++)
        {
            NPRINTF(stderr, "%s : c (%.3lf,%.3lf) e (%.3lf,%.3lf) epot %.3lf\n",
                    lens[i].n, lens[i].C.x, lens[i].C.y, lens[i].emass, lens[i].theta*RTD, lens[i].epot);

            if ( lens[i].type == 12 )
            {
                NPRINTF(stderr, "c %.2lf rhos %.2leMsol/Mpc3 M200 %.3leMsol\n", lens[i].beta, lens[i].pmass, lens[i].masse);

                NPRINTF(stderr, "rs %.2lf\"(%.2lfkpc) s0 %.2lfkm/s alpha %.3lf\n",
                        lens[i].rc, lens[i].rckpc, lens[i].sigma, lens[i].alpha);

                if ( lens[i].rcut != DBL_MAX )
                    NPRINTF(stderr, "R200 %.2lf\"(%.2lfkpc)\n", lens[i].rcut, lens[i].rcutkpc);

            }
            else if ( lens[i].type == 13 )
            {
                NPRINTF(stderr, "re %.2lf(%.2lfkpc) sigma_e %.2le Msol/kpc^2 n %.3lf\n",
                        lens[i].rc, lens[i].rckpc, lens[i].sigma, lens[i].alpha);
            }
            else if ( lens[i].type == 16 )
            {
                NPRINTF(stderr, "rs %.2lf\"(%.2lfkpc) s0 %.2lfkm/s\n",
                        lens[i].rc, lens[i].rckpc, lens[i].sigma);
            }
            else
            {
                NPRINTF(stderr, "rc %.2lf(%.2lfkpc) s0 %.2lfkm/s alpha %.3lf\n",
                        lens[i].rc, lens[i].rckpc, lens[i].sigma, lens[i].alpha);

                if ( lens[i].rcut != DBL_MAX )
                    NPRINTF(stderr, "rcut %.2lf(%.2lfkpc)\n", lens[i].rcut, lens[i].rcutkpc);
            }

            if (lens[i].type == 6 || lens[i].type == 89 )
                NPRINTF(stderr, "beta %.3lf\n", lens[i].beta);

            if (lens[i].type == 7)
                NPRINTF(stderr, "masse %.3lf 10^12 Msol\n", lens[i].masse);

        }
    }
    else /* spline mapping */
    {
        sp_set(map_p, G.nx, G.ny, map_axx, map_ayy);
        NPRINTF(stderr, "WRITE: absolute potential map -> pot.best\n");
        wr_pot("pot.best", map_p);
        NPRINTF(stderr, "WRITE: absolute mass map -> mass.best\n");
        wr_mass("mass.best", map_axx, map_ayy);
    }

    // Print the optimized redshifts
    for (i = 0; i < I.nzlim; i++)
    {
        splitzmlimit(zlim[i].n, limages);
        j = 0;
        while ( indexCmp( multi[j][0].n, limages[0] ) ) j++;
        z_dlsds = multi[j][0].dr;
        multi[j][0].z = zero(0.1, 100., fz_dlsds);
        NPRINTF(stderr, "#%s z:%.3lf dlsds:%.3lf \n",
                zlim[i].n, multi[j][0].z, multi[j][0].dr);
    }

    NPRINTF(stderr, "***********************************************************\n");

    /* ecriture des resultats dans best.par */

    best = fopen("best.par", "w");

    time( &rawtime );
    fprintf(best, "#%s\n", asctime( localtime ( &rawtime ) ) );
    if ( I.forme >= 0 )
        fprintf(best, "#Source plane optimization\n");
    else
        fprintf(best, "#Image plane optimization\n");

    fprintf(best, "#Chi2tot(dof=%d): %.4lf\n", getNConstraints() - getNParameters(), chi0);
    fprintf(best, "#Chi2pos: %.3lf\n", chip);
	fprintf(best, "#Chi2_vel: %.3lf\n",chi_vel);      //THIS IS MINE
	fprintf(best, "#Chi2_mass: %.3lf\n",chi_mass);      //THIS IS MINE
    fprintf(best, "#Chi2formex: %.3lf\n", chix);
    fprintf(best, "#Chi2formey: %.3lf\n", chiy);
    fprintf(best, "#Chi2l: %.3lf\n", chil);
    if ( M.inverse == 3 )
        fprintf( best, "#log(Evidence) : %.3lf\n", evidence);

    // REDSHIFTS
    for (i = 0; i < I.nzlim; i++)
    {
        splitzmlimit( zlim[i].n, limages);
        j = 0;
        while ( indexCmp( multi[j][0].n, limages[0] ) ) j++;
        fprintf(best, "#%s z:%.3lf dlsds:%.3lf \n",
                zlim[i].n, multi[j][0].z, multi[j][0].dr);
    };

    fprintf(best, "#n_Warning: %d\n", nwarn);

    // RUNMODE
    fprintf(best, "runmode\n");
    fprintf(best, "\treference     %d %lf %lf\n", M.iref, M.ref_ra, M.ref_dec);
    fprintf(best, "\timage     %d %s\n", M.image, M.imafile);
    fprintf(best, "\tsource    %d %s\n", M.source, M.sourfile);
    if (M.study)
        fprintf(best, "\tstudy     %d %s\n", M.study, M.studyfile);
    if (M.seeing)
        fprintf(best, "\timseeing  %.3lf\n", M.seeing);
    if (M.imass)
        fprintf(best, "\tmass      %d %d %.3lf %s\n", M.imass, M.nmass, M.zmass, M.massfile);
    if (M.iampli)
        fprintf(best, "\tampli	   %d %d %.3lf %s\n", M.iampli, M.nampli, M.zampli, M.amplifile);
    if (M.ishear)
        fprintf(best, "\tshear  %d %d %.3lf %s\n", M.ishear, M.nshear, M.zshear, M.shearfile);
    if (M.ishearf)
        fprintf(best, "\tshearfield  %d %.3lf %s %d\n", M.ishearf, M.zshearf, M.shearffile, M.nshearf);
    if (M.grille)
        fprintf(best, "\tgrille   %d %d %lf\n", M.grille, M.ngrille, M.zgrille);
    if (M.pixel)
        fprintf(best, "\tpixel     %d %d %s\n", M.pixel, M.npixel, M.pixelfile);

    fprintf(best, "\tend\n");

    // GRILLE
    fprintf(best, "grille\n");
    fprintf(best, "\tnombre      %d\n", G.ngrid);
    fprintf(best, "\tpolaire     %d\n", 0);
    fprintf(best, "\tnlentille   %ld\n", G.nlens);
    if ( strcmp(CL.algorithm, "MARCHINGSQUARES") )
        fprintf(best, "\tnlens_crit   %ld\n", G.nlens_crit);
    if ( sigposAs.bk != 0 )
    {
        for ( i = 0; i < I.n_mult; i++)
            for ( j = 0; j < I.mult[i]; j++)
                fprintf( best, "\tsigposAs  %s  %lf\n", multi[i][0].n, I.sig2pos[i][j]);
    }
    fprintf(best, "\tend\n");

    // SOURCE
    if (M.image != 0 || M.source != 0)
    {
        fprintf(best, "source\n");
        fprintf(best, "\tz_source     %.3lf\n", S.zs);
        fprintf(best, "\tend\n");
    };
    
    // IMAGE 
    if (I.nzlim != 0 || I.zarclet > 0)
    {
        fprintf(best, "image\n");
        for( i = 0; i < I.nzlim; i++ )
        {
            splitzmlimit( zlim[i].n, limages);
            j = 0; while ( indexCmp( multi[j][0].n, limages[0] ) ) j++;
            fprintf(best, "\tz_m_limit %d %s %d %.3lf  %.3lf  %.4lf \n", 1, zlim[i].n, 0, multi[j][0].z, 0., 0.);
        }

        if( I.zarclet > 0 )
            fprintf(best, "\tz_arclet %lf\n", I.zarclet);

        fprintf(best, "\tend\n");
    }

    // CLEANLENS
    if (M.iclean != 0)
    {
        fprintf(best, "cleanlens\n");
        fprintf(best, "\tcleanset  %d  %f\n", M.iclean, M.zclean);
        if(strcmp(imFrame.pixfile, ""))
            fprintf(best, "\timframe %d  %s\n", imFrame.format, imFrame.pixfile);
        if(strcmp(ps.pixfile, "")) 
            fprintf(best, "\tsframe  %s\n", ps.pixfile);
        if(strcmp(M.centerfile, ""))
            fprintf(best, "\tc_image  %s\n", M.centerfile);
        if(imFrame.ncont > 0)
        {
            fprintf(best, "\tncont  %d  %s\n", imFrame.ncont, imFrame.outfile);
            for(i = 0; i < imFrame.ncont; i++)
                fprintf(best, "\tcontour\t%d %s\n", i+1, imFrame.contfile[i]);
        }
        if( imFrame.column != 1 )
           fprintf(best, "\tcolumn  %d\n", imFrame.column); 
        fprintf(best, "\techant\t%d\n", imFrame.ech);
        fprintf(best, "\ts_echant\t%d\n", ps.ech);
        fprintf(best, "\ts_n\t%d\n", ps.nx);
        if( imFrame.header != 0 )
            fprintf(best, "\theader\t%d\n", imFrame.header);
        fprintf(best, "\tpixelx\t%lf\n", imFrame.pixelx);
        fprintf(best, "\tpixely\t%lf\n", imFrame.pixely);
        fprintf(best, "\txmin\t%lf\n", imFrame.xmin);
        fprintf(best, "\tymin\t%lf\n", imFrame.ymin);
        fprintf(best, "\ts_xmin\t%lf\n", ps.xmin);
        fprintf(best, "\ts_ymin\t%lf\n", ps.ymin);
        fprintf(best, "\ts_xmax\t%lf\n", ps.xmax);
        fprintf(best, "\ts_ymax\t%lf\n", ps.ymax);
        fprintf(best, "\tend\n");
    }

    // Write all the potentials with arcsec and kpc values
    for (i = 0; i < G.nlens; i++)
        writePotentiel(best, i, 3);

    // CLINE
    fprintf(best, "cline\n");
    fprintf(best, "\tnplan    %d", CL.nplan);
    for (i = 0; i < CL.nplan; i++)
        fprintf(best, " %.3lf ", CL.cz[i]);

    fprintf(best, "\n");
    fprintf(best, "\tdmax     %.3lf\n", CL.dmax);
    fprintf(best, "\talgorithm   %s\n", CL.algorithm);
    if ( !strcmp(CL.algorithm, "MARCHINGSQUARES") )
    {
        fprintf(best, "\tlimitHigh   %.1lf\n", CL.limitHigh);
        fprintf(best, "\tlimitLow    %.3lf\n", CL.cpas);
    }
    else
        fprintf(best, "\tpas      %.3lf\n", CL.cpas);
    fprintf(best, "\tend\n");

	//DYNFILE          
    if ( Dy.dyntype != 0 )
    {
    	fprintf(best,"dynfile\n");
    	fprintf(best,"\tdyntype    %d\n",Dy.dyntype);
    	fprintf(best,"\tdynnumber    %d\n",Dy.dynnumber);
    	fprintf(best,"\tvelocity    %lf  \n",Dy.dynvel);
    	fprintf(best,"\te_velocity    %lf  \n",Dy.dynevel);
    	fprintf(best,"\tindependent mass    %.3le  \n",Dy.indmass);
    	fprintf(best,"\tindependent e_mass    %.3le  \n",Dy.indemass);
    	fprintf(best,"\treference radius kpc    %lf  \n",Dy.refradius);
    	fprintf(best,"\tend\n");
    }
	
    // GRANDE
    fprintf(best, "grande\n");
    fprintf(best, "\tiso         %d %d %.3lf %.3lf %.3lf\n", L.iso, L.nmaxiso, L.scale, L.zonex, L.zoney);
    fprintf(best, "\tname        best\n");
    fprintf(best, "\tprofil      %d %d\n", L.profil, L.pt);
    fprintf(best, "\tcontour     %d %d\n", L.ncourbe, L.pt);
    fprintf(best, "\tlarge_dist  %.3lf\n", L.dlarge);
    fprintf(best, "\tend\n");

    // OBSERVATION
    if (M.pixel || M.iclean)
    {
        fprintf(best, "observation\n");
        if( O.setseeing == 1 )
            fprintf(best, "\tseeing       %d %lf\n", O.setseeing, O.seeing);
        else if( O.setseeing == 2)
            fprintf(best, "\tseeing_e      %d %lf %lf %lf\n", O.setseeing, O.seeing_a, O.seeing_b, O.seeing_angle);
        else if( O.setseeing == 3)
            fprintf(best, "\tpsf      %d %s\n", O.setseeing, O.psffile);

        fprintf(best, "\tbinning      %d %d\n", O.setbin, O.bin);
        fprintf(best, "\tbruit        %d\n", O.bruit);
        fprintf(best, "\tSKY          %.3lf\n", O.SKY);
        if( O.gain > 0 )
            fprintf(best, "\tdispersion   %.3lf\n", sqrt(O.SKY / O.gain));

        fprintf(best, "\tidum         %d\n", O.idum);
        fprintf(best, "\tend\n");
    };
 
    // VELOCITY FIELD
    if(M.cube || (M.iclean==3))
    {
      fprintf(best, "vfield\n");
      fprintf(best, "\tprofile   %d\n", vf.profile);
      fprintf(best, "\tx_centre  %.3lf\n", vf.C.x);
      fprintf(best, "\ty_centre  %.3lf\n", vf.C.x);
      fprintf(best, "\tvt        %.3lf\n", vf.vt);
      fprintf(best, "\trt        %.3lf\n", vf.rt);
      fprintf(best, "\ti         %.3lf\n", vf.i);
      fprintf(best, "\ttheta     %.3lf\n", vf.theta);
      fprintf(best, "\tlcent     %.3lf\n", vf.lcent);
      fprintf(best, "\tsigma     %.3lf\n", vf.sigma);
      fprintf(best, "\tend\n");
    }

    // COSMOLOGY
    fprintf(best, "cosmologie\n");
    fprintf(best,"\tmodel       %d\n",C.model);
    fprintf(best, "\tH0        %.3lf\n", C.H0);
    fprintf(best, "\tomegaM    %.3lf\n", C.omegaM);
    fprintf(best, "\tomegaX    %.3lf\n", C.omegaX);
    if ( C.kcourb == 0. ) fprintf(best, "\tomegaK    0.\n");
    fprintf(best, "\twX        %.3lf\n", C.wX);
    fprintf(best, "\twa        %.3lf\n", C.wa);
    fprintf(best, "\tend\n");

    // CHAMP
    fprintf(best, "champ\n");
    fprintf(best, "\txmin     %.3lf\n", F.xmin);
    fprintf(best, "\txmax     %.3lf\n", F.xmax);
    fprintf(best, "\tymin     %.3lf\n", F.ymin);
    fprintf(best, "\tymax     %.3lf\n", F.ymax);
    if(F.lmin>0)
    {
       fprintf(best, "\tlmin     %.3lf\n", F.lmin);
       fprintf(best, "\tlmax     %.3lf\n", F.lmax);
    }
    fprintf(best, "\tend\n");


    fprintf(best, "fini\n");
    fclose(best);


    /******************************************************************
     *  Write the bestopt.par file
     ******************************************************************/

    if ( M.inverse < 3 )
    {
        for (i = 0; i < G.no_lens; i++)
        {
            lmin[i] = lmin_s[i];
            lmax[i] = lmax_s[i];
            prec[i] = prec_s[i];
            for (j = 0; j < NPAMAX; j++)
                block[i][j] = block_s[i][j];
        };

        for (i = 0; i < I.nzlim; i++)
            zlim_s[i] = zlim[i];
    }

    besto = fopen("bestopt.par", "w");

    fprintf(besto, "#Chi2tot: %.3lf\n", chi0);
    fprintf(besto, "#Chi2pos: %.3lf\n", chip);
	fprintf(besto, "#Chi2_vel: %.3lf\n",chi_vel);      //THIS IS MINE     TV
	fprintf(besto, "#Chi2_mass: %.3lf\n",chi_mass);      //THIS IS MINE   TV
    fprintf(besto, "#Chi2formex: %.3lf\n", chix);
    fprintf(besto, "#Chi2formey: %.3lf\n", chiy);
    fprintf(besto, "#Chi2l: %.3lf\n", chil);
    if ( M.inverse == 3 )
        fprintf( besto, "#log(Evidence) : %.3lf\n", evidence);

    // RUNMODE
    fprintf(besto, "runmode\n");
    fprintf(besto, "\treference     %d %lf %lf\n", M.iref, M.ref_ra, M.ref_dec);
    fprintf(besto, "\timage     %d %s\n", M.image, M.imafile);
    if( strcmp(M.sourfile, "source.best") )
        fprintf(besto, "\tsource    %d %s\n", M.source, M.sourfile);

    if (M.study)
        fprintf(besto, "\tstudy     %d %s\n", M.study, M.studyfile);
    if (M.seeing)
        fprintf(besto, "\timseeing  %.3lf\n", M.seeing);
    if ( M.inverse < 3 )
        fprintf(besto, "\tinverse   %d %d\n", M.inverse, M.itmax);
    else
        fprintf(besto, "\tinverse   %d %lf %d\n", M.inverse, M.rate, M.itmax );

    if (M.imass)
        fprintf(besto, "\tmass      %d %d %.3lf %s\n", M.imass, M.nmass, M.zmass, M.massfile);
    if (M.iampli)
        fprintf(besto, "\tampli	   %d %d %.3lf %s\n", M.iampli, M.nampli, M.zampli, M.amplifile);
    if (M.ishear)
        fprintf(besto, "\tshear  %d %d %.3lf %s\n", M.ishear, M.nshear, M.zshear, M.shearfile);
    if (M.ishearf)
        fprintf(besto, "\tshearfield  %d %.3lf %s %d\n", M.ishearf, M.zshearf, M.shearffile, M.nshearf);
    if (M.grille)
        fprintf(besto, "\tgrille   %d %d %lf\n", M.grille, M.ngrille, M.zgrille);
    if (M.pixel)
        fprintf(besto, "\tpixel     %d %d %s\n", M.pixel, M.npixel, M.pixelfile);

    fprintf(besto, "\tend\n");

    // IMAGE
    fprintf(besto, "image\n");
    fprintf(besto, "\tmultfile    %d %s\n", I.n_mult, I.multfile);
    fprintf(besto, "\tforme       %d\n", I.forme);
    if (I.stat > 0)
        fprintf(besto, "\tarcletstat  %d %d %s\n", I.stat, I.statmode, I.arclet);

    if (I.srcfit > 0)
        fprintf(besto, "\tsourcefit   %d %s %s\n", I.srcfit, I.srcfitFile, I.srcfitMethod);

    if (I.npcl > 0)
        for (i = 0; i < I.npcl; i++)
            fprintf(besto, "\tcritic    %d %.3lf %.3lf %.3lf %.3lf %.3lf\n",
                    cl[i].n, cl[i].C.x, cl[i].C.y, cl[i].phi / DTR, cl[i].dl, cl[i].z);

    if (I.nzlim > 0)
        for (i = 0; i < I.nzlim; i++)
            fprintf(besto, "\tz_m_limit  %d %s %d %.3lf %.3lf %.4lf \n", i + 1,
                    zlim[i].n, zlim[i].bk, zlim[i].min, zlim[i].max, zlim[i].dderr);

    if (I.zarclet > 0)
        fprintf(besto, "\tz_arclet    %lf\n", I.zarclet);

    if (zalim.bk > 0)
        fprintf(besto, "\tz_a_limit    %d %lf %lf\n", zalim.bk, zalim.min, zalim.max);

    if (I.mult_abs > 0)
        fprintf(besto, "\tmult_wcs    %d\n", I.mult_abs);

    if ( sigposAs.bk == 0 )
        fprintf(besto, "\tsigposArcsec   %lf\n", sigposAs.min);
    else
        fprintf(besto, "\tsigposArcsec   %d %lf %lf\n", sigposAs.bk, sigposAs.min, sigposAs.max);

    fprintf(besto, "\tend\n");

    // GRILLE
    fprintf(besto, "grille\n");
    fprintf(besto, "\tnombre      %d\n", G.ngrid);
    fprintf(besto, "\tpolaire     %d\n", 0);
    fprintf(besto, "\tnlentille   %ld\n", G.nlens);
    if ( strcmp(CL.algorithm, "MARCHINGSQUARES") )
        fprintf(besto, "\tnlens_crit   %ld\n", G.nlens_crit);

    fprintf(besto, "\tnlens_opt   %ld\n", G.no_lens);
    fprintf(besto, "\tend\n");

    // SOURCE
    if (M.image != 0 || M.source != 0)
    {
        fprintf(besto, "source\n");
        fprintf(besto, "\tz_source     %.3lf\n", S.zs);
        fprintf(besto, "\tend\n");
    };

    // CLEANLENS
    if (M.iclean != 0)
    {
        fprintf(besto, "cleanlens\n");
        fprintf(besto, "\tcleanset  %d  %f\n", M.iclean, M.zclean);
        if(strcmp(imFrame.pixfile, ""))
            fprintf(besto, "\timFrame %d  %s\n", imFrame.format, imFrame.pixfile);
        if(strcmp(ps.pixfile, "")) 
            fprintf(besto, "\tsframe  %s\n", ps.pixfile);
        if(strcmp(M.centerfile, ""))
            fprintf(besto, "\tc_image  %s\n", M.centerfile);
        if(imFrame.ncont > 0)
        {
            fprintf(besto, "\tncont  %d  %s\n", imFrame.ncont, imFrame.outfile);
            for(i = 0; i < imFrame.ncont; i++)
                fprintf(besto, "\tcontour\t%d %s\n", i+1, imFrame.contfile[i]);
        }
        if( imFrame.column != 1 )
           fprintf(besto, "\tcolumn  %d\n", imFrame.column); 
        fprintf(besto, "\techant\t%d\n", imFrame.ech);
        fprintf(besto, "\ts_echant\t%d\n", ps.ech);
        fprintf(besto, "\ts_n\t%d\n", ps.nx);
        if( imFrame.header != 0 )
            fprintf(besto, "\theader\t%d\n", imFrame.header);
        fprintf(besto, "\tpixelx\t%lf\n", imFrame.pixelx);
        fprintf(besto, "\tpixely\t%lf\n", imFrame.pixely);
        fprintf(besto, "\txmin\t%lf\n", imFrame.xmin);
        fprintf(besto, "\tymin\t%lf\n", imFrame.ymin);
        fprintf(besto, "\ts_xmin\t%lf\n", ps.xmin);
        fprintf(besto, "\ts_ymin\t%lf\n", ps.ymin);
        fprintf(besto, "\ts_xmax\t%lf\n", ps.xmax);
        fprintf(besto, "\ts_ymax\t%lf\n", ps.ymax);
        fprintf(besto, "\tend\n");
    }

    if ( M.iclean == 2 )
        for( i = 0; i < S.ns; i++ )
        {
            // SHAPE MODEL
            fprintf(besto, "shapemodel\n");
            fprintf(besto, "\tid   %s\n",  source[i].n);
            fprintf(besto, "\ts_center_x   %.3lf\n",  source[i].C.x);
            fprintf(besto, "\ts_center_y   %.3lf\n",  source[i].C.y);
            fprintf(besto, "\ts_angle      %.3lf\n",  source[i].E.theta * RTD);
            fprintf(besto, "\ts_sigx       %.3lf\n",  source[i].E.a);
            fprintf(besto, "\ts_sigy       %.3lf\n",  source[i].E.b);
            fprintf(besto, "\tmag          %.3lf\n",  source[i].mag);
            fprintf(besto, "\tend\n");

            // SHAPE LIMITS
            fprintf(besto, "shapelimit\n");
            if( sblock[i][SCX] )
                fprintf(besto, "\ts_center_x   %d %.3lf %.3lf\n", sblock[i][SCX], smin[i].C.x, smax[i].C.x);
            if( sblock[i][SCY] )
                fprintf(besto, "\ts_center_y   %d %.3lf %.3lf\n", sblock[i][SCY], smin[i].C.y, smax[i].C.y);
            if( sblock[i][STHETA] )
                fprintf(besto, "\ts_angle      %d %.3lf %.3lf\n", sblock[i][STHETA], smin[i].E.theta * RTD, smax[i].E.theta * RTD);
            if( sblock[i][SA] )
                fprintf(besto, "\ts_sigx       %d %.3lf %.3lf\n", sblock[i][SA], smin[i].E.a, smax[i].E.a);
            if( sblock[i][SB] )
             fprintf(besto, "\ts_sigy       %d %.3lf %.3lf\n", sblock[i][SB], smin[i].E.b, smax[i].E.b);
            if( sblock[i][SFLUX] )
                fprintf(besto, "\tmag          %d %.3lf %.3lf\n", sblock[i][SFLUX], smin[i].mag, smax[i].mag);
            fprintf(besto, "\tend\n");
        }

    // Write the optimized potentials and limits
    for ( i = 0 ; i < G.nplens[0] ; i++ )
    {
        writePotentiel(besto, i, 2);
        writeLimit(besto, i);
    }

    // POTFILE
    for( i = 0; i < G.npot; i++ )
        if ( P[i].ftype != 0)
        {
            fprintf(besto, "potfile%d\n", i);
            fprintf(besto, "\tfilein  %d %s\n", P[i].ftype, P[i].potfile );
            fprintf(besto, "\tzlens   %lf\n", P[i].zlens );
            fprintf(besto, "\ttype    %d\n", P[i].type );
            if ( P[i].corekpc != -1 )
                fprintf(besto, "\tcorekpc %lf\n", P[i].corekpc );
            else
                fprintf(besto, "\tcore %lf\n", P[i].core );
    
            fprintf(besto, "\tmag0    %lf\n", P[i].mag0 );
            fprintf(besto, "\tsigma   %d %lf %lf\n", P[i].isigma, P[i].sigma1, P[i].sigma2 );
            if ( P[i].cutkpc1 != DBL_MAX )
                fprintf(besto, "\tcutkpc  %d %lf %lf\n", P[i].ircut, P[i].cutkpc1, P[i].cutkpc2 );
            else
                fprintf(besto, "\tcut     %d %lf %lf\n", P[i].ircut, P[i].cut1, P[i].cut2 );
    
            if ( P[i].ftype == 62 )
            {
                fprintf(besto, "\tm200slope   %d %lf %lf\n", P[i].islope, P[i].slope1, P[i].slope2 );
                fprintf(besto, "\tc200slope %d %lf %lf\n", P[i].ivdslope, P[i].vdslope1, P[i].vdslope2 );
                fprintf(besto, "\tm200   %d %lf %lf\n", P[i].ia, P[i].a1, P[i].a2 );
                fprintf(besto, "\tc200   %d %lf %lf\n", P[i].ib, P[i].b1, P[i].b2 );
            }
            else
            {
                fprintf(besto, "\tslope   %d %lf %lf\n", P[i].islope, P[i].slope1, P[i].slope2 );
                fprintf(besto, "\tvdslope %d %lf %lf\n", P[i].ivdslope, P[i].vdslope1, P[i].vdslope2 );
                fprintf(besto, "\tvdscatter %d %lf %lf\n", P[i].ivdscat, P[i].vdscat1, P[i].vdscat2 );
                fprintf(besto, "\trcutscatter %d %lf %lf\n", P[i].ircutscat, P[i].rcutscat1, P[i].rcutscat2 );
            }
            fprintf(besto, "\tend\n");
        }

    // Write the grid potentials
    for( i = G.nmsgrid; i < G.nlens; i++ )
    {
        writePotentiel(besto, i, 2);
    //    writeLimit(besto, i);
    }

    // CLINE
    fprintf(besto, "cline\n");
    fprintf(besto, "\tnplan    %d", CL.nplan);
    for (i = 0; i < CL.nplan; i++)
        fprintf(besto, " %.3lf ", CL.cz[i]);

    fprintf(besto, "\n");
    fprintf(besto, "\tdmax     %.3lf\n", CL.dmax);
    fprintf(besto, "\talgorithm   %s\n", CL.algorithm);
    if ( !strcmp(CL.algorithm, "MARCHINGSQUARES") )
    {
        fprintf(besto, "\tlimitHigh   %.1lf\n", CL.limitHigh);
        fprintf(besto, "\tlimitLow    %.3lf\n", CL.cpas);
    }
    else
        fprintf(besto, "\tpas      %.3lf\n", CL.cpas);

    fprintf(besto, "\tend\n");

	//DYNFILE   
    if ( Dy.dyntype != 0 )
    {
	    fprintf(besto,"dynfile\n");
	    fprintf(besto,"\tdyntype    %d\n",Dy.dyntype);
	    fprintf(besto,"\tdynnumber    %d\n",Dy.dynnumber);
	    fprintf(besto,"\tvelocity    %lf  \n",Dy.dynvel);
	    fprintf(besto,"\te_velocity    %lf  \n",Dy.dynevel);
	    fprintf(besto,"\tindependent mass    %.3le  \n",Dy.indmass);
	    fprintf(besto,"\tindependent e_mass    %.3le  \n",Dy.indemass);
	    fprintf(besto,"\treference radius kpc    %lf  \n",Dy.refradius);
	    fprintf(besto,"\tend\n");
    }

    // GRANDE
    fprintf(besto, "grande\n");
    fprintf(besto, "\tiso         %d %d %.3lf %.3lf %.3lf\n", L.iso, L.nmaxiso, L.scale, L.zonex, L.zoney);
    fprintf(besto, "\tname        best\n");
    fprintf(besto, "\tprofil      %d %d\n", L.profil, L.pt);
    fprintf(besto, "\tcontour     %d %d\n", L.ncourbe, L.pt);
    fprintf(besto, "\tlarge_dist  %.3lf\n", L.dlarge);
    fprintf(besto, "\tend\n");

    // OBSERVATIONS
    if (M.pixel || M.iclean)
    {
        fprintf(besto, "observation\n");
        if( O.setseeing == 1 )
            fprintf(besto, "\tseeing       %d %lf\n", O.setseeing, O.seeing);
        else if( O.setseeing == 2)
            fprintf(besto, "\tseeing_e      %d %lf %lf %lf\n", O.setseeing, O.seeing_a, O.seeing_b, O.seeing_angle);
        else if( O.setseeing == 3)
            fprintf(besto, "\tpsf      %d %s\n", O.setseeing, O.psffile);

        fprintf(besto, "\tbinning      %d %d\n", O.setbin, O.bin);
        fprintf(besto, "\tbruit        %d\n", O.bruit);
        fprintf(besto, "\tSKY          %.3lf\n", O.SKY);
        if( O.gain > 0 )
            fprintf(besto, "\tdispersion   %.3lf\n", sqrt(O.SKY / O.gain));

        fprintf(besto, "\tidum         %d\n", O.idum);
        fprintf(besto, "\tend\n");
    };
    // VELOCITY FIELD

    if(M.cube || (M.iclean==3))
    {
      fprintf(besto, "vfield\n");
      fprintf(besto, "\tprofile   %d\n", vf.profile);
      fprintf(besto, "\tx_centre  %.3lf\n", vf.C.x);
      fprintf(besto, "\ty_centre  %.3lf\n", vf.C.x);
      fprintf(besto, "\tvt        %.3lf\n", vf.vt);
      fprintf(besto, "\trt        %.3lf\n", vf.rt);
      fprintf(besto, "\ti         %.3lf\n", vf.i*RTD);
      fprintf(besto, "\ttheta     %.3lf\n", vf.theta*RTD);
      fprintf(besto, "\tlcent     %.3lf\n", vf.lcent);
      fprintf(besto, "\tsigma     %.3lf\n", vf.sigma);
      fprintf(besto, "\tend\n");

      // VELOCITY FIELD LIMITS
      fprintf(besto, "vfieldlimit\n");
      fprintf(besto, "\tx_centre  %d %.3lf %.3lf\n", vfblock[VFCX],vfmin.C.x,vfmax.C.x);
      fprintf(besto, "\ty_centre  %d %.3lf %.3lf\n", vfblock[VFCY],vfmin.C.y,vfmax.C.y);
      fprintf(besto, "\tvt        %d %.3lf %.3lf\n", vfblock[VFVT],vfmin.vt,vfmax.vt);
      fprintf(besto, "\trt        %d %.3lf %.3lf\n", vfblock[VFRT],vfmin.rt,vfmax.rt);
      fprintf(besto, "\ti         %d %.3lf %.3lf\n", vfblock[VFI],vfmin.i*RTD,vfmax.i*RTD);
      fprintf(besto, "\ttheta     %d %.3lf %.3lf\n", vfblock[VFTHETA],vfmin.theta*RTD,vfmax.theta*RTD);
      fprintf(besto, "\tlcent     %d %.3lf %.3lf\n", vfblock[VFLCENT],vfmin.lcent,vfmax.lcent);
      fprintf(besto, "\tsigma     %d %.3lf %.3lf\n", vfblock[VFSIGMA],vfmin.sigma,vfmax.sigma);
      fprintf(besto, "\tend\n");
    }

    // COSMOLOGY
    fprintf(besto, "cosmologie\n");
    fprintf(besto,"\tmodel       %d\n",C.model);
    fprintf(besto, "\tH0         %.3lf\n", C.H0);
    fprintf(besto, "\tomegaM     %.3lf\n", C.omegaM);
    fprintf(besto, "\tomegaX     %.3lf\n", C.omegaX);
    if ( C.kcourb == 0. ) fprintf(besto, "\tomegaK     0.\n");
    fprintf(besto, "\twX         %.3lf\n", C.wX);
    fprintf(besto, "\twa         %.3lf\n", C.wa);
    fprintf(besto, "\tend\n");

    // COSMOLIMITS
    fprintf(besto, "cosmolimit\n");
    fprintf(besto, "\tomegaM		%d %.3lf %.3lf\n",
            cblock[OMEGAM], clmin.omegaM, clmax.omegaM);
    fprintf(besto, "\tomegaX		%d %.3lf %.3lf\n",
            cblock[OMEGAX], clmin.omegaX, clmax.omegaX);
    fprintf(besto, "\twX		%d %.3lf %.3lf\n",
            cblock[WX], clmin.wX, clmax.wX);
    fprintf(besto, "\twa 	%d %.3lf %.3lf\n",
            cblock[WA], clmin.wa, clmax.wa);
    fprintf(besto, "\tend\n");

    // CHAMP
    fprintf(besto, "champ\n");
    fprintf(besto, "\txmin     %.3lf\n", F.xmin);
    fprintf(besto, "\txmax     %.3lf\n", F.xmax);
    fprintf(besto, "\tymin     %.3lf\n", F.ymin);
    fprintf(besto, "\tymax     %.3lf\n", F.ymax);
    fprintf(besto, "\tend\n");
    if(F.lmin>0)
    {
       fprintf(besto, "\tlmin     %.3lf\n", F.lmin);
       fprintf(besto, "\tlmax     %.3lf\n", F.lmax);
    }
    fprintf(besto, "\tend\n");

    fprintf(besto, "fini\n");
    fclose(besto);

    // Write the arclet.best file
    if (I.stat == 1)
    {
        for (i = 0; i < S.ns; i++)
        {
            z_dlsds = arclet[i].dr;
            arclet[i].z = zero(lens[0].z, 100., fz_dlsds);
        };
        ecrire_r(0, S.ns, arclet, "arclet.best", 1);
    }
}


/* Write a potentiel section in best file 
 * flag = 1 : print values in arcsec
 * flag = 2 : print values in kpc
 * flag = 3 : print values in arcsec and kpc
 */
static void writePotentiel(FILE *best, long int i, int flag)
{
    extern struct pot lens[];

    fprintf(best, "potentiel %s\n", lens[i].n);
    fprintf(best, "\tprofil       %d\n", lens[i].type);
    if (lens[i].type == 9 )
    {
        fprintf(best, "\trhos     %.4lf\n", lens[i].pmass);
        fprintf(best, "\tz_lens     %.4lf\n", lens[i].z);
        fprintf(best, "\tend\n");
        return;
    }
    if (lens[i].type == 14 )
    {
        fprintf(best, "\tgamma     %.4lf\n", lens[i].emass);
        fprintf(best, "\tangle_pos     %.4lf\n", lens[i].theta*RTD);
        fprintf(best, "\tz_lens     %.4lf\n", lens[i].z);
        fprintf(best, "\tend\n");
        return;
    }
    fprintf(best, "\tx_centre     %.3lf\n", lens[i].C.x);
    fprintf(best, "\ty_centre     %.3lf\n", lens[i].C.y);
    if (lens[i].type != 0 && lens[i].type != 2 && lens[i].type != 7)
    {
        fprintf(best, "\tellipticite     %.3lf\n", lens[i].emass);
        if ( lens[i].type == 121 )
        {
		    fprintf(best,"\ttheta       %.3lf\n",lens[i].theta*RTD);  // triaxial NFW
		    fprintf(best,"\tphi       %.3lf\n",lens[i].phi*RTD);
        }
        else
			fprintf(best,"\tangle_pos       %.3lf\n",lens[i].theta*RTD);
    }

    if ( lens[i].type == 12 )
    {
        if( flag & 1 ) fprintf(best, "\tscale_radius      %.3lf\n", lens[i].rc);
        if( flag & 2 ) fprintf(best, "\tscale_radius_kpc  %.3lf\n", lens[i].rckpc);

        if ( lens[i].rcut != DBL_MAX )
        {
            if( flag & 1 ) fprintf(best, "\tr200               %.3lf\n", lens[i].rcut);
            if( flag & 2 ) fprintf(best, "\tr200_kpc           %.3lf\n", lens[i].rcutkpc);
        }

        fprintf(best, "\tv_disp          %.3lf\n", lens[i].sigma);

        fprintf(best, "\tconcentration   %.3lf\n", lens[i].beta);
        fprintf(best, "\tm200            %.3le\n", lens[i].masse);
        fprintf(best, "\trhos            %.3le\n", lens[i].pmass);
        fprintf(best, "\trc_slope        %.3lf\n", lens[i].rcslope);/////////////////
        fprintf(best, "\talpha           %.3lf\n", lens[i].alpha);
    }
    else if( lens[i].type == 16 )
    {
        if( flag & 1 ) fprintf(best, "\tscale_radius      %.3lf\n", lens[i].rc);
        if( flag & 2 ) fprintf(best, "\tscale_radius_kpc  %.3lf\n", lens[i].rckpc);
        fprintf(best, "\tv_disp          %.3lf\n", lens[i].sigma);
    }
    else
    {
        if (lens[i].type != 0 && fabs(lens[i].type) != 1 && lens[i].type != 7)
        {
            if( flag & 1 ) fprintf(best, "\tcore_radius         %.3lf\n", lens[i].rc);
            if( flag & 2 ) fprintf(best, "\tcore_radius_kpc     %.3lf\n", lens[i].rckpc );
        }

        if ( lens[i].rcut != DBL_MAX )
        {
            if( flag & 1 ) fprintf(best, "\tcut_radius         %.3lf\n", lens[i].rcut);
            if( flag & 2 ) fprintf(best, "\tcut_radius_kpc     %.3lf\n", lens[i].rcutkpc );
        }

        if ( lens[i].type == 13 )
        {
            fprintf(best, "\tsigma_e     %.3le\n", lens[i].sigma);
            fprintf(best, "\tn     %.3lf\n", lens[i].alpha);
        }
        else
            fprintf(best, "\tv_disp     %.3lf\n", lens[i].sigma);

        if ( lens[i].type == 3 || lens[i].type == 6 || lens[i].type == 84 ||
             lens[i].type == 87 || lens[i].type == 88 )
            fprintf(best, "\talpha     %.3lf\n", lens[i].alpha);

        if (lens[i].type == 6 || lens[i].type == 89 )
        {
            fprintf(best, "\tbeta     %.3lf\n", lens[i].beta);
            fprintf(best, "\trc_slope     %.3lf\n", lens[i].rcslope);
        }

        if (lens[i].type == 7)
            fprintf(best, "\tmasse     %.3lf\n", lens[i].masse);
    }

    if (lens[i].mag != 0)
        fprintf(best, "\tmag		  %.3lf\n", lens[i].mag);

    fprintf(best, "\tz_lens     %.4lf\n", lens[i].z);
    fprintf(best, "\tend\n");
}

/* Write a Limit section */
static void writeLimit(FILE *best, long int i)
{
    extern struct g_grille   G;
    extern int block[][NPAMAX];
    extern struct pot lens[], lmin[], lmax[], prec[];

    fprintf(best, "limit %s\n", lens[i].n);
    if ( block[i][CX] != 0 )
        fprintf(best, "\tx_centre     %d %.3lf %.3lf %.3lf\n", block[i][CX],
                lmin[i].C.x, lmax[i].C.x, prec[i].C.x);
    if ( block[i][CY] != 0 )
        fprintf(best, "\ty_centre     %d %.3lf %.3lf %.3lf\n", block[i][CY],
                lmin[i].C.y, lmax[i].C.y, prec[i].C.y);
    if (lens[i].type != 0 && lens[i].type != 2 && lens[i].type != 7)
    {
        if ( block[i][EPOT] != 0 )
            fprintf(best, "\tellip_pot     %d %.3lf %.3lf %.3lf\n", block[i][EPOT], lmin[i].epot, lmax[i].epot, prec[i].epot);
        if ( block[i][EMASS] != 0 )
            fprintf(best, "\tellipticite     %d %.3lf %.3lf %.3lf\n", block[i][EMASS],
                    lmin[i].emass, lmax[i].emass, prec[i].epot);
        if ( block[i][THETA] != 0 )
            fprintf(best, "\tangle_pos     %d %.3lf %.3lf %.3lf\n", block[i][THETA],
                    lmin[i].theta*RTD, lmax[i].theta*RTD, prec[i].theta*RTD);
        if ( block[i][PHI] != 0 )
            fprintf(best, "\tphi     %d %.3lf %.3lf %.3lf\n", block[i][PHI],
                    lmin[i].phi*RTD, lmax[i].phi*RTD, prec[i].phi*RTD);
    }

    if ( lens[i].type == 12 || lens[i].type == 16 )
    {

        if ( block[i][RC] != 0 )
            fprintf(best, "\tscale_radius_kpc     %d %.3lf %.3lf %.3lf\n",
                    block[i][RC], lmin[i].rckpc, lmax[i].rckpc, prec[i].rckpc);

        if ( lens[i].rcut != DBL_MAX && block[i][RCUT] != 0 )
            fprintf(best, "\tvirial_radius_kpc     %d %.3lf %.3lf %.3lf\n",
                    block[i][RCUT], lmin[i].rcutkpc, lmax[i].rcutkpc, prec[i].rcutkpc );

        if ( block[i][BETA] != 0 )
            fprintf(best, "\tconcentration    %d %.3lf %.3lf %.3lf\n",
                    block[i][BETA], lmin[i].beta, lmax[i].beta, prec[i].beta);

        if ( block[i][MASSE] != 0 )
            fprintf(best, "\tvirial_mass    %d %.3le %.3le %.3le\n",
                    block[i][MASSE], lmin[i].masse, lmax[i].masse, prec[i].masse);


        if ( block[i][PMASS] != 0 )
            fprintf(best, "\trhos    %d %.3le %.3le %.3le\n",
                    block[i][PMASS], lmin[i].pmass, lmax[i].pmass, prec[i].pmass);
		
		if ( block[i][RCSLOPE] != 0 )
            fprintf(best, "\trc_slope    %d %.3le %.3le %.3le\n",
                    block[i][RCSLOPE], lmin[i].rcslope, lmax[i].rcslope, prec[i].rcslope);    ////
		
		
    }
    else
    {
        if ( block[i][RC] != 0 )
            fprintf(best, "\tcore_radius_kpc     %d %.3lf %.3lf %.3lf\n",
                    block[i][RC], lmin[i].rckpc, lmax[i].rckpc, prec[i].rckpc);

        if ( lens[i].rcut != DBL_MAX && block[i][RCUT] != 0 )
            fprintf(best, "\tcut_radius_kpc     %d %.3lf %.3lf %.3lf\n",
                    block[i][RCUT], lmin[i].rcutkpc, lmax[i].rcutkpc, prec[i].rcutkpc );


        if (lens[i].type == 6 || lens[i].type == 89 )
            if ( block[i][BETA] != 0 )
                fprintf(best, "\tbeta    %d %.3lf %.3lf %.3lf\n", block[i][BETA],
                        lmin[i].beta, lmax[i].beta, prec[i].beta);

        if (lens[i].type == 7 )
            if ( block[i][MASSE] != 0 )
                fprintf(best, "\tmasse    %d %.3le %.3le %.3le\n", block[i][MASSE],
                        lmin[i].masse, lmax[i].masse, prec[i].masse);
    }

    if ( block[i][B0] != 0 )
    {
        if ( lens[i].type == 13 )
            fprintf(best, "\tsigma_e    %d %.3le %.3le %.3le\n",
                    block[i][B0], lmin[i].sigma, lmax[i].sigma, prec[i].sigma);
        else
            fprintf(best, "\tv_disp    %d %.3lf %.3lf %.3lf\n",
                    block[i][B0], lmin[i].sigma, lmax[i].sigma, prec[i].sigma);
    }

    if ( i >= G.nmsgrid )
    {
        if ( block[i][B0] != 0 )
            fprintf(best, "\tv_disp    %d %.3lf %.3lf %.3lf\n",
                    block[i][B0], lmin[i].sigma, lmax[i].sigma, prec[i].sigma);
        else
            fprintf(best, "\trhos    %d %.3lf %.3lf %.3lf\n",
                    block[i][PMASS], lmin[i].pmass, lmax[i].pmass, prec[i].pmass);
    }


    if ( block[i][ALPHA] != 0 )
        fprintf(best, "\talpha     %d %.3lf %.3lf %.3lf\n",
                block[i][ALPHA], lmin[i].alpha, lmax[i].alpha, prec[i].alpha);

    if ( block[i][ZLENS] != 0 )
        fprintf(best, "\tz_lens     %d %.3lf %.3lf %.3lf\n",
                block[i][ZLENS], lmin[i].z, lmax[i].z, prec[i].z);

    if ( block[i][PMASS] != 0 && lens[i].type == 9 ) 
        fprintf(best, "\trhos     %d %.3lf %.3lf %.3lf\n",
                block[i][PMASS], lmin[i].pmass, lmax[i].pmass, prec[i].pmass);

    fprintf(best, "\tend\n");
}
