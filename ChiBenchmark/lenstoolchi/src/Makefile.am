# Program Makefile for lenstool
noinst_LIBRARIES = liblenstool.a

AM_CFLAGS 	= @CFITSIO_HDR@ @WCS_HDR@ @GSL_HDR@ -I@top_srcdir@/include \
		  -I@top_srcdir@/liblt 
#AM_CFLAGS 	= @CFITSIO_HDR@ @WCS_HDR@ -I@top_srcdir@/include \
#		  -I@top_srcdir@/liblt -I@top_srcdir@/libnest

liblenstool_a_SOURCES =  \
	al_sq_point.c      amplif.c            amplif_mat.c    amplif_matinv.c \
	bayesapp.c         bayesys3.c          bicubic_coef.c  bicubic_int.c \
	brent.c             chi_invim.c     chsigne.c \
	classer.c          cleanlens.c         comp_chi_osv.c \
	comp_dchi_osv.c    copyright.c         cor_seeing.c    cor_shear.c \
	cosmoratio.c       cp_errdiff.c        cp_images.c     cp_im.c \
	crea_filtre.c      critic_an.c         criticinv.c     criticnew.c \
	cv_cpsf.c          d_binning.c         d_bruiter.c     d_gauss.c \
	diag.c             diff_mag.c          d_integrer.c    dist.c     \
	dist_min.c         distor.c            do_itos.c       d_poisson.c \
	d_profil.c         d_rndschechter.c    d_rndtype.c     d_rndz.c   d_rndzsmail.c \
	d_seeing.c         e_amp.c             ecrire_r.c      e_dpl.c \
        dyn_mass.c         e_cube.c                                     \
	e_g2cpx.c          e_gcpx.c            e_giant.c       e_grad2.c \
	e_grad.c           e_hernquist.c e_im_prec.c         e_inthere.c     e_lens.c \
	e_lensing.c        e_lens_P.c          e_mag.c         e_mass.c \
	e_nfw.c            e_nfwgt.c           e_nfwtri.c      e_pcpx.c        e_pixel.c \
	e_pot.c            err_invim.c         e_sersic.c      e_einasto.c	e_tau.c \
	e_test.c           e_testg.c           e_test_P.c      e_time.c \
	e_unlens.c         e_unlens_fast.c     e_unlensgrid.c  e_unmag.c \
	e_zeroamp.c        fft.c               fftcc_im.c      fftc_im.c \
	fmin_ell.c         follow.c            followi.c       formeli.c \
	frprmn.c           fr_sq_point.c       f_shape2.c      f_shape3.c \
	f_shape4.c         f_shape_abs.c       f_shape.c       f_shmap.c \
	f_source.c         fz_dlsds.c          g_ampli.c       g_amplif.c \
	g_curv.c           g_dpl.c             g_grid.c        g_mass.c \
	g_poten.c          g_profil.c          g_prop.c        g_radial.c \
	grid.c             gridp.c             gs_ecrire.c     gs_ecrire_t.c \
	g_shear.c          g_shearf.c          g_time.c        hilbert.c \
	ic_product.c       im_add_const.c      i_marker.c      im_bin.c \
	im_convolve.c      im_stat.c           imtosou.c       inconvexe.c \
	init_grille.c      inside.c            interpol.c      inverse.c \
	isoima.c           iter.c              keep_cl.c       lire.c \
	local.c            lt_math.c           min_parabol.c \
	min_slo_par.c      multiscale_grid.c   nimage.c        norm_angle.c \
	o_big_slope.c      o_chi.c             o_chi_flux.c    \
	o_chires.c         o_dpl.c             o_flux.c        o_global.c \
	o_invim.c          o_keep_min.c        o_keepz_min.c   o_kp_sp_min.c \
	o_mag.c            o_mag_m.c           o_min_loc.c     o_min_slope.c \
	o_prep.c           o_prep_mult.c       o_print.c       o_print_res.c \
	opt_source.c       o_rescale.c         o_run1.c        o_run2.c \
	o_run_bayes.c      o_run.c             o_run_mc0.c     o_run_mc.c \
	o_runpot1.c        o_runpot2.c         o_runz0.c       o_runz2.c \
	o_runz.c           o_runz_mc0.c        o_runz_mc.c     o_set_err.c \
	o_set_exc.c        o_set_lens_bayes.c  o_set_lens.c    o_set_lmax.c \
	o_set_lmin.c       o_set_map.c         o_set_map_mc.c  o_set_map_z.c o_set_source.c \
	o_set_start.c      o_shape.c           o_slope_sp.c    o_stat.c \
	o_step.c           o_swi_big.c         polxy.c         pro_arclet.c \
	ptau.c              random.c        random_lt.c \
	r_cleanlens.c      r_cline.c           r_cosmolimit.c  r_cosmologie.c \
	readBayesModels.c  readimage.c     readcube.c read_lenstable.c \
	r_frame.c          r_grille.c          r_image.c       r_large.c \
	r_limit.c          r_msgrid.c          r_observ.c      rotation.c \
	rotmatrix.c        r_potentiel.c       r_potfile.c     r_vfieldlimit.c \
	r_runmode.c        r_source.c          r_shapelimit.c r_shapemodel.c \
	r_vfield.c         r_dynfile.c        read_table_einasto.c 	s_compmag.c         set_default.c \
	set_lens.c         set_lens_par.c      set_potfile.c \
	set_res_par.c      sig2posS4j.c        sig2posS.c      sig2posSe.c \
	sig2posSj.c        signe.c             slope.c         sort.c \
	sortf.c            spline.c            split_image.c   sp_main.c \
	sp_set.c           s_sof.c             s_sourcebox.c   s_source.c \
	st_opt.c           study.c             tirer.c         tracepot.c \
	unlens1.c          unlens_bc.c         updatecut.c     update_emass.c \
	update_epot.c      w_critic.c          weight_baryc.c  w_gianti.c \
	w_prop.c           wr_arclet.c         wr_mass.c       wr_pot.c \
	w_sicat.c          zero.c              zero_t.c        zonemult.c \
	omp_funs.c d_seeing_omp.c d_bruiter_omp.c#\ o_run_nest.c

bin_PROGRAMS 	= lenstool

lenstool_SOURCES = main.c

lenstool_LDFLAGS = @CFITSIO_LIB@ @WCS_LIB@ @GSL_LIB@

lenstool_LDADD 	= liblenstool.a ../liblt/liblt.a -lcfitsio -lwcs -lgsl -lgslcblas -lm  
#lenstool_LDADD 	= liblenstool.a ../liblt/liblt.a ../libnest/libnest3.a -L/usr/local/lib -lgfortran -L/Users/ejullo/softs/lapack-3.2.1 -llapack -lblas -lcfitsio -lwcs -lm
