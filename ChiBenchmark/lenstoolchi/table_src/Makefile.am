# Program Makefile for lenstool
bin_PROGRAMS 	= lenstool_tab
AM_CPPFLAGS 	= -I@top_builddir@include @CFITSIO_HDR@

lenstool_tab_SOURCES 	= 	\
		e_nfwg.c lenstool_tab.c midpnt.c nrutil.c polint.c qromo.c \
		read_bin.c lens1jnew_func.h  lenstool_tab_func.h  nrutil.h \
		vd2dlens_func.h

lenstool_tab_LDFLAGS	= @CFITSIO_LIB@ -lcfitsio -lm

