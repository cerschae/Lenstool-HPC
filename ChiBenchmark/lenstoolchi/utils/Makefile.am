bin_PROGRAMS = bayesChires bayesImage bayesSource bayesBest bayesDlsDs bayesGrad bayesMap bayesAmp 
if PGPLOT
bin_PROGRAMS += Histogram Histogram2D 
endif

AM_FFLAGS =  -fPIC -w -O3 -fno-automatic -finit-local-zero --fixed-line-length-none 
AM_CFLAGS = -I@top_srcdir@/include -I@top_srcdir@/liblt @GSL_HDR@

bayesChires_SOURCES = bayesChires.c
bayesImage_SOURCES = bayesImage.c
bayesSource_SOURCES = bayesSource.c
bayesBest_SOURCES = bayesBest.c
bayesDlsDs_SOURCES = bayesDlsDs.c
bayesGrad_SOURCES = bayesGrad.c
bayesAmp_SOURCES = bayesAmp.c
bayesMap_SOURCES = bayesMap.c


bayesChires_LDFLAGS = @CFITSIO_LIB@ @WCS_LIB@ @GSL_LIB@
bayesChires_LDADD = @top_builddir@/src/liblenstool.a \
		    @top_builddir@/liblt/liblt.a -lcfitsio -lwcs -lm -lgsl -lgslcblas      

bayesImage_LDFLAGS = @CFITSIO_LIB@ @WCS_LIB@ @GSL_LIB@
bayesImage_LDADD = @top_builddir@/src/liblenstool.a \
		    @top_builddir@/liblt/liblt.a -lcfitsio -lwcs -lm -lgsl -lgslcblas      

bayesSource_LDFLAGS = @CFITSIO_LIB@ @WCS_LIB@ @GSL_LIB@
bayesSource_LDADD = @top_builddir@/src/liblenstool.a \
		    @top_builddir@/liblt/liblt.a -lcfitsio -lwcs -lm -lgsl -lgslcblas      

bayesBest_LDFLAGS = @CFITSIO_LIB@ @WCS_LIB@ @GSL_LIB@
bayesBest_LDADD = @top_builddir@/src/liblenstool.a \
		    @top_builddir@/liblt/liblt.a -lcfitsio -lwcs -lm -lgsl -lgslcblas      

bayesDlsDs_LDFLAGS = @CFITSIO_LIB@ @WCS_LIB@ @GSL_LIB@
bayesDlsDs_LDADD = @top_builddir@/src/liblenstool.a \
		    @top_builddir@/liblt/liblt.a -lcfitsio -lwcs -lm -lgsl -lgslcblas      

bayesGrad_LDFLAGS = @CFITSIO_LIB@ @WCS_LIB@ @GSL_LIB@
bayesGrad_LDADD = @top_builddir@/src/liblenstool.a \
		    @top_builddir@/liblt/liblt.a -lcfitsio -lwcs -lm -lgsl -lgslcblas      

bayesMap_LDFLAGS = @CFITSIO_LIB@ @WCS_LIB@ @GSL_LIB@
bayesMap_LDADD = @top_builddir@/src/liblenstool.a \
		    @top_builddir@/liblt/liblt.a -lcfitsio -lwcs -lm -lgsl -lgslcblas      

bayesAmp_LDFLAGS = @CFITSIO_LIB@ @WCS_LIB@ @GSL_LIB@
bayesAmp_LDADD = @top_builddir@/src/liblenstool.a \
		    @top_builddir@/liblt/liblt.a -lcfitsio -lwcs -lm -lgsl -lgslcblas      

if PGPLOT
Histogram_SOURCES = histogram.f 

Histogram2D_SOURCES = histogram2D.f

Histogram_LDFLAGS =  @CFITSIO_LIB@ @PGPLOT_LIB@ @X_LIBS@ @X_PRE_LIBS@ @X_EXTRA_LIBS@
Histogram_LDADD = -lpgplot -lX11 -lpng -lcfitsio
Histogram2D_LDFLAGS =  @CFITSIO_LIB@ @PGPLOT_LIB@ @X_LIBS@ @X_PRE_LIBS@ @X_EXTRA_LIBS@
Histogram2D_LDADD = -lpgplot -lX11 -lpng -lcfitsio
endif

EXTRA_DIST = bayesChires.h histogram.inc \
	     sort.f smooth1D.f regrid.f  regrid1D.f   rfft.f    smooth.f

