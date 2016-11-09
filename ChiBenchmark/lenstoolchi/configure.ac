# configure.ac file for Lenstool
# Eric Jullo (Jan 2014)
AC_INIT(Lenstool, 6.8, [Jean-Paul Kneib <jean-paul.kneib@epfl.ch>, Eric Jullo <eric.jullo@lam.fr>])
AC_CONFIG_SRCDIR([src/e_pot.c],[liblt/rdf_fits.c])
AC_CONFIG_AUX_DIR(autoconf)
AM_INIT_AUTOMAKE([no-define])

AC_PROG_CC([icc gcc cc cl])
AC_PROG_RANLIB

# Check debug option
AC_ARG_WITH([debug], [AS_HELP_STRING([--with-debug],
        [Compile in debug mode])],[],[])

if test "x$with_debug" = xyes ; then
	CFLAGS=$(echo "$CFLAGS -g" | sed -e "s/-O./-O0/")
    AC_MSG_NOTICE([Enable debug mode])
fi	

# Check OpenMP support
AC_OPENMP
CFLAGS="$CFLAGS $OPENMP_CFLAGS"

# Checks for m library
AC_CHECK_LIB([m],[sin])

# Check CFITSIO library
AC_ARG_WITH([cfitsio-include-path], [AS_HELP_STRING([--with-cfitsio-include-path=DIR],
            [directory containing cfitsio directory, defaults to /usr/include])], [])

if test "x$with_cfitsio_include_path" = x ; then
    with_cfitsio_include_path="/usr/include"
fi
AC_CHECK_FILE([${with_cfitsio_include_path}/fitsio.h], [], 
              AC_MSG_ERROR([Include file fitsio.h not found]))
CFITSIO_HDR="-I$with_cfitsio_include_path"

AC_ARG_WITH([cfitsio-lib-path], [AS_HELP_STRING([--with-cfitsio-lib-path],
            [path to the  cfitsio libraries, defaults to /usr/lib])], [])
if test "x$with_cfitsio_lib_path" = x ; then
    with_cfitsio_lib_path="/usr/lib"
fi
CFITSIO_LIB="-L$with_cfitsio_lib_path"
LIBS="$CFITSIO_LIB -lm"
AC_CHECK_LIB([cfitsio], [ffdkopn], [], 
             AC_MSG_ERROR([Function ffdkopn in cfitsio library not found]))

AC_SUBST(CFITSIO_HDR)
AC_SUBST(CFITSIO_LIB)

# Check WCSLIB library
AC_ARG_WITH([wcslib-include-path], [AS_HELP_STRING([--with-wcslib-include-path=DIR],
            [directory containing wcslib directory, defaults to /usr/include])], [])

if test "x$with_wcslib_include_path" = x ; then
    with_wcslib_include_path="/usr/include"
fi
AC_CHECK_FILE([${with_wcslib_include_path}/wcs.h], [], 
               AC_MSG_ERROR([Include file wcs.h not found]))
WCS_HDR="-I$with_wcslib_include_path"

AC_ARG_WITH([wcslib-lib-path], [AS_HELP_STRING([--with-wcslib-lib-path],
            [path to the  wcslib libraries, defaults to /usr/lib])], [])
if test "x$with_wcslib_lib_path" = x ; then
    with_wcslib_lib_path="/usr/lib"
fi
WCS_LIB="-L$with_wcslib_lib_path"
LIBS="$WCS_LIB -lm"
AC_CHECK_LIB([wcs], [pix2wcs], [], 
             AC_MSG_ERROR([Function pix2wcs in wcs library not found]))

AC_SUBST(WCS_HDR)
AC_SUBST(WCS_LIB)

# Check PGPLOT library
AC_ARG_WITH([pgplot], [AS_HELP_STRING([--with-pgplot],
            [Enable PGPLOT tools])], [])

if test "x$with_pgplot" = x; then
    AC_MSG_WARN([PGPLOT disabled])
else

    AC_PROG_F77
    AC_PATH_XTRA

    AC_ARG_WITH([pgplot-include-path], [AS_HELP_STRING([--with-pgplot-include-path=DIR],
                [directory containing pgplot header files])], [])
    if test "x$with_pgplot_include_path" = x ; then
        with_pgplot_include_path="/usr/include"
    fi
    AC_CHECK_FILE([${with_pgplot_include_path}/cpgplot.h], [],
                  AC_MSG_ERROR([File cpgplot.h not found]))
    PGPLOT_HDR="-I$with_pgplotinc"

    AC_ARG_WITH([pgplot-lib-path], [AS_HELP_STRING([--with-pgplot-lib-path=DIR],
              [directory containing pgplot library, defaults to /usr/lib])], [])

    if test "x$with_pgplot_lib_path" = x ; then
        with_pgplot_lib_path="/usr/lib"   
    fi
    PGPLOT_LIB="-L$with_pgplot_lib_path"

    AC_MSG_NOTICE([PGPLOT enabled])
fi 
AM_CONDITIONAL([PGPLOT], [test x$with_pgplot = xyes])

AC_SUBST(PGPLOT_HDR)
AC_SUBST(PGPLOT_LIB)

# Check GSL library
AC_ARG_WITH([gsl-include-path], [AS_HELP_STRING([--with-gsl-include-path=DIR],
            [path to the gsl directory, defaults to /usr/include])], [])

if test "x$with_gsl_include_path" = x ; then
    with_gsl_include_path="/usr/include"
fi

AC_CHECK_FILE([${with_gsl_include_path}/gsl/gsl_cblas.h], [],
              AC_MSG_ERROR([Invalid path to gsl/gsl_cblas.h]))
GSL_HDR="-I$with_gsl_include_path"

AC_ARG_WITH([gsl-lib-path], [AS_HELP_STRING([--with-gsl-lib-path],
            [path to the  gsl libraries, defaults to /usr/lib])], [])
if test "x$with_gsl_lib_path" = x ; then
    with_gsl_lib_path="/usr/lib"
fi
GSL_LIB="-L$with_gsl_lib_path"
LIBS="$GSL_LIB"
AC_CHECK_LIB([gslcblas],[cblas_dgemm], [],
              AC_MSG_ERROR([Function cblas_dgemm in gslcblas library not found]))
AC_CHECK_LIB([gsl],[gsl_blas_dgemm], [],
              AC_MSG_ERROR([Function gsl_blas_dgemm in gsl library not found]))

AC_SUBST(GSL_HDR)
AC_SUBST(GSL_LIB)
LIBS=

AC_CONFIG_FILES([Makefile src/Makefile liblt/Makefile table_src/Makefile utils/Makefile])
AC_OUTPUT
