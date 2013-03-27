# Filename: configure.ac
# Author: David Mattli <dmattli@usgs.gov
# License: PUBLIC DOMAIN
# Portions from to

AC_INIT([prasterblaster], [0.01], [dmattli@usgs.gov])  
	 	
AM_INIT_AUTOMAKE([-Wall -Werror foreign])  
AC_PROG_LIBTOOL
AC_CONFIG_MACRO_DIR([m4])

LT_INIT([shared static])

# Check for programs
AC_PROG_CXX([mpiCC mpiCXX mpic++])
AC_PROG_CPP
AC_PROG_INSTALL

AC_LANG([C++])
	 	
AC_ARG_WITH(gdal-incdir,
  AS_HELP_STRING( [--with-gdal-incdir=path],
  [Specify the path to the GDAL includes]),
  [GDAL_CPPFLAGS="-I$withval"], [])

AC_ARG_WITH(gdal-libdir,
  AS_HELP_STRING( [--with-gdal-libdir=path],
  [Specify the path to the GDAL library]),
  [GDAL_CXXFLAGS="-L$withval -Wl,-R$withval"], [])

AC_ARG_WITH(tiff-incdir,
  AS_HELP_STRING( [--with-tiff-incdir=path],
  [Specify the path to the tiff includes]),
  [TIFF_CPPFLAGS="-I$withval"], [])

AC_ARG_WITH(tiff-libdir,
  AS_HELP_STRING( [--with-tiff-libdir=path],
  [Specify the path to the tiff shared library file]),
  [TIFF_CXXFLAGS="-L$withval -Wl,-R$withval"], [])
AC_SUBST(GDAL_CPPFLAGS)
AC_SUBST(GDAL_CFLAGS)
CPPFLAGS="${CPPFLAGS} ${GDAL_CPPFLAGS} ${TIFF_CPPFLAGS}"
CXXFLAGS="${CXXFLAGS} ${GDAL_CXXFLAGS} ${TIFF_CXXFLAGS}"

LIBS="$LIBS -lm"

dnl Check for sincos
AC_MSG_CHECKING([for sincos])
AC_LINK_IFELSE([
AC_LANG_SOURCE([[
#include <stdio.h>
#include <math.h>
int main (void) {
double s, c;
scanf ("%f", &s);
sincos (s, &s, &c);
printf ("%f %f", s, c);
return 0;
}
]])
],[
  AC_MSG_RESULT([yes])
  AC_DEFINE(HAVE_SINCOS, 1, [Define if the sincos() function is available])
],[
  AC_MSG_RESULT([no])
])


dnl Check for shared_ptr
AC_MSG_CHECKING([for std::shared_ptr])
AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([[#include <memory>]],
                     [[std::shared_ptr<int> test_shared_ptr;]])],
    [AC_MSG_RESULT([yes])
     AC_DEFINE_UNQUOTED([HAVE_STD_SHARED_PTR], 1, [Define to 1 if you have std::shared_ptr.])],
    [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <tr1/memory>]],
                                        [[std::tr1::shared_ptr<int> test_shared_ptr;]])],
         [AC_MSG_RESULT([yes])
          AC_DEFINE_UNQUOTED([HAVE_TR1_SHARED_PTR], 1, [Define to 1 if you have std::tr1::shared_ptr.])],
         [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <boost/shared_ptr.hpp>]],
                                             [[boost::shared_ptr<int> test_shared_ptr;]])],
            [AC_MSG_RESULT([yes])
             AC_DEFINE_UNQUOTED([HAVE_BOOST_SHARED_PTR], 1, [Define to 1 if you have boost::shared_ptr.])],
            [AC_MSG_RESULT([no])
             AC_MSG_ERROR([configure cannot find shared_ptr.])])])])


# Check for GDAL library and headers
AC_CHECK_HEADER([gdal.h], [], [echo 'gdal.h not found! Set --with-gdal-incdir=path'; exit])
AC_SEARCH_LIBS([GDALAllRegister], [gdal], [], [AC_MSG_ERROR(['libgdal not found! Set --with-gdal-libdir=path'])], [])

# Check for tiff headers
AC_CHECK_HEADER([tiff.h], [], [echo 'tiff.h not found! Set --with-tiff-incdir=path'; exit])
AC_SEARCH_LIBS([TIFFOpen], [tiff], [], [AC_MSG_ERROR(['libtiff not found! Set --with-tiff-libdir=path'])], [])
AC_CONFIG_HEADERS([config.h]) 
AC_CONFIG_FILES([  
  Makefile  
  src/Makefile
  src/gctp_cpp/Makefile
  src/demos/Makefile
  src/gtest/Makefile
  tests/Makefile
])  
	 	
AC_OUTPUT  
