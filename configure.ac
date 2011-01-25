# Filename: configure.ac
# Author: David Mattli <dmattli@usgs.gov
# License: PUBLIC DOMAIN


AC_INIT([prasterblaster], [0.01], [dmattli@usgs.gov])  
#m4_include([m4/ax_boost_python.m4])
#m4_include([m4/ax_python.m4])
	 	
AM_INIT_AUTOMAKE([-Wall -Werror foreign])  
AC_PROG_LIBTOOL
AC_CONFIG_MACRO_DIR([m4])
	 	
AC_ARG_WITH(gdal-incdir,
  AS_HELP_STRING( [--with-gdal-incdir=path],
  [Specify the path to the GDAL includes]),
  [GDAL_CPPFLAGS="-I$withval"], [])

AC_ARG_WITH(gdal,
  AS_HELP_STRING( [--with-gdal-libdir=path],
  [Specify the path to the GDAL library]),
  [GDAL_CFLAGS="-L $withval"], [])

AC_ARG_WITH(mpi-c-wrapper,  
  AS_HELP_STRING( [--with-mpi-c-wrapper=path],
  [Specify the path to the mpi c compiler wrapper]),
  [CC=$withval], [])
AC_ARG_WITH(mpi-cxx-wrapper,
  AS_HELP_STRING( [--with-mpi-cxx-wrapper=path],
  [Specify the path to the mpi c++ compiler wrapper]),
  [CXX=$withval], [])

AC_SUBST(GDAL_CPPFLAGS)
AC_SUBST(GDAL_CFLAGS)
CPPFLAGS="${CPPFLAGS} ${GDAL_CPPFLAGS}"
CFLAGS="$CFLAGS $GDAL_CFLAGS"

AC_PROG_CC([mpicc])
AC_PROG_CXX([mpiCC mpiCXX mpic++])
AC_PROG_CPP

# Check for GDAL library and headers
AC_CHECK_HEADER(gdal.h, [], [echo 'gdal.h not found! Set --with-gdal-incdir=path'; exit])
AC_CHECK_LIB(gdal, GDALAllRegister, [], 
             [echo "GDAL library not found! Set --with-gdal-libdir=path"; exit])
	 	
AC_CONFIG_HEADERS([config.h])  
AC_CONFIG_FILES([  
  Makefile  
  src/Makefile
  src/gctp_cpp/Makefile
])  
	 	
AC_OUTPUT  