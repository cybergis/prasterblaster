AC_INIT([prasterblaster], [0.01], [dmattli@usgs.gov])  
m4_include([m4/ax_boost_python.m4])
m4_include([m4/ax_python.m4])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])  
AC_PROG_CC([mpicc])
AC_PROG_CXX([mpiCC mpiCXX mpic++])


AC_PROG_LIBTOOL
AM_PATH_CHECK()

#AC_CHECK_LIB(mpi, MPI_Init, [], [echo "MPI not found!"; exit 1])
AX_BOOST_PYTHON


AC_ARG_ENABLE(mpi-c-wrapper,
	AS_HELP_STRING( [--with-mpi-c-wrapper=path],
	[Specify the path to the mpi c compiler wrapper]),
	[CC=$withval], [])
AC_ARG_ENABLE(mpi-cxx-compiler,
	AS_HELP_STRING( [--with-cxx-mpi-wrapper=path],
	[Specify the path to the mpi c++ compiler wrapper]),
	[CXX=$withval], [])

AC_CHECK_LIB(gdal, GDALAllRegister, [], [echo "GDAL not found!"; exit])

AC_CONFIG_HEADERS([config.h])  
AC_CONFIG_FILES([  
 Makefile  
 src/Makefile
 src/gctp_cpp/Makefile
 src/python/Makefile
 tests/Makefile
])  
AC_OUTPUT  
