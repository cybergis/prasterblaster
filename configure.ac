AC_INIT([prasterblaster], [0.01], [dmattli@usgs.gov])  
AM_INIT_AUTOMAKE([-Wall -Werror foreign])  
AC_PROG_CC([mpiCC mpicc])
AC_PROG_CXX([mpic++ mpiCXX mpicxx])

AC_PROG_RANLIB

AC_CHECK_LIB(mpi, MPI_Init, [], [echo "MPI not found!"; exit 1])
AC_CHECK_LIB(gdal, GDALAllRegister, [], [echo "GDAL not found!"; exit])

AC_CONFIG_HEADERS([config.h])  
AC_CONFIG_FILES([  
 Makefile  
 src/Makefile
 src/pRPL/Makefile
 src/gctp_cpp/Makefile
])  
AC_OUTPUT  
