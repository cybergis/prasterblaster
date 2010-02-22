AC_INIT([prasterblaster], [0.01], [dmattli@usgs.gov])  
AM_INIT_AUTOMAKE([-Wall -Werror foreign])  
AC_PROG_CC([mpiCC mpicc])
AC_PROG_CXX([mpic++ mpiCXX mpicxx])

AC_CHECK_LIB(mpi, MPI_Init)
AC_CHECK_LIB(gctp, adjust_lon)

AC_CONFIG_HEADERS([config.h])  
AC_CONFIG_FILES([  
 Makefile  
 src/Makefile  
])  
AC_OUTPUT  