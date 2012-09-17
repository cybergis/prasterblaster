

noinst_LTLIBRARIES = libsptw.la
libsptw_la_SOURCES = sptw.cc
libsptw_la_CPPFLAGS = -I../../src/

bin_PROGRAMS = top/prasterblasterpio
top_prasterblasterpio_SOURCES = prasterblaster-pio.cc
top_prasterblasterpio_CPPFLAGS = -I..
top_prasterblasterpio_LDFLAGS = $(all_libraries)  -L../gctp_cpp -L..
top_prasterblasterpio_LDADD = ../librasterblaster.la ../gctp_cpp/libgctp_cpp.la libsptw.la
top_prasterblasterpio_CXXFLAGS = -lgdal -lrasterblaster 


AM_CXXFLAGS = -W -Wall -Wextra -Wcast-align \
	-Wpointer-arith -Wsign-compare -g\
	-Wformat=2 -Wno-format-y2k -Wmissing-braces -Wparentheses \
	-Wtrigraphs -Wstrict-aliasing=2
